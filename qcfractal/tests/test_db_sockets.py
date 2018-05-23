"""
Tests the database wrappers

All tests should be atomic, that is create and cleanup their data
"""

import pytest
import numpy as np

# Import the DQM collection
import qcfractal as dserver
import qcfractal.interface as qp


@pytest.fixture(scope="module", params=["mongo"])
def db_socket(request):
    print("")
    db_name = "dqm_local_values_test"

    # IP/port/drop table is specific to build
    if request.param == "mongo":
        db = dserver.db_socket_factory("127.0.0.1", 27017, db_name, db_type=request.param)

        # Clean and re-init the databse
        db.client.drop_database(db._project_name)
        db.init_database()
    else:
        raise KeyError("DB type %s not understood" % request.param)

    yield db

    if request.param == "mongo":
        db.client.drop_database(db_name)
    else:
        raise KeyError("DB type %s not understood" % request.param)



def test_molecules_add(db_socket):

    water = qp.data.get_molecule("water_dimer_minima.psimol")

    # Add once
    ret1 = db_socket.add_molecules({"new_water": water.to_json()})
    assert ret1["meta"]["success"] is True
    assert ret1["meta"]["n_inserted"] == 1

    # Try duplicate adds
    ret2 = db_socket.add_molecules({"new_water2": water.to_json()})
    assert ret1["meta"]["success"] is True
    assert ret2["meta"]["n_inserted"] == 0
    assert ret2["meta"]["duplicates"][0] == "new_water2"

    # Assert the ids match
    assert ret1["data"]["new_water"] == ret2["data"]["new_water2"]

    # Pull molecule from the DB for tests
    db_json = db_socket.get_molecules(water.get_hash(), index="hash")[0]
    water.compare(db_json)

    # Cleanup adds
    ret = db_socket.del_molecules(water.get_hash(), index="hash")
    assert ret == 1


def test_molecules_add_many(db_socket):
    water = qp.data.get_molecule("water_dimer_minima.psimol")
    water2 = qp.data.get_molecule("water_dimer_stretch.psimol")

    ret = db_socket.add_molecules({"water1": water.to_json(), "water2": water2.to_json()})
    assert ret["meta"]["n_inserted"] == 2

    # Cleanup adds
    ret = db_socket.del_molecules([water.get_hash(), water2.get_hash()], index="hash")
    assert ret == 2

    ret = db_socket.add_molecules({"water1": water.to_json(), "water2": water2.to_json()})
    assert ret["meta"]["n_inserted"] == 2

    # Cleanup adds
    ret = db_socket.del_molecules(list(ret["data"].values()), index="id")
    assert ret == 2


def test_molecules_get(db_socket):

    water = qp.data.get_molecule("water_dimer_minima.psimol")

    # Add once
    ret = db_socket.add_molecules({"water": water.to_json()})
    assert ret["meta"]["n_inserted"] == 1
    water_id = ret["data"]["water"]

    # Pull molecule from the DB for tests
    db_json = db_socket.get_molecules(water_id, index="id")[0]
    water_db = qp.Molecule.from_json(db_json)
    water_db.compare(water)

    # Cleanup adds
    ret = db_socket.del_molecules(water_id, index="id")
    assert ret == 1


def test_options_add(db_socket):

    opts = qp.data.get_options("psi_default")

    ret = db_socket.add_options([opts, opts])
    assert ret["n_inserted"] == 1

    ret = db_socket.add_options(opts)
    assert ret["n_inserted"] == 0

    del opts["_id"]
    assert opts == db_socket.get_options({"name": opts["name"], "program": opts["program"]})[0]

    assert 1 == db_socket.del_option(opts["program"], opts["name"])

def test_options_error(db_socket):
    opts = qp.data.get_options("psi_default")

    del opts["name"]
    ret = db_socket.add_options(opts)
    assert ret["n_inserted"] == 0
    assert len(ret["validation_errors"]) == 1

def test_databases_add(db_socket):

    db = {"category": "OpenFF", "name": "Torsion123", "something": "else", "array": ["54321"]}

    ret = db_socket.add_database(db)
    del db["_id"]
    assert ret["n_inserted"] == 1

    new_db = db_socket.get_database(db["category"], db["name"])
    assert db == new_db

    ret = db_socket.del_database(db["category"], db["name"])
    assert ret == 1


def test_results_add(db_socket):

    # Add two waters
    water = qp.data.get_molecule("water_dimer_minima.psimol")
    water2 = qp.data.get_molecule("water_dimer_stretch.psimol")
    mol_insert = db_socket.add_molecules({"water1": water.to_json(), "water2": water2.to_json()})

    page1 = {
        "molecule_id": mol_insert["data"]["water1"],
        "method": "M1",
        "basis": "B1",
        "option": "default",
        "program": "P1",
        "other_data": 5
    }

    page2 = {
        "molecule_id": mol_insert["data"]["water2"],
        "method": "M1",
        "basis": "B1",
        "option": "default",
        "program": "P1",
        "other_data": 10
    }

    ret = db_socket.add_results([page1, page2])
    assert ret["n_inserted"] == 2

    ret = db_socket.del_results(ret["ids"], index="id")
    assert ret == 2

    ret = db_socket.del_molecules(list(mol_insert["data"].values()), index="id")
    assert ret == 2

### Build out a set of query tests

@pytest.fixture(scope="module")
def db_results(db_socket):
    # Add two waters
    water = qp.data.get_molecule("water_dimer_minima.psimol")
    water2 = qp.data.get_molecule("water_dimer_stretch.psimol")
    mol_insert = db_socket.add_molecules({"water1": water.to_json(), "water2": water2.to_json()})

    page1 = {
        "molecule_id": mol_insert["data"]["water1"],
        "method": "M1",
        "basis": "B1",
        "option": "default",
        "program": "P1",
        "return_result": 5
    }

    page2 = {
        "molecule_id": mol_insert["data"]["water2"],
        "method": "M1",
        "basis": "B1",
        "option": "default",
        "program": "P1",
        "return_result": 10
    }

    page3 = {
        "molecule_id": mol_insert["data"]["water1"],
        "method": "M1",
        "basis": "B1",
        "option": "default",
        "program": "P2",
        "return_result": 15
    }

    page4 = {
        "molecule_id": mol_insert["data"]["water1"],
        "method": "M2",
        "basis": "B1",
        "option": "default",
        "program": "P2",
        "return_result": 15
    }

    page5 = {
        "molecule_id": mol_insert["data"]["water2"],
        "method": "M2",
        "basis": "B1",
        "option": "default",
        "program": "P1",
        "return_result": 20
    }

    pages_insert = db_socket.add_results([page1, page2, page3, page4, page5])

    yield db_socket

    # Cleanup
    ret = db_socket.del_results(list(pages_insert["data"].values()), index="id")
    assert ret == pages_insert["meta"]["n_inserted"]

    ret = db_socket.del_molecules(list(mol_insert["data"].values()), index="id")
    assert ret == mol_insert["meta"]["n_inserted"]


def test_results_query_total(db_results):

    assert 5 == len(db_results.get_results({}))


def test_results_query_method(db_results):

    assert 5 == len(db_results.get_results({"method": ["M2", "M1"]}))
    assert 2 == len(db_results.get_results({"method": ["M2"]}))
    assert 2 == len(db_results.get_results({"method": "M2"}))


def test_results_query_dual(db_results):

    assert 5 == len(db_results.get_results({"method": ["M2", "M1"], "program": ["P1", "P2"]}))
    assert 1 == len(db_results.get_results({"method": ["M2"], "program": "P2"}))
    assert 1 == len(db_results.get_results({"method": "M2", "program": "P2"}))


def test_results_query_project(db_results):
    tmp = db_results.get_results({"method": "M2", "program": "P2"}, projection={"return_result": True})[0]
    assert set(tmp.keys()) == {"return_result"}
    assert tmp["return_result"] == 15


