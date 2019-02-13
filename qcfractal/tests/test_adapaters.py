"""
Explicit tests for queue manipulation.
"""

from numpy import allclose
import pytest

import qcfractal.interface as portal
from qcfractal import testing
from qcfractal.testing import (
    reset_server_database, queue_manager, managed_compute_server, parameterizable_fractal_compute_server)
import json


@testing.using_rdkit
def test_adapter_single(managed_compute_server):
    client, server, manager = managed_compute_server
    reset_server_database(server)

    # Add compute
    hooh = portal.data.get_molecule("hooh.json")
    ret = client.add_compute("rdkit", "UFF", "", "energy", None, [hooh.json_dict()], tag="other")

    # Force manager compute and get results
    manager.await_results()
    ret = client.get_results()
    assert len(ret) == 1


@pytest.mark.parametrize("cores_per_task,memory_per_task", [
    (None, None),
    (1, 1),
    (2, 1.9)
])
@testing.using_psi4
def test_keyword_args_passing(parameterizable_fractal_compute_server, cores_per_task, memory_per_task):
    psi4_mem_buffer = 0.95  # Memory consumption buffer on psi4
    client, server, adapter_client = parameterizable_fractal_compute_server
    reset_server_database(server)
    tasks = [  # Emulate the QueueManager test function
        json.loads(json.dumps({
            "id": "123456789012345678901234",  # Placeholder ID
            "spec": {
                "function":
                    "qcengine.compute",
                "args": [{
                    "molecule": portal.data.get_molecule("hooh.json").json(as_dict=True),
                    "driver": "energy",
                    "model": {"method": "HF",
                              "basis": "sto-3g"},
                    "keywords": {},
                    "return_output": True,
                    'qcfractal_tags': {'program': 'psi4', 'options': None}
                }, "psi4"],
                "kwargs": {}
            },
            "parser": "single",
            "hooks": [],
            "tag": "other"
        }))
    ]
    with queue_manager(client,
                       adapter_client,
                       cores_per_task=cores_per_task,
                       memory_per_task=memory_per_task) as manager:
        manager.queue_adapter.submit_tasks(tasks)
        manager.await_results()
        ret = client.get_results()
        assert len(ret) == 1
        provenance = ret[0]['provenance']
        if cores_per_task is not None:
            assert provenance['nthreads'] == cores_per_task
        if memory_per_task is not None:
            assert allclose(provenance['memory'], memory_per_task*psi4_mem_buffer)


@testing.using_rdkit
def test_adapter_error_message(managed_compute_server):
    client, server, manager = managed_compute_server
    reset_server_database(server)

    # HOOH without connectivity, RDKit should fail
    hooh = portal.data.get_molecule("hooh.json").json_dict()
    del hooh["connectivity"]
    mol_ret = client.add_molecules({"hooh": hooh})

    ret = client.add_compute("rdkit", "UFF", "", "energy", None, mol_ret["hooh"])
    queue_id = ret.submitted[0]

    # Nothing should have happened yet
    assert len(manager.list_current_tasks()) == 0

    # Pull out a special iteration on the queue manager
    manager.update()
    assert len(manager.list_current_tasks()) == 1

    manager.await_results()
    assert len(manager.list_current_tasks()) == 0

    db = server.objects["storage_socket"]
    ret = db.get_queue(status="ERROR")["data"]

    assert len(ret) == 1
    assert "connectivity graph" in ret[0]["error"]["error_message"]
    server.objects["storage_socket"].queue_mark_complete([queue_id])


@testing.using_rdkit
def test_adapter_raised_error(managed_compute_server):
    client, server, manager = managed_compute_server
    reset_server_database(server)

    # HOOH without connectivity, RDKit should fail
    hooh = portal.data.get_molecule("hooh.json").json_dict()
    del hooh["connectivity"]
    mol_ret = client.add_molecules({"hooh": hooh})

    ret = client.add_compute("something_bad", "UFF", "", "energy", None, mol_ret["hooh"])
    queue_id = ret.submitted[0]

    manager.await_results()

    db = server.objects["storage_socket"]
    ret = db.get_queue(status="ERROR")["data"]

    assert len(ret) == 1
    assert "QCEngine Call Error" in ret[0]["error"]["error_message"]
    server.objects["storage_socket"].queue_mark_complete([queue_id])
