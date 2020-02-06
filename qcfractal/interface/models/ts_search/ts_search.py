"""
A module for performing reaction transition state searches

Todo:
  - Think whether GSM should have harnesses. probably yes, and it should be treated like we treat ESS.
  - Still need to add Gaussian harnesses to run ESS-driven optimizations
  - conda-package ARC's zmats module and import here (add dependency)
"""

from enum import Enum
import numpy as np

# from sqlalchemy.orm import column_property, relationship

from .rmg_db import determine_reaction_family
import qcelemental as qcel
# from qcfractal.storage_sockets.models.sql_models import KeywordsORM, KVStoreORM, MoleculeORM

from rmgpy.data.rmg import RMGDatabase
from rmgpy.reaction import Reaction
from rmgpy.species import Species


class TSSearch(object):
    def __init__(
            self,
            ts_search_id: int,
            well_ids_1: list,
            well_ids_2: list,
            methods: list,
            levels: dict,
            user_guesses: list = None,
            rmg_db: RMGDatabase = None
    ) -> None:
        """Initializes a TSSearch instance.

        Parameters
        ----------
        ts_search_id : int
            The TSSearch ID.
        well_ids_1 : list
            Entries are Molecule IDs, collectively describing a well (e,g,. reactants).
        well_ids_1 : list
            Entries are Molecule IDs, collectively describing a well (e,g,. products).
        methods : list
            The TS Search methods to carry out.
            Optional values: 'autotst', 'gsm', 'heuristics', 'kinbot', 'ml', 'neb_ase', 'neb_terachem', 'qst2', 'user'.
        levels: dict
            Keys are job types (allowed values are 'opt', 'freq', 'sp', 'irc'),
            values are the corresponding levels of theory.
            Note the IRC should be done at the same level as the geometry optimization,
            if a different level is specified an additional optimization will be spawned prior to the IRC calculation.
        user_guesses: list, optional
            Entries are string representations of Cartesian coordinate.
        rmg_db: RMGDatabase, optional
            The RMG database object, mandatory for the following methods: 'autotst', 'heuristics', 'kinbot'.
        """

        if any([not isinstance(well, list) for well in [well_ids_1, well_ids_1]]):
            raise TypeError(f'well_ids_1 and well_ids_2 must be lists, got {type(well_ids_1)}, {type(well_ids_2)}')
        if any([not isinstance(molecule_id, int) for molecule_id in well_ids_1 + well_ids_1]):
            raise TypeError(f'All molecule ID entries of a well must be integers, '
                            f'got {[type(mol_id) for mol_id in well_ids_1]} '
                            f'and {[type(mol_id) for mol_id in well_ids_1]}')
        if not isinstance(ts_search_id, int):
            raise TypeError(f'id should be an integer, got {ts_search_id} which is a {type(ts_search_id)}')

        self.ts_search_id = ts_search_id

        self.methods = [TSMethodsEnum(method) for method in methods]
        self.levels = {TSJobTypesEnum(key): val for key, val in levels.items()}

        self.well_ids_1 = well_ids_1
        self.well_ids_2 = well_ids_2
        self.well_1 = [get_molecule(id=molecule_id) for molecule_id in self.well_ids_1]  # Todo: I haven't figured out get_molecule() yet
        self.well_2 = [get_molecule(id=molecule_id) for molecule_id in self.well_ids_2]

        self.user_guesses = user_guesses
        self.rmg_db = rmg_db
        self.ts_guesses = dict()

        self.determine_rmg_reaction_family()

        self.search()

    def __repr__(self) -> str:
        """A short representation of the current TSSearch.

        Returns
        -------
        str
            The desired representation.
        """
        return f"TSSearch(ts_search_id='{self.ts_search_id}', " \
            f"well_ids_1='{self.well_ids_1}', well_ids_2='{self.well_ids_2}', " \
            f"methods='{self.methods}', levels={self.levels}, " \
            f"user_guesses={self.user_guesses}, rmgdb = {self.rmg_db})"

    def search(self) -> None:
        """
        Execute the selected TS search methods.
        Populates self.ts_guesses
        """

        for method in self.methods:
            if found(method):

                if method == TSMethodsEnum.autotst:
                    pass

                elif method == TSMethodsEnum.gsm:
                    pass

                elif method == TSMethodsEnum.heuristics:
                    pass

                elif method == TSMethodsEnum.kinbot:
                    pass

                elif method == TSMethodsEnum.ml:
                    pass

                elif method == TSMethodsEnum.neb_ase:
                    pass

                elif method == TSMethodsEnum.neb_terachem:
                    pass

                elif method == TSMethodsEnum.qst2:
                    pass

                elif method == TSMethodsEnum.user and self.user_guesses is not None:
                    for coords in self.user_guesses:
                        symbols, geometry = str_to_geometry(coords)  # Todo: How is this done in QCX?
                        if not len(qcel.molutil.guess_connectivity(symbols, geometry, threshold=0.9)):
                            # no colliding atoms, proceed
                            ts_guess_id = len(list(self.ts_guesses.keys()))  # 0-indexed
                            self.ts_guesses[ts_guess_id] = {'method': method,
                                                            'symbols': symbols,
                                                            'guess geometry': geometry,
                                                            'optimized geometry': None,
                                                            'energy': None,
                                                            'irc isomorphism': [None, None],
                                                            }

    def determine_rmg_reaction_family(self) -> None:
        """
        Create an RMG Reaction object corresponding to the given wells,
        and determine the RMG family.
        """

        if any(method in self.methods
               for method in [TSMethodsEnum.autotst, TSMethodsEnum.heuristics, TSMethodsEnum.kinbot]):
            rmg_reaction = Reaction(reactants=[Species(smiles=well.to_smiles()) for well in self.well_1],
                                    products=[Species(smiles=well.to_smiles()) for well in self.well_2])
            self.rmg_family = determine_reaction_family(self.rmg_db, rmg_reaction)
        else:
            self.rmg_family = None


class TSMethodsEnum(str, Enum):
    """
    The supported methods for a TS search. The methods which are available are a finite set.
    """

    autotst = 'autotst'  # AutoTST
    gsm = 'gsm'  # double ended growing string method (DE-GSM)
    heuristics = 'heuristics'  # brute force heuristics
    kinbot = 'kinbot'  # KinBot
    ml = 'ml'  # machine learning  Todo: will we have more than one ML module? probably yes... expand
    neb_ase = 'neb_ase'  # NEB in ASE
    neb_terachem = 'neb_terachem'  # NEB in TeraChem
    qst2 = 'qst2'  # Synchronous Transit-Guided Quasi-Newton (STQN) implemented in Gaussian
    user = 'user'  # user guesses


class TSJobTypesEnum(str, Enum):
    """
    The available job types in a TS search. The job types which are available are a finite set.
    """

    opt = 'opt'  # geometry optimization
    freq = 'freq'  # frequency calculation
    sp = 'sp'  # single point calculation
    irc = 'irc'  # internal redundant coordinate calculation


def found(method: TSMethodsEnum, raise_error: bool = False) -> bool:
    """
    Check whether a TSSearch method exists.

    Parameters
    ----------
    method: TSMethodsEnum
        The method to check.
    raise_error: bool, optional
        Whether to raise an error if the module was not found.

    Returns
    -------
    bool
        Whether the module was found,
    """

    # Hard coding heuristics and user guesses, since they always exist
    if method in [TSMethodsEnum.heuristics, TSMethodsEnum.user]:
        return True
    # Hard coding ML to return False, this repository is not up yet
    if method in [TSMethodsEnum.ml]:
        # The ML method is not implemented yet
        return False

    software_dict = {TSMethodsEnum.autotst: 'autotst',
                     TSMethodsEnum.gsm: 'gsm',
                     TSMethodsEnum.kinbot: 'kinbot',
                     TSMethodsEnum.neb_ase: 'ase',
                     TSMethodsEnum.neb_terachem: 'terachem',
                     TSMethodsEnum.qst2: 'gaussian',
                    }
    url_dict = {TSMethodsEnum.autotst: 'https://github.com/ReactionMechanismGenerator/AutoTST',
                TSMethodsEnum.gsm: 'https://github.com/ZimmermanGroup/molecularGSM',
                TSMethodsEnum.kinbot: 'https://github.com/zadorlab/KinBot',
                TSMethodsEnum.neb_ase: 'https://wiki.fysik.dtu.dk/ase/install.html',
                TSMethodsEnum.neb_terachem: 'http://www.petachem.com/products.html',
                TSMethodsEnum.qst2: 'https://gaussian.com/',
                }

    return qcel.util.which(software_dict[method],
                 return_bool=True,
                 raise_error=raise_error,
                 raise_msg=f'Please install {software_dict[method]} to use the {method} method, '
                           f'see {url_dict[method]} for more information',
                 )

