"""
A TS adapter for user guesses.
"""

import itertools

# import zmats

from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.data.rmg import RMGDatabase
from rmgpy.exceptions import ActionError
from rmgpy.reaction import Reaction

from arc.common import colliding_atoms, get_logger, key_by_val
from arc.species.converter import zmat_from_xyz, zmat_to_xyz
from arc.species.zmat import compare_zmats, get_parameter_from_atom_indices, is_angle_linear, up_param

from .ts_adapter import TSAdapter


class UserAdapter(TSAdapter):

    def __init__(self):



    def __repr__(self) -> str:
        """A short representation of the current UserAdapter.

        Returns
        -------
        str
            The desired representation.
        """
        return f"UserAdapter()"




