"""
Assign final compact object types and masses based on SPISEA (Hosek et al 2020; Rose et al 2022).
"""

__all__ = ["SpiseaIfmr", ]
__author__ = "M.J. Huston"
__date__ = "2025-12-10"

import pandas
import numpy as np
from ._initial_final_mass_relation import InitialFinalMassRelation
import pdb
from typing import Set, Tuple, Dict, Union
from spisea import ifmr as spisea_ifmr

class SpiseaIfmr(InitialFinalMassRelation):
    """
    Post-processing to account for dim compact objects, based on PopSyCLE (Rose et al 2022).

    Attributes
    ----------
    spisea_ifmr_name='IFMR_N20_Sukhbold' : string
        selected initial-final mass relation, must exist in SPISEA
    """

    def __init__(self, logger, spisea_ifmr_name='IFMR_N20_Sukhbold', **kwargs):
        super().__init__(logger, **kwargs)
        #: initial-final mass relation name to determine compact object masses within SPISEA.
        self.name='SpiseaIfmr'
        self.spisea_ifmr = getattr(spisea_ifmr, spisea_ifmr_name)()

    def process_compact_objects(self, m_init: Union[np.ndarray, float],
                                      feh_init: Union[np.ndarray, float]):
        raise ValueError("SpiseaIfmr is only compatible with SpiseaGenerator, not StarGenerator")
        return
