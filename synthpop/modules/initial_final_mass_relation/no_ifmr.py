"""
Filler module for no IFMR handling.
"""

__all__ = ["NoIfmr", ]
__author__ = "M.J. Huston"
__date__ = "2025-10-24"

import pandas
import numpy as np
from ._initial_final_mass_relation import InitialFinalMassRelation
from typing import Set, Tuple, Dict, Union

class NoIfmr(InitialFinalMassRelation):
    """
    Post-processing to nan out dim compact objects

    Attributes
    ----------
    """

    def __init__(self, logger, remove=True, **kwargs):
        super().__init__(logger, **kwargs)
        self.name='NoIfmr'

    def process_compact_objects(self, m_init: Union[np.ndarray, float],
                                      feh_init: Union[np.ndarray, float]):
        """
        Get the final masses and types for compact objects
        """
        n_objs = len(m_init)
        return np.repeat(0.0, n_objs), np.repeat(np.nan, n_objs)
