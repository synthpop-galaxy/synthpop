"""
This file contains the base class for the initial-final mass relation
"""

__all__ = ["InitialFinalMassRelation"]
__author__ = "M.J. Huston"
__credits__ = ["M.J. Huston"]
__date__ = "2025-10-15"

from typing import Union, Callable
from types import ModuleType

import numpy as np
from scipy import integrate, interpolate
from abc import ABC, abstractmethod


class InitialFinalMassRelation(ABC):
    """
    The initial-final mass relation (IFMR) class for Population class.
    A keyword name is given upon initialization to select the form
    of the IFMR. 
    
    Methods:
    --------
    process_compact_objects(mass, metallicity, age) - returns final mass + phase (remnant type)
    """

    def __init__(self, logger: ModuleType = None, **kwargs):
        """
        Initialize the IFMR class for a Population class
        """
        self.logger = logger

    # This is only a placeholder. The function should be defined in a subclass
    @abstractmethod
    def process_compact_objects(self, m_init: Union[np.ndarray, float],
                                      feh_init: Union[np.ndarray, float]):
        raise NotImplementedError('No IFMR specified')
