"""
This file contains the base class for multiplicity
"""

__all__ = ["Multiplicity"]
__author__ = "M.J. Huston"
__credits__ = ["M.J. Huston"]
__date__ = "2025-10-15"

from typing import Union, Callable
from types import ModuleType

import numpy as np
from scipy import integrate, interpolate
from abc import ABC, abstractmethod


class Multiplicity(ABC):
    """
    Multiplicity base class
    
    Methods:
    --------
    generate_companions
    """

    def __init__(self, logger: ModuleType = None, **kwargs):
        """
        Initialize the IFMR class for a Population class
        """
        self.logger = logger

    # This is only a placeholder. The function should be defined in a subclass
    @abstractmethod
    def generate_companions(self, m_init: Union[np.ndarray, float],
                                      feh_init: Union[np.ndarray, float]):
        raise NotImplementedError('No Multiplicity specified')
