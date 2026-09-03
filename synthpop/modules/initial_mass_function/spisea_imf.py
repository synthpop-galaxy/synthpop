"""
Initial mass function for a piecewise power law, e.g.:
    
        For M ≤ m1:    ξ(M) ∝ M^-a0
        
        For m1 < M ≤ m2:    ξ(M) ∝ M^-a1
        
        For m2 < M:           ξ(M) ∝ M^-a2
        
The number of splitpoints is modifiable.
"""

__all__ = ["SpiseaImf", ]
__date__ = "2025-12-10"
__author__ = "M.J. Huston"

import numpy as np

try:
    from ._initial_mass_function import InitialMassFunction
except ImportError:
    from _initial_mass_function import InitialMassFunction

from typing import Callable
from spisea.imf import imf as spisea_imf
from spisea.imf import multiplicity as spisea_multiplicity

class SpiseaImf(InitialMassFunction):
    """
    Initial mass function generator for a piecewise power law
    
    Attributes:
    -----------
    min_mass : float
        lower initial mass limit
    max_mass : float
        upper initial mass limit
    alphas : ndarray [float]
        power law indices for the mass chunks from lower mass to higher
    splitpoints : ndarray [float]
        mass values where pieces split;
        length must be length alphas minus 1
    """

    def __init__(self, min_mass=None, max_mass=None, spisea_imf_name='Kroupa_2001', spisea_multiplicity_name=None,
                    spisea_imf_kwargs={}, **kwargs):
        super().__init__(min_mass, max_mass)
        self.imf_name = 'SpiseaImf'
        self.spisea_imf_kwargs = spisea_imf_kwargs
        self.spisea_imf_name = spisea_imf_name
        self.spisea_multiplicity = None
        if spisea_multiplicity_name is not None:
            self.spisea_multiplicity = getattr(spisea_multiplicity, spisea_multiplicity_name)
        self.spisea_imf = getattr(spisea_imf, spisea_imf_name)(
                            multiplicity=self.spisea_multiplicity, **spisea_imf_kwargs)
        if self.min_mass is None:
            self.min_mass = self.spisea_imf._mass_limits[0]
        if self.max_mass is None:
            self.max_mass = self.spisea_imf._mass_limits[-1]

    def imf(self, m: Union[np.ndarray, float]) -> Union[np.ndarray, float]:
        raise NotImplementedError("SpiseaImf cannot be used by StarGenerator.")
        #return self.spisea_imf.xi(m)
