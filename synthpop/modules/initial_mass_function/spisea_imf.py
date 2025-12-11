"""
Initial mass function for a piecewise power law, e.g.:
    
        For M ≤ m1:    ξ(M) ∝ M^-a0
        
        For m1 < M ≤ m2:    ξ(M) ∝ M^-a1
        
        For m2 < M:           ξ(M) ∝ M^-a2
        
The number of splitpoints is modifiable.
"""

__all__ = ["SpiseaImf", ]
__date__ = "2025-12-10"

import numpy as np

try:
    from ._initial_mass_function import InitialMassFunction
except ImportError:
    from _initial_mass_function import InitialMassFunction

from typing import Callable
from spisea.imf import imf as spisea_imf
from spisea.imf import spisea_multiplicity as spisea_multiplicity

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
                    spisea_imf_kwargs={}):
        super().__init__(min_mass, max_mass)
        self.imf_name = 'SpiseaImf'
        self.min_mass=min_mass
        self.max_mass=max_mass
        self.spisea_imf_kwargs = spisea_imf_kwargs
        self.spisea_imf_name = spisea_imf_name
        self.spisea_multiplicity = None
        if spisea_multiplicity_name is not None:
            self.spisea_multiplicity = getattr(spisea_multiplicity, spisea_multiplicity_name)
        self.spisea_imf = getattr(spisea_imf, spisea_imf_name)(massLimits=np.array([min_mass, max_mass]), 
                            multiplicity=self.spisea_multiplicity, **spisea_imf_kwargs)

    # returns the number of stars of that initial mass
    def imf(self, m_in):
        """
        Initial mass function

        Parameters
        ----------
        m_in: initial mass

        Returns
        -------
        prob: probability at the initial mass

        """
        if not isinstance(m_in, np.ndarray):
            m = np.array([m_in])
        else:
            m = m_in

        prob = np.sum([i(m) for i in self.imf_parts], axis=0)

        if not isinstance(m_in, np.ndarray):
            return prob[0]
        return prob

