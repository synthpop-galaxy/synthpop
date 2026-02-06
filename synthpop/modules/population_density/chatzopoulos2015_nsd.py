"""
NSD Density profile from Chatzopoulos et al. (2015)
"""

__all__ = ["Chatzopoulos2015Nsd", ]
__author__ = "M.J. Huston"
__date__ = "2026-02-06"

import numpy as np
import scipy.special
from ._population_density import PopulationDensity

class Chatzopoulos2015Nsd(PopulationDensity):
    """
    NSD density profile
    
    Attributes
    ----------
    gamma : float
    q : float
    a0 : float [kpc]
    M : float [Msun]
    """
        
    def __init__(
            self, gamma=0.07, q=0.28, a0=0.182, M=6.2e9,
            **kwargs):
        super().__init__(**kwargs)
        self.gamma = gamma 
        self.q = q 
        self.a0 = a0 
        self.M = M
        
    @staticmethod
    def a_func(r,z):
        return np.sqrt(r**2 + z**2/self.q**2)

    def density(self, r, phi_rad, z):
        a = self.a_func(r,z)
        rho = (3-self.gamma)*self.M / (4*np.pi*self.q) \
              * self.a0/(a**self.gamma * (a+self.a0)**(4-self.gamma))
        return rho
