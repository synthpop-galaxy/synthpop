"""
NSD Density profile from Gallego-Cano et al. (2020)
"""

__all__ = ["GallegoCano2020Nsd", ]
__author__ = "M.J. Huston"
__date__ = "2026-02-06"

import numpy as np
import scipy.special
from ._population_density import PopulationDensity

class GallegoCano2020Nsd(PopulationDensity):
    """
    NSD density profile
    
    Attributes
    ----------
    rho2 : float [Msun per kpc^3]
        central density  of the second component
    rho1_rho2 : float 
        ratio of rho1 to rho2, the central density of the second component
    R1, R2 : float [kpc]
        radial scale length of each component
    q : float
        uhh
    n1 : float
        exponent for radial scaling
    n2 : float 
        exponent for height scaling
    """
        
    def __init__(
            self, rho2=170e10, rho1_rho2=1.311, 
            R1=0.00506, R2=0.0246, 
            q=0.37, n1=0.72, n2=0.79,
            **kwargs):
        super().__init__(**kwargs)
        self.rho1 = rho1_rho2*rho2
        self.rho2 = rho2
        self.R1 = R1 
        self.R2 = R2 
        self.q = q
        self.n1 = n1 
        self.n2 = n2

        
    @staticmethod
    def a_func(r,z):
        return np.sqrt(r**2 + z**2/self.q**2)

    def density(self, r, phi_rad, z):
        a = a_func(r,z)
        rho = self.rho1*np.exp(-(a/self.R1)**self.n1) + self.rho2*np.exp(-(a/self.R2)**self.n2)

        return rho
