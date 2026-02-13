"""
NSD Density profile from Launhardt et al. (2002)
"""

__all__ = ["Launhardt2002Nsd", ]
__author__ = "M.J. Huston"
__date__ = "2026-02-06"

import numpy as np
import scipy.special
from ._population_density import PopulationDensity

class Launhardt2002Nsd(PopulationDensity):
    """
    NSD density profile
    
    Attributes
    ----------
    rho1 : float [Msun per kpc^3]
        central density  of the first component
    rho1_rho2 : float 
        ratio of rho1 to rho2, the central density of the second component
    R1, R2 : float [kpc]
        radial scale length of each component
    z0 : float [kpc]
        scale height
    n_R : float
        exponent for radial scaling
    n_z : float 
        exponent for height scaling
    """
        
    def __init__(
            self, rho1=15.2e10, rho1_rho2=3.9, 
            R1=0.120, R2=0.220, z0=0.45,
            n_R=5, n_z=1.4,
            **kwargs):
        super().__init__(**kwargs)
        self.rho1 = rho1
        self.rho2 = rho1/rho1_rho2
        self.R1 = R1 
        self.R2 = R2 
        self.z0 = z0
        self.n_R = n_R 
        self.n_z = n_z

    def density(self, r, phi_rad, z):
        rho = self.rho1 * np.exp(-0.693* ((r/self.R1)**self.n_R + (np.abs(z)/self.z0)**self.n_z)) \
             +self.rho2 * np.exp(-0.693* ((r/self.R2)**self.n_R + (np.abs(z)/self.z0)**self.n_z)) 

        return rho
