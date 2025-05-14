"""
Thick Disk density model from the Gaia Unvierse Model Snapshot
"""

__all__ = ["GumsThickdisk", ]
__date__ = "2023-04-03"
__author__ = "J. Klüter"

import numpy as np
import scipy.special
from ._population_density import PopulationDensity

class GumsThickdisk(PopulationDensity):
    """
    Thick Disk density model from the Gaia Unvierse Model Snapshot

    Attributes
    ----------
    rho0 : float
    hz : float
    hr : float
    flare_flag : boolean
    """

    def __init__(self, rho0, hr, hz, flare_flag=True, **kwargs):
        super().__init__(**kwargs)
        self.name = "GumsThickdisk"
        self.density_unit = 'mass'
        self.rho0 = rho0
        self.hz = hz
        self.hr = hr
        self.flare_flag = flare_flag

    def density(self, r: np.ndarray, phi_rad: np.ndarray, z: np.ndarray) -> np.ndarray:
        """

        Estimates the density at the given position

        Parameters
        ----------
        r : ndarray ['kpc']
            Distance to z axis
        phi_rad : ndarray ['rad']
            azimuth angle of the stars. phi_rad = 0 is pointing towards sun.
        z : height above the galactic plane (corrected for warp of the galaxy)

        Returns
        -------
        rho : ndarray [M_sun/kpc^3 or #/kpc^3]
            density at the given location, either in number density evolved
            mass density or initial mass density should be specified in density_unit.

        """

        if self.flare_flag:
            k_flare = self.get_kappa_flare(r)
            k_flare0 = self.get_kappa_flare(self.sun.r)
        else:
            k_flare = 1
            k_flare0 = 1
        # density normalization at the Sun [solar mass per kpc^3]
        d0 = np.cosh(self.sun.z/(2*self.hz*k_flare0))**(-2)
        rho = self.rho0/d0 * np.exp(-(r - self.sun.r) / self.hr) \
              * np.cosh(z / (2*self.hz*k_flare)) ** (-2)

        return rho
