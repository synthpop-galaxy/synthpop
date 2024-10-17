""" Thin disk density profile from 
Koshimoto et al 2021 """
__all__ = ["Koshimoto2021Thindisk", ]
__author__ = "M.J. Huston"
__date__ = "2022-02-02"
__license__ = "GPLv3"
__version__ = "1.0.0"

import numpy as np
from .. import const
from ._population_density import PopulationDensity


class Koshimoto2021Thindisk(PopulationDensity):
    def __init__(self, R0, z0, z45, rho0, Rbreak=5.3, **kwargs):
        """
        r, theta, and z are Galactocentric cylindrical coordinates. [r], [z] = kpc, kpc: [phi_rad] = radians

        Based on Koshimoto et al. 2021

        """
        super().__init__(**kwargs)
        self.density_unit = 'mass'
        self.R0 = R0
        self.z0 = z0
        self.z45 = z45
        self.rho0 = rho0
        self.Rbreak = Rbreak

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
        r_greater_45 = r >= 4.5
        zR = r_greater_45 * (self.z0 - (self.z0 - self.z45) * (self.sun.r - r) / (
                    self.sun.r - 4.5)) + (1 - r_greater_45) * self.z45

        r_greater_break = r >= self.Rbreak 
        rho = self.rho0 * self.z0 / zR * (
                    r_greater_break * np.exp(-(r - self.sun.r) / self.R0) * (
                        1 / np.cosh(-abs(z) / zR)) ** 2 + (1 - r_greater_break) * np.exp(
                -(self.Rbreak - self.sun.r) / self.R0) * (1 / np.cosh(-abs(z) / zR)) ** 2)

        return rho
