""" Thick disk density profile from
Koshimoto et al. 2021 """
__all__ = ["Koshimoto2021Thickdisk", ]
__author__ = "M.J. Huston"
__date__ = "2022-02-02"
__license__ = "GPLv3"
__version__ = "1.0.0"

import numpy as np
from .. import const
from ._population_density import PopulationDensity


# TODO is there a better name?
class Koshimoto2021Thickdisk(PopulationDensity):
    def __init__(
            self, rho0=(1.7e-3 + 4.4e-4) * 10 ** 9, z0=0.903, R0=2.200, Rbreak=5.300, **kwargs
            ):
        """
        set up the thick disc density profile
        default values corresponds to the ThickDisk density profile in Koshimoto et al. 2021
        """
        super().__init__(**kwargs)
        self.density_unit = 'mass'
        self.rho0 = rho0
        self.z0 = z0
        self.R0 = R0
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
        r_past_break = (np.array(r) > self.Rbreak).astype(int)
        rho_past = self.rho0 * np.exp(-(r - self.sun.r) / self.R0) * np.exp(-abs(z) / self.z0)
        rho_pre = self.rho0 * np.exp(-(self.Rbreak - self.sun.r) / self.R0) * np.exp(
            -abs(z) / self.z0)
        return rho_past * r_past_break + rho_pre * (1 - r_past_break)
