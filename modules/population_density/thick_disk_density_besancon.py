""" Density Subclass to describe the thick disk density
from the  Besancon model robin et al. (2003)"""
__all__ = ["ThickDiskDensityBesancon", ]
__author__ = "M.J. Houston"
__date__ = "2022-07-12"
__license__ = "GPLv3"
__version__ = "1.0.0"

import numpy as np
from .. import const
from ._population_density import PopulationDensity


class ThickDiskDensityBesancon(PopulationDensity):
    """
    This is the Parent Class for population densities functions

    Attributes
    ----------
    population_density_name : str
        name of the population density
    density_unit : str
        specifier if density profile returns "number"-, "mass"- or "init_mass"-density
    (more attributes are specified in the subclasses)

    Methods
    -------
    __init__() :
        Initialize the PopulationDensity class
    density(r: ndarray, theta, ndarray, z: ndarray) : ndarray
        estimate the density at (r,theta,z)
        (specified in the subclasses)
    """

    def __init__(self, rho0, hr, hz, xl, flare_flag=False, **kwargs):
        """
        initializing

        Parameters
        ----------
        rho0 : float [Msun/kpc**3]
            density at the location of the sun
        hr ; float [kpc]
            scale height in radii
        hz  : float [kpc]
                scale height in height
        xl : float
            transit point
        flare_flag: bool
            flag if flare is included or not
        """
        self.population_density_name = "ThickDisk"
        self.density_unit = 'mass'
        self.rho0 = rho0
        self.xl = xl
        self.hr = hr
        self.hz = hz
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
            k_flare = self.get_flare(r)
            k_flare0 = self.get_flare(const.X_SUN)
        else:
            k_flare = 1
            k_flare0 = 1

        # estimate d0 by passing the coordinates of the sun

        # r_sun == -x_sun
        if abs(const.Z_SUN) < self.xl:
            # most likely
            d0 = 1 / k_flare0 * (1 - const.Z_SUN ** 2 / (self.xl * (2 * self.hz + self.xl)))
        else:
            d0 = np.exp((self.xl - np.abs(const.Z_SUN)) / self.hz) / (1 + self.xl / (2 * self.hz))

        # used so numpy arrays, and floats works the same without checking the type of the input
        abs_z_lt_xl = np.abs(z) <= self.xl  # 1 if True, 0 if False

        rho = self.rho0 / d0 * np.exp(-(r - abs(const.X_SUN)) / self.hr)
        rho *= (1 / k_flare * (1 - z ** 2 / (self.xl * (2 * self.hz + self.xl))
                               )) ** abs_z_lt_xl  # only if abs(z) <= xl

        rho *= (np.exp((self.xl - np.abs(z)) / self.hz) / (1 + self.xl / 2 / self.hz)
                ) ** (1 - abs_z_lt_xl)  # only if abs(z) > xl
        return rho
