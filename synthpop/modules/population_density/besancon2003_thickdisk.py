"""
Thick disk density from the Besancon model (Robin et al., 2003)
"""

__all__ = ["Besancon2003Thickdisk", ]
__author__ = "M.J. Huston"
__date__ = "2022-07-12"

import numpy as np
from .. import const
from ._population_density import PopulationDensity

class Besancon2003Thickdisk(PopulationDensity):
    """
    Besancon Model (Robin et al., 2003) thick disk
    
    Attributes
    ----------
    rho0 : float [Msun/kpc**3]
        density at the location of the sun
    hr : float [kpc]
        scale height in radii
    hz  : float [kpc]
            scale height in height
    xl : float
        transit point
    flare_flag: boolean
        flag if flare is included or not
    """

    def __init__(self, rho0, hr, hz, xl, flare_flag=False, **kwargs):
        super().__init__(**kwargs)
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
            k_flare = self.get_kappa_flare(r)
            k_flare0 = self.get_kappa_flare(self.sun.r)
        else:
            k_flare = 1
            k_flare0 = 1

        # estimate d0 by passing the coordinates of the sun

        # r_sun == -x_sun
        if abs(self.sun.z) < self.xl:
            # most likely
            d0 = (1 - self.sun.z ** 2 / (self.xl * (2 * self.hz/k_flare0 + self.xl)))
        else:
            d0 = np.exp((self.xl - np.abs(self.sun.z)) / self.hz/k_flare0) \
                    / (1 + self.xl / (2 * self.hz/k_flare0))

        # used so numpy arrays, and floats works the same without checking the type of the input
        abs_z_lt_xl = np.abs(z) <= self.xl  # 1 if True, 0 if False

        rho = self.rho0 / d0 * np.exp(-(r - self.sun.r) / self.hr)
        if self.xl>0:
            rho *= (1 - z ** 2 / (self.xl * (2 * self.hz/k_flare + self.xl))
                               ) ** abs_z_lt_xl  # only if abs(z) <= xl

        rho *= (np.exp((self.xl - np.abs(z)) / self.hz / k_flare)
                / (1 + self.xl / 2 / self.hz/k_flare)
                ) ** (1 - abs_z_lt_xl)  # only if abs(z) > xl
        return rho
