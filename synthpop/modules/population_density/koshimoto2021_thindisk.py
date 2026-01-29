"""
Thin disk density profile from Koshimoto et al (2021)
"""

__all__ = ["Koshimoto2021Thindisk", ]
__author__ = "M.J. Huston"
__date__ = "2022-02-02"

import numpy as np
from .. import const
from ._population_density import PopulationDensity

class Koshimoto2021Thindisk(PopulationDensity):
    """
    Thin disk density profile from Koshimoto et al. (2021)
    
    Attributes
    ----------
    R : float [kpc]
        disk scale length
    z_sun : float [kpc]
        disk scale height at the solar position for linear scale height model 
        and everywhere for the flat scale height model
    z_45 : float [kpc]
        disk scale height at 4.5 kpc for linear scale height model
    rho0 : float [M_sum/kpc^3]
        mass density at the solar position
    R_break : float [kpc]
        distance within which surface density is flat
    linear_z : boolean
        if True, use linear scale height; if False, use flat
    """

    def __init__(self, R, z_sun, z_45, rho_sun, R_break=5.3, linear_z=False, **kwargs):
        super().__init__(**kwargs)
        self.density_unit = 'mass'
        self.R = R
        self.z_sun = z_sun
        self.z_45 = z_45
        self.rho_sun = rho_sun
        self.R_break = R_break
        self.linear_z = linear_z

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
        r_greater_45 = (r > 4.5)
        if self.linear_z:
            zR = r_greater_45 * (self.z_sun - (self.z_sun - self.z_45) * (self.sun.r - r) / (
                    self.sun.r - 4.5)) + (1 - r_greater_45) * self.z_45
        else:
            zR = self.z_sun

        r_greater_break = (r > self.R_break)
        rho = self.rho_sun * self.z_sun / zR * (
                    r_greater_break * np.exp(-(r - self.sun.r) / self.R) * (
                        1 / np.cosh(-abs(z) / zR)) ** 2 + (1 - r_greater_break) * np.exp(
                -(self.R_break - self.sun.r) / self.R) * (1 / np.cosh(-abs(z) / zR)) ** 2)

        return rho
