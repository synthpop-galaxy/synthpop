"""
Density Subclass to describe the bulge density
from the Besancon model (Robin et al., 2003)
"""
__all__ = ["Besancon2003Bulge", ]
__author__ = "M.J. Huston"
__date__ = "2022-07-12"

import numpy as np
from .. import const
from ._population_density import PopulationDensity

class Besancon2003Bulge(PopulationDensity):
    """
    Besancon model (Robin et al. 2003) bulge density

    Attributes
    ----------
    x0, y0, z0 : float [kpc]
        semi major axis of the bulge
    Rc  : float [kpc]
        critical Radius
    n0 : float [stars/kpc^3]
        central number density
    bar_angle : float [°]
        angle of the bar,
    """
        
    def __init__(self, x0=1.59, y0=0.424, z0=0.424, Rc=2.54, n0=1.37e10, bar_angle=11.1, **kwargs):
        super().__init__()
        self.population_density_name = "Besancon2003Bulge"
        self.density_unit = 'number'
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.Rc = Rc
        self.n0 = n0
        self.bar_ang = bar_angle * np.pi / 180

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
        # align coordinates with the bar,
        # x is pointing away from us.
        xb = -r * np.cos(phi_rad - self.bar_ang)
        yb = r * np.sin(phi_rad - self.bar_ang)
        zb = z

        # calculations
        rs_squared = np.sqrt(((xb / self.x0) ** 2 + (yb / self.y0) ** 2) ** 2 + (zb / self.z0) ** 4)

        rho = self.n0 * np.exp(-0.5 * rs_squared)

        # 1 if sqrt(xb**2+yb**2) < Rc
        rho *= np.exp(-2 * np.maximum(np.sqrt(xb ** 2 + yb ** 2) - self.Rc, 0) ** 2)

        return rho
