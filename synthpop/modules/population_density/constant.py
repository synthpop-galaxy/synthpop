"""
Constant density model, for testing purposes only.
"""

__all__ = ["Constant", ]
__author__ = "M.J. Huston"
__date__ = "2026-01-29"

import numpy as np
import scipy.special
from .. import const
from ._population_density import PopulationDensity

class Constant(PopulationDensity):
    """
    Constant density
    
    Attributes
    ----------
    rho : float
        density at all positions in [density_unit]/kpc^2
    """

    def __init__(self, density_unit: str, rho = 1e8, **kwargs):
        super().__init__(**kwargs)
        self.population_density_name = "Constant"
        self.density_unit = density_unit
        self.rho=rho

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
        return np.ones(r.shape)*self.rho
