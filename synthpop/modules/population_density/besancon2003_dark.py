"""
Density subclass for dark matter from Besancon model
(Robin et al., 2003)
"""
__all__ = ["DarkDensity", ]
__date__ = "2021-10-18"

import numpy as np
from ._population_density import PopulationDensity

class Besancon2003Dark(PopulationDensity):
    """
    Dark halo density distributions from the Besancon model (Robin et al. 2003).

    Attributes
    ----------
    e : float
        flattening ratio
    pc : float
        central density
    Rc : float
        cut-off radius
    """

    def __init__(self, e=1, pc=1.079e8, Rc=2.697, **kwargs):
        super().__init__(**kwargs)
        self.density_unit = 'mass'
        self.e = e
        self.pc = pc
        self.Rc = Rc

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

        a = np.sqrt(r ** 2 + (z / self.e) ** 2)
        rho = self.pc / (1. + (a / self.Rc) ** 2)

        return rho
