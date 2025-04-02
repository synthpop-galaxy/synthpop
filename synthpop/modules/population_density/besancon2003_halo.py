"""
Subclass for a spheroidal density distribution, used e.g.
for the Halo in the Besancon Model (Robin et al., 2003).
"""
__all__ = ["Besancon2003Halo", ]
__author__ = "J. Klüter"
__date__ = "2022-07-12"

import numpy as np
from .. import const
from ._population_density import PopulationDensity

class Besancon2003Halo(PopulationDensity):
    """
    Spheroidal density profile:
    the density is constant within a radius of ac, and then drops using a power law

    Attributes
    ----------
    p0:
        density at the position of the sun
    e : float
        Flattening ratio between x (or y) and z axis
    ac: float
        Radius until the spheroid has a constant density
    power: float
        Slope of the decay. Should be less than 0
    """

    def __init__(self, p0=9.32e3, e=0.76, ac=0.5, power=-2.44, **kwargs):
        super().__init__(**kwargs)
        self.population_density_name = "Spheroid"
        self.density_unit = 'mass'
        self.p0 = p0
        self.e = e
        self.ac = ac
        self.power = power

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
        a0 = np.sqrt(self.sun.r ** 2 + (self.sun.z / self.e) ** 2)

        return self.p0 * (np.maximum(a, self.ac) / a0) ** self.power
