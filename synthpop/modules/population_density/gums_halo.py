"""
Spheroidal density distribution, used for the Halo in the Gaia Universe Model Snapshot
and other models.
"""

__all__ = ["GumsHalo", ]
__author__ = "J. Klüter"
__date__ = "2023-04-03"


import numpy as np
from .. import const
from ._population_density import PopulationDensity


class GumsHalo(PopulationDensity):
    def __init__(self, rho0, e=0.774, rc=2.175827881, ah=1, dh=2.777231, **kwargs):
        """
        Spheroidal density profile
        e.g. To describe the Halo density in the Besancon Gaia Universe Model

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
        super().__init__(**kwargs)
        self.population_density_name = "GumsHalo"
        self.density_unit = 'mass'
        self.rho0 = rho0
        self.e = e
        self.rc = rc
        self.ah = ah
        self.dh = dh

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

        d0 = a0 ** (-self.ah) * (self.rc+a0) ** (-self.dh)
        rho = self.rho0 / d0 * a ** (-self.ah) * (self.rc + a) ** (-self.dh)

        return rho
