"""
This file contains the base class for the population density distributions.
"""
__all__ = ["PopulationDensity", ]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__license__ = "GPLv3"
__date__ = "2022-07-12"
__version__ = '1.0.0'

from abc import ABC, abstractmethod
import numpy as np
from .. import const

class PopulationDensity(ABC):
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

    def __init__(self):
        self.population_density_name = 'None'
        self.density_unit = "'mass', 'init_mass', or 'number'"

    @abstractmethod
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
        raise NotImplementedError('No density profile is implemented!')

    @staticmethod
    def get_flare(
            r: np.ndarray or float, gamma_flare: float = None, radius_flare: float = None
            ) -> np.ndarray or float:
        """
        Estimates the correction factor for the Warp.
        The scale height should then be multiplied by kappa_flare

        Parameters
        ----------
        r : float or ndarray
            galactocentric radii
        gamma_flare: float or None
            slope of the flare
            if None use the default value from const

        radius_flare: float or None
            radius when the flare starts
            if None use the default value from const

        Returns
        -------
        kappa_flare : ndarray:
            correction factor for the scale height
        """

        gamma_flare = const.GAMMA_FLARE if gamma_flare is None else gamma_flare
        radius_flare = const.R_FLARE if radius_flare is None else radius_flare

        return 1 + gamma_flare * np.maximum(r - radius_flare, 0)

    @staticmethod
    def warp_correction(r_kpc, phi_rad):
        """Correction for the WARP following Chen X. et al 2019 """
        # only if r >= rWarp, else hw=0
        hw = const.AMP_WARP * np.maximum(r_kpc - const.R_WARP, 0) ** const.ALPHA_WARP
        return hw * np.sin(phi_rad - const.PHI_WARP)