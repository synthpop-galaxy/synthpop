"""
This file contains the base class for the population density distributions.
"""
__all__ = ["PopulationDensity", ]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-07-12"

from abc import ABC, abstractmethod
from types import ModuleType
from typing import Tuple
import numpy as np

from .. import const, default_sun

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
    # placeholder for average mass, and emass/imass correction
    # are set when generating a field
    average_mass  = None
    average_mass_coor = None

    def __init__(self,
            sun: ModuleType = None,
            coord_trans: ModuleType = None,
            gamma_flare: float = None,
            h_flare: float = None,
            radius_flare: float = 0,
            logger: ModuleType = None,
            **kwargs):
        """
        Initialize the Population Density class
        SubClasses MUST define a density_unit!

        Parameters
        ----------
        sun : SunInfo
            location and velocity of the sun and lsr
            see synthpop_utils/sun_info
        coord_trans: ModuleType
            the coordinate transformation package
            see synthpop_utils/coordinate_transformations
        gamma_flare, radius_flare: float
            parameters to implement the flare of the milky way

        """
        # sun sun sun, here it comes
        self.logger = logger
        self.sun = sun if sun is not None else default_sun

        self.population_density_name = 'None'
        self.density_unit = "one of 'mass', 'init_mass', or 'number'"
        self.coord_trans = coord_trans

        self.gamma_flare = 0 if gamma_flare is None and h_flare is None else gamma_flare
        self.h_flare = h_flare
        self.radius_flare = radius_flare
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

    def get_kappa_flare(self,
            r_kpc: np.ndarray or float,
            gamma_flare: float = None,
            h_flare: float = None, radius_flare: float = None
            ) -> np.ndarray or float:
        """
        Estimates the correction factor for the Warp.
        The scale height should then be multiplied by kappa_flare

        Parameters
        ----------
        r_kpc : float or ndarray
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

        if gamma_flare is None and h_flare is None:
            gamma_flare = self.gamma_flare
            h_flare = self.h_flare

        if radius_flare is None:
            radius_flare = self.radius_flare

        if gamma_flare is not None:
            return 1 + gamma_flare * np.maximum(r_kpc - radius_flare, 0)

        return np.exp(np.maximum(r_kpc - radius_flare, 0)/h_flare)

    def gradient(self, r_kpc: np.ndarray or float,
                 phi_rad: np.ndarray or float,
                 z_kpc: np.ndarray or float,
                 eps: Tuple[float] = (1e-3, 1e-4, 1e-3)
                 ) -> Tuple[np.ndarray or float, np.ndarray or float, np.ndarray or float]:
        """
        return the gradient at the given location

        Parameters
        ----------
        r_kpc :  float, ndarray [kpc]
            Radius in kpc
        phi_rad : float, ndarray [rad]
            polar angle follows the rotation of the galaxy
            zero point is at the sun.
        z_kpc : float, ndarray [kpc]
            height above/below the galactic plane
        eps: Tuple(float,float,float)
            difference used to estimate the gradient

        Returns
        -------
        dRho_dR :  float, ndarray
            Gradient in Radius
        dRho_dPhi : float, ndarray
            Gradient in polar angle direction
        dRho_dz : float, ndarray
            Gradient in z direction
        """

        dRho_dR = (self.density(r_kpc + eps[0], phi_rad, z_kpc)
                   - self.density(r_kpc - eps[0], phi_rad, z_kpc)) / (2 * eps[0])
        dRho_dPhi = (self.density(r_kpc, phi_rad + eps[1], z_kpc)
                     - self.density(r_kpc, phi_rad - eps[1], z_kpc)) / (2 * eps[1])
        dRho_dz = (self.density(r_kpc, phi_rad, z_kpc + eps[2])
                   - self.density(r_kpc, phi_rad, z_kpc - eps[2])) / (2 * eps[2])

        return dRho_dR, dRho_dPhi, dRho_dz
