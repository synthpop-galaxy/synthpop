"""
This file contains the base class for the Metallicity distributions.
"""
__all__ = ['Metallicity']
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-06-29"

import numpy as np
from abc import ABC, abstractmethod
from types import ModuleType
from .. import default_sun

class Metallicity(ABC):
    """
    Metallicity class for a Population class. he appropriate subclass is
    assigned based on the metallicity_func_kwargs through the "get_subclass" factory.

    Attributes
    ----------
    metallicity_func_name : str
            name of the metallicity Class
    (more attributes are specified in the subclasses)
        
    Methods
    -------
    __init__(self,**kwargs) : None
        initialize the metallicity class
    draw_random_metallicity(self, N: int or None = None,
            x: np.ndarray or float = None, y: np.ndarray or float = None,
            z: np.ndarray or float = None,
            ) : float, ndarray  [[Fe/H]]
        returns one or more random metallicity values in [Fe/H].
        (specified in the subclasses)
    average_metallicity(self) : float [[Fe/H]]
        returns the average metallicity for the distribution in [Fe/H]
        (specified in the subclasses)
    """
    
    def __init__(self,
            sun: ModuleType = None,
            coord_trans: ModuleType = None,
            logger: ModuleType = None,
            **kwargs):
        """
        Initialize the Metallicity class for a Population class.

        Parameters
        ----------
        sun : SunInfo
            location and velocity of the sun and lsr
            see synthpop_utils/sun_info
        coord_trans: ModuleType
            the coordinate transformation package
            see synthpop_utils/coordinate_transformations
        **kwargs : dict, optional
            control keywords for the metallicity class read from the population.json files
        """
        # sun sun sun, here it comes
        self.logger = logger
        self.sun = sun if sun is not None else default_sun
        self.coord_trans = coord_trans
        self.metallicity_func_name = 'None'

    @abstractmethod
    def draw_random_metallicity(self,
                                N: int or None = None, x: np.ndarray or float = None,
                                y: np.ndarray or float = None, z: np.ndarray or float = None,
                                age: np.ndarray or float = None,
                                ) -> np.ndarray or float:
        """
        Generate a random metallicity from the distribution

        Parameters
        ----------
        N : int, None, optional
            if N is set to an integer, an array with N random metallicities is returned
        x, y, z, : float, ndarray [kpc]
            galactocentric cartesian coordinates
            has length N, or is a float if N == None
        age : float, ndarray [Gyr]
            age of the stars
            has length N, or is a float if N == None

        Returns
        -------
        random_metallicity : float, ndarray [[Fe/H]]
            array of random metallicities
        """

        raise NotImplementedError('No Metallicity subclass set')

    def average_metallicity(self) -> float:
        """
        Returns the average metallicity
        """
        raise NotImplementedError('Average_metallicity not set in  subclass set')

    def likelyhood_distribution(self, met: np.ndarray) -> np.ndarray:
        """ analytic version  of likelyhood_distribution. only used for the validating"""
        return 0
