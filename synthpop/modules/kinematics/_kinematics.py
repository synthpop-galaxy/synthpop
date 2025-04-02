"""
This file contains the base class for the velocity distributions.
"""

__all__ = ["Kinematics", ]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-06-29"

from abc import ABC, abstractmethod
from typing import Tuple, Type
from types import ModuleType
import numpy as np
from .. import default_sun

class Kinematics(ABC):
    """
    The Kinematics base class for a Population class. The appropriate subclass is
    assigned based on the kinematics_func_kwargs through the "get_subclass" factory.

    Attributes
    ----------
    kinematics_func_name : str
        name of the Kinematics Class
    (more attributes are specified in the subclasses)

    Methods
    -------
    __init__(self, Population) : None
        initialize the Kinematics class
    draw_random_velocity(self, x: ndarray, y: ndarray, z: ndarray, mass: ndarray = None,
                all_x: ndarray = None, all_y: ndarray = None, all_z: ndarray = None,
                all_mass: ndarray = None
                ) : ndarray [km/s]
        returns a random velocity of a star in km/s.


    """

    def __init__(self,
            sun: ModuleType = None,
            coord_trans: ModuleType = None,
            density_class: ModuleType = None,
            logger: ModuleType = None,

                 **kwargs):
        """
        Initialize the Kinematics class
        
        Parameters
        ----------
        sun : SunInfo
            location and velocity of the sun and lsr
            see synthpop_utils/sun_info
        coord_trans: ModuleType
            the coordinate transformation package
            see synthpop_utils/coordinate_transformations
        density_class: ModuleType
            the density profile for the current population

        additional Parameters can be set by each subclass.
        """

        self.logger = logger
        self.kinematics_func_name = 'None'

        # sun sun sun, here it comes
        self.sun = sun if sun is not None else default_sun

        self.coord_trans = coord_trans
        self.density_class = density_class

    @abstractmethod
    def draw_random_velocity(
            self, x: np.ndarray or float, y: np.ndarray or float,
            z: np.ndarray or float, mass: np.ndarray or float = None,
            all_x: np.ndarray or float = None, all_y: np.ndarray or float = None,
            all_z: np.ndarray or float = None, all_mass: np.ndarray or float = None,
            all_density_classes: Tuple[Type] = None,
            **kwargs
            ) -> Tuple[np.ndarray, np.ndarray, np.ndarray] or Tuple[float, float, float]:

        """
        Generate a random velocities u,v,w  by using a velocity dispersion
        arguments starting with 'all_' are only provided if the "kinematics_at_the_end"
        keyword in the configuration is set to True

        Parameters
        ----------
        x, y, z : nparray, float [kpc]
            galactocentric coordinates
        mass : ndarray, float , optional [Msun]
            mass of the stars. Has to be the same format and shape as x,y & z
        all_x, all_y, all_z : nparray, optional  [kpc]
            galactocentric coordinates of all stars generated in all populations
        all_mass : nparray, optional  [kpc]
            mass of all the generated stars. Has to have the same length as all_x
        all_density_classes : Tuple[Type]
            the population density class of all population
            see _population_density/_population_density for available functions.

        Returns
        -------
        u, v, w : ndarray [km/s]
            velocity in galactocentric x,y,and z direction.
        """

        raise NotImplementedError('draw_random_velocity is not specified')

    # @abstractmethod
    # def average_velocity(self):
    #     """ Returns the average age of the population """
    #     raise NotImplementedError('draw_random_velocity is not specified')
