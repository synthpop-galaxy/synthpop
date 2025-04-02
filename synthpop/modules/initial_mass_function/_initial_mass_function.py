"""
This file contains the base class for the initial mass function.
"""

__all__ = ["InitialMassFunction"]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-06-29"

from typing import Union, Callable
from types import ModuleType

import numpy as np
from scipy import integrate, interpolate
from abc import ABC, abstractmethod


class InitialMassFunction(ABC):
    """
    The initial mass function (IMF) class for Population class.
    A keyword name is given upon initialization to select the form
    of the initial mass function. Also, can be initialized with
    preselected minimum and maximum masses
    
    Methods:
    --------
    imf(mass) - returns fraction of stars at the initial mass?
    random_mass() - returns a random mass drawn from the selected imf
    average_mass(min_mass, max_mass) - calculate the average mass from
            the selected imf with options for restricting range
    positional_kwargs
    """
    # place to save grid for multiple usage
    F_imf: Callable = None
    F_imf_inverse: Callable = None
    grid_min: float = np.inf
    grid_max: float = -np.inf

    mass_grid: np.ndarray = None
    prob_dens: np.ndarray = None

    def __init__(self, min_mass=None, max_mass=None, logger: ModuleType = None):
        """
        Initialize the IMF class for a Population class
        """
        self.logger = logger
        # default mass limits
        if min_mass is  None: min_mass = 0.01
        if max_mass is  None: max_mass = 100
        self.min_mass = min_mass
        self.max_mass = max_mass

    # returns the number of stars of that initial mass
    # This is only a placeholder. The function should be defined in a subclass

    @abstractmethod
    def imf(self, m: Union[np.ndarray, float]) -> Union[np.ndarray, float]:
        raise NotImplementedError('No IMF specified')

    def interpolate_F(self, min_mass: float, max_mass: float) -> None:
        """
        integrate and interpolator the imf 
        This allows that any function can be given as IMF
        F_imf and F_inf_inverse can be given by the child class, this can speed up the process
        """
        self.logger.info('Start interpolating the IMF to get F_imf and F_inf_inverse')
        step = 1e-4
        # estimate grid size
        grid_min = self.grid_min if self.grid_min is not None else min_mass
        grid_max = self.grid_max if self.grid_max is not None else max_mass
        self.grid_min = np.floor(min(self.min_mass, min_mass, grid_min) / step) * step
        self.grid_max = np.ceil(max(self.max_mass, max_mass, grid_max) / step + 1) * step
        mass_grid = np.arange(self.grid_min, self.grid_max, step)

        # Evaluate probability density function
        prob_dens = self.imf(mass_grid)
        prob_dens = prob_dens / np.sum(prob_dens)

        # use a cumulative sum as integral (all intervals have the same length)
        F_mass_grid = np.zeros(len(mass_grid))
        F_mass_grid[1:] = np.cumsum((prob_dens[1:] + prob_dens[:-1]) / 2)

        # remove points from grid which have a zero likelihood density function
        mass_grid = mass_grid[prob_dens > 0]
        F_mass_grid = F_mass_grid[prob_dens > 0]
        self.F_imf_inverse = interpolate.interp1d(F_mass_grid, mass_grid, kind='cubic')
        self.F_imf = interpolate.interp1d(mass_grid, F_mass_grid, kind='cubic')
        self.logger.info('Done interpolating the IMF')

    def average_mass(self,
            min_mass: Union[np.ndarray, float, None] = None,
            max_mass: Union[float, None] = None
            ) -> float:
        """
        Returns the average mass for the IMF
        """
        # if no mass range is specified, use mass range from __init__
        # still want to keep ability to custom mass range?
        min_mass = min_mass if min_mass is not None else self.min_mass
        max_mass = max_mass if max_mass is not None else self.max_mass
        # just returning the total mass divided by the total number of stars to get the average
        # total mass is \int_min_mass^max_mass m*imf(m)dm
        # total number is \int_min_mass^max_mass imf(m)dm
        # made numerical integral

        total_mass = integrate.quad(lambda m: m * self.imf(m), min_mass, max_mass)[0]

        if self.F_imf is None:
            total_num = integrate.quad(self.imf, min_mass, max_mass)[0]
        else:
            total_num = self.F_imf(max_mass) - self.F_imf(min_mass)
        if total_mass == 0:
            return 0.
        return total_mass / total_num

    def draw_random_mass(
            self,
            min_mass: Union[np.ndarray, float, None] = None,
            max_mass: Union[float, None] = None,
            N: Union[float, None] = None
            ) -> Union[np.ndarray, float]:
        """
        Draw N random masses between minimum mass min_mass and maximum mass max_mass
        if N is None returns the mass as a float
        Evaluate the density profile along a grid of step 1e-4 
        use the trapezoids along the grid as integral
        use a cubic interpolation to get the inverse function
    
        """
        min_mass = min_mass if min_mass is not None else self.min_mass
        max_mass = max_mass if max_mass is not None else self.max_mass
        if (self.F_imf is None) or (np.min(min_mass) < self.grid_min) or (max_mass > self.grid_max):
            self.interpolate_F(np.min(min_mass), max_mass)

        rand = np.random.uniform(self.F_imf(min_mass), self.F_imf(max_mass), N)
        return self.F_imf_inverse(rand)
