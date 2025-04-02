"""
Metallicity class for a Gaussian distribution,
given a mean, standard deviation, and
upper and lower limits.
"""

__all__ = ['Gaussian']
__author__ = "S. Johnson"
__date__ = "2022-07-06"

import numpy as np
import math
from .. import const
from ._metallicity import Metallicity

class Gaussian(Metallicity):
    """
    Gaussian metallicity distribution

    Attributes
    ----------
    mean : float [[Fe/H]]
        the mean metallicity in [Fe/H] for the Gaussian distribution
    std : float [[Fe/H]]
        the standard deviation of metallicity in [Fe/H] for the Gaussian distribution
    lower_bound : float [[FE/H]] 
        lower limit for truncation of the distribution
    upper_bound : float [[FE/H]] 
        upper limit for truncation of the distribution

    Methods
    -------
    __init__(self,Population) : None
        initialize the metallicity class, and set the control parameters.
    draw_random_metallicity(self, N: int or None = None) : np.ndarray, float [[Fe/H]]
        return a random metallicity drawn from a Gaussian distribution
    average_metallicity(self) : float [[Fe/H]]
        return the average metallicity
    """

    def __init__(
            self, mean: float, std: float, low_bound: float = -4, high_bound: float = 0.5, gradient=0.0, **kwargs
            ):
        super().__init__(**kwargs)
        self.metallicity_func_name = 'gaussian'
        self.mean = mean
        self.std = std
        self.lower = low_bound
        self.upper = high_bound
        self.gradient = gradient

        if self.lower >= self.upper:
            raise ValueError(f"{low_bound = } has to be strictly smaller than {high_bound = }")

    def draw_random_metallicity(self, N: int or None = None, x=None,y=None,z=None, **kwargs) -> np.ndarray or float:
        """
        Returns one or more metallicities in [Fe/H] from a Gaussian distribution.

        Parameters
        ----------
        N : int, None, optional
            if N is set to an integer, an array with N random ages is returned
            
        Returns
        -------
        val : ndarray, float [Gyr]
            single metallicities or ndarray of N metallicities in [Fe/H]
        """
        if not ((x is None) or (y is None) or (z is None)):
            r_kpc = np.sqrt(x**2 + y**2)
        else:
            r_kpc = self.sun.r
        mean = self.mean + self.gradient * (r_kpc - self.sun.r)
        
        if N is None:
            # generate a single value
            while True:
                val = np.random.normal(mean, self.std)
                if self.lower < val < self.upper:
                    return val

        else:
            # generate multiple values
            val = np.random.normal(mean, self.std, N)
            while True:
                outside = (self.lower > val) | (val > self.upper)
                if not any(outside):
                    return val
                val[outside] = np.random.normal(self.mean, self.std, sum(outside))

    def average_metallicity(self) -> float:
        """Determine the average metallicity of the population"""
        return self.mean

    def likelyhood_distribution(self, met: np.ndarray) -> np.ndarray:
        """ analytic version  of likelyhood_distribution. only used for the validating"""

        properbility = np.exp(-0.5*((met-self.mean)/self.std)**2) / np.sqrt(2*np.pi) / self.std

        properbility[met>self.upper] = 0
        properbility[met<self.lower] = 0

        norm = (math.erf((self.upper-self.mean) / self.std / np.sqrt(2))
                + math.erf((self.mean-self.lower) / self.std / np.sqrt(2)) )/2
        return properbility / norm
