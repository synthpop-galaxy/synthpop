"""
Metallicity class for a Gaussian distribution with a gradient,
given a mean, standard deviation, gradient, and
upper and lower limits.
"""

__all__ = ['GaussianGradient']
__author__ = "J. Klüter"
__date__ = "2022-07-06"

import numpy as np
from ._metallicity import Metallicity
from .. import const

class GaussianGradient(Metallicity):
    """
    Gaussian subclass of Metallicity base class. This subclass is for populations that
    have metallicity characterized by a gaussian metallicity distribution plus a radius gradient.

    Attributes
    ----------
    mean : float [[Fe/H]]
        the mean metallicity in [Fe/H] for the Gaussian distribution at the Solar position
    std : float [[Fe/H]]
        the standard deviation of metallicity in [Fe/H] for the Gaussian distribution
    lower_bound : float [[Fe/H]]
        lower limit for truncation of the distribution
    upper_bound : float [[Fe/H]]
        upper limit for truncation of the distribution
    radial_gradient : float [[Fe/H]/kpc]
        metallicity gradient with Galactocentric radius

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
            self, mean: float, std: float, low_bound: float = -4, high_bound: float = 0.5,
            radial_gradient: float = 0.0,  **kwargs):
        super().__init__(**kwargs)
        self.metallicity_func_name = 'gaussian_gradient'
        self.mean = mean
        self.std = std
        self.lower = low_bound
        self.upper = high_bound
        self.radial_gradient = radial_gradient
        if self.lower >= self.upper:
            raise ValueError(f"{low_bound = } has to be strictly smaller than {high_bound = }")

    def draw_random_metallicity(self, N, x,y,z, **kwargs) -> np.ndarray or float:
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
        radius = np.sqrt(x**2+y**2)

        if N is None:
            # generate a single value
            while True:
                val = np.random.normal(self.mean, self.std) + self.radial_gradient * (radius - self.sun.r)
                if self.lower < val < self.upper:
                    return val
        else:
            # generate multiple values
            val = np.random.normal(self.mean, self.std, N) + self.radial_gradient * (radius - self.sun.r)
            while True:
                outside = (self.lower > val) | (val > self.upper)
                if not any(outside):
                    return val
                val[outside] = np.random.normal(self.mean, self.std, sum(outside)) + self.radial_gradient * (radius[outside] - self.sun.r)
        return val

    def average_metallicity(self) -> float:
        """Determine the average metallicity of the population"""
        return self.mean
