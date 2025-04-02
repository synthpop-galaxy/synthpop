"""
Metallicity class for a double Gaussian distribution, where the two
distributions have a relative weight, upper and lower limits,
and each Gaussian is given a mean and standard deviation.
"""

__all__ = ['DoubleGaussian']
__author__ = "S. Johnson"
__date__ = "2022-07-06"

import math
import numpy as np
from ._metallicity import Metallicity

class DoubleGaussian(Metallicity):
    """
    Double Gaussian metallicity class

    Attributes
    ----------
    metallicity_func_name : string
        A class attribute for the name of the _MetallicityBase subclass that this is.
    weight  : float
        percentage of stars belonging to the first Gaussian distribution ( A1/(A1+A2))
    mean1 : float [[Fe/H]]
        the mean value of the first Gaussian distribution
    std1 : float [[Fe/H]]
        the standard deviation of the first Gaussian distribution
    mean2 : float [[Fe/H]]
        the mean value of the second Gaussian distribution
    std2 : float [[Fe/H]]
        the standard deviation of the second Gaussian distribution
    lower_bound : float [[FE/H]]
        lower limit for truncation of the distribution
    upper_bound : float [[FE/H]] 
        upper limit for truncation of the distribution

    Methods
    -------
    __init__(self,Population) : None
        check that required attributes of Population are present:
            -Population.metallicity_func_kwargs : float
    random_metallicity(self) : float
        return a random metallicity drawn from the weighted double Gaussian
    average_metallicity(self) : float
        return the average metallicity
    """

    def __init__(
            self, weight: float, mean1: float, std1: float, mean2: float, std2: float,
            low_bound: float = -4, high_bound: float = 0.5, **kwargs
            ):
        super().__init__(**kwargs)
        self.metallicity_func_name = 'double-gaussian'

        # check if weight is in the correct range
        if (weight > 1) or (weight < 0):
            raise ValueError(f'{weight = }" has to be between 0 and 1')

        self.weight = weight
        self.mean1 = mean1
        self.std1 = std1
        self.mean2 = mean2
        self.std2 = std2
        self.lower = low_bound
        self.upper = high_bound

        if self.lower >= self.upper:
            raise ValueError(f"{low_bound = } has to be strictly smaller than {high_bound = }")

    def _gen_met(self, N: int = 1) -> np.ndarray:
        """
        generate random metallicities using a mixture methode
        """
        # randomly assign if a star belongs to the first or second distribution
        first_pop = np.random.random(N) < self.weight
        second_pop = np.logical_not(first_pop)
        met_array = np.ones(N)

        # generate stars which belongs to the first gaussian distribution
        met_array[first_pop] = np.random.normal(self.mean1, self.std1, sum(first_pop))

        # generate stars which belongs to the second gaussian distribution
        met_array[second_pop] = np.random.normal(self.mean2, self.std2, sum(second_pop))

        return met_array

    def draw_random_metallicity(self, N: int or None = None, **kwargs) -> np.ndarray or float:
        """
        Returns one or more metallicities in [Fe/H] from a double Gaussian distribution.

        Parameters
        ----------
        N : int, None, optional
            if N is set to an integer, an array with N random ages is returned

        Returns
        -------
        val : ndarray, float [Gyr]
            single metallicities or ndarray of N metallicities in [Fe/H]
        """

        if N is None:
            # generate a single value

            while True:
                val = self._gen_met(1)
                if self.lower < val < self.upper:
                    return val[0]
        else:
            # generate multiple values
            val = self._gen_met(N)
            while True:
                outside = (self.lower > val) | (val > self.upper)
                if not any(outside):
                    return val
                val[outside] = self._gen_met(sum(outside))

    def average_metallicity(self) -> float:
        """Determine the average metallicity of the population"""

        return self.weight * self.mean1 + (1 - self.weight) * self.mean2


    def likelyhood_distribution(self, met: np.ndarray) -> np.ndarray:
        """ analytic version  of likelyhood_distribution. only used for the validating"""

        properbility = self.weight / np.sqrt(2*np.pi)/self.std1 \
                       * np.exp(-0.5*((met-self.mean1)/self.std1)**2) \
                       + (1 - self.weight) / np.sqrt(2*np.pi)/self.std2 \
                       * np.exp(-0.5 * ((met - self.mean2) / self.std2) ** 2)

        # include upper and lower cut
        properbility[met>self.upper] = 0
        properbility[met<self.lower] = 0
        # renormalize
        norm = self.weight * (math.erf((self.upper-self.mean1)/self.std1/np.sqrt(2))
               + math.erf((self.mean1-self.lower)/self.std1/np.sqrt(2)))/2 \
               + (1-self.weight) * (math.erf((self.upper - self.mean2) / self.std2 / np.sqrt(2))
               + math.erf((self.mean2 - self.lower) / self.std2 / np.sqrt(2))) / 2

        return properbility / norm
