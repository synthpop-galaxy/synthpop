"""
Metallicity module for a single value
"""

__all__ = ["SingleValue", ]
__author__ = "S. Johnson"
__date__ = "2022-07-06"

import numpy as np
from ._metallicity import Metallicity

class SingleValue(Metallicity):
    """
    Metallicity distribution for a single value

    Attributes
    ----------
    metallicity_value : float [[Fe/H]]
        [Fe/H] value for all stars
    """

    def __init__(self, met_value: float, **kwargs):
        super().__init__(**kwargs)
        self.metallicity_func_name = 'single_value'
        self.metallicity_value = met_value

    def draw_random_metallicity(
            self,
            N: int or None = None,
            **kwargs
            ) -> np.ndarray or float:
        """
        Generate a "random" metallicity from a single value distribution

        Parameters
        ----------
        N : int, None, optional
            if N is set to an integer, an array with N random metallicities is returned
        Returns
        -------
        met : float, ndarray [[Fe/H]]
            single metallicities or numpy array of N metallicities in [[Fe/H]]
        """
        if N is None:
            return self.metallicity_value
        else:
            return np.ones(N) * self.metallicity_value

    def average_metallicity(self) -> float:
        """Determine the average metallicity of the population"""
        return self.metallicity_value
