""" This file includes a SingleValue Metallicity Distribution """

__all__ = ["SingleValue", ]
__author__ = "S. Johnson"
__date__ = "2022-07-06"
__license__ = "GPLv3"
__version__ = "1.0.0"

import numpy as np
from ._metallicity import Metallicity


class SingleValue(Metallicity):
    """
    Single value subclass of _MetallicityBase class. This subclass is for Populations that
    have metallicity characterized by a single metallicity value.

    Attributes
    ----------
    metallicity_func_name : str
        name of the metallicity subclass set to single_value.
    metallicity_value : float [[Fe/H]]
        The value of the metallicity all stars in this Population are assumed to have
        in [Fe/H].

    Methods
    -------
    __init__(self,**kwargs) : None
        initialize the metallicity class
    draw_random_metallicity(self,  N: int or None = None) :  ndarray or float  [[Fe/H]]
        returns one or more values for the metallicities of a star in [[Fe/H]].
    average_metallicity(self) : float [[Fe/H]]
        returns the average metallicity of the distribution in [[Fe/H]].
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
