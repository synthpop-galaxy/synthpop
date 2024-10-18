__all__ = ['CustomMetallicity', ]

import numpy as np
from ._metallicity import Metallicity


class CustomMetallicity(Metallicity):
    """
    A template CustomMetallicity subclass that is a copy of the Gaussian subclass.
    Gaussian subclass of _MetallicityBase class. This sublass is for Populations that 
    have metallicity characterized by a single metallicity value. 

    Attributes
    ----------
    metallicity_func_name : string
        A class attribute for the name of the _MetallicityBase subclass that this is.
    mean : float [[Fe/H]]
        the mean metallicity in [Fe/H] for the Gaussian distribution
    std : float [[Fe/H]]
        the standard deviation of metallicity in [Fe/H] for the Gaussian distribution
    
    Methods
    -------
    __init__(self,Population) : None
        check that required attributes of Population are present:
            -Population.metallicity_func_kwargs : float
    random_metallicity(self) : float
        return a random metallicity drawn from a Gaussian distribution
    average_metallicity(self) : float
        return the average metallicity
    """

    def __init__(self, mean, std, **kwargs):
        """Init"""
        super().__init__(**kwargs)
        self.mean = mean
        self.std = std

    def draw_random_metallicity(self, N=1, x=None, y=None, z=None, age=None):
        """Returns a random metallicity for a star in [Fe/H] from a Gaussian distriution"""
        return np.random.normal(self.mean, self.std, N)

    def average_metallicity(self, ):
        """Determine the average metallicity of the population"""
        return self.mean
