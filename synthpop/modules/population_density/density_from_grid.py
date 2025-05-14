"""
Density function that interpolates over a grid
"""

__all__ = ["density_from_grid"]
__author__ = "M.J. Huston"
__date__ = "2024-04-17"

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from ._population_density import PopulationDensity
from .. import const

class density_from_grid(PopulationDensity):
    """
    Generic PopulationDensity subclass to interpolate over a grid

    Attributes
    ----------
    moment_file : str
        name of the file containing a grid with 3 required columns: r, z, and rho
        units must be kpc, kpc, and Msun/pc^3
        file must be whitespace delimited and have comments marked with '#'
    density_unit : str
        "mass" or "number" to specify units for the provided density
    """
    
    def __init__(
            self, moment_file=None, density_unit='mass',abs_z=True,
            **kwargs
            ):
        super().__init__(**kwargs)
        dat = pd.read_csv(const.MOMENTS_DIR + '/' + moment_file,
            sep='\s+', comment='#')
        self.interpolate_rho = LinearNDInterpolator(list(zip(dat['r'],dat['z'])), 
            dat['rho'], fill_value=0.0, rescale=False)
        self.density_unit = density_unit
        self.abs_z=abs_z

    def density(self, r, theta, z):
        if self.abs_z:
            z = np.abs(z)

        return self.interpolate_rho(list(zip(r,z))) * 1e9
