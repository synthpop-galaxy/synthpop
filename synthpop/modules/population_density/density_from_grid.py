"""
Density function that interpolates over a grid. Created for NSD and NSC models
tabulated from AGAMA, but can use other files with the same format.
"""

__all__ = ["density_from_grid"]
__author__ = "M.J. Huston"
__date__ = "2024-04-17"

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator
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
    abs_z : str
        if True, take the absolute value of the z coordinate before evaluating
    """
    
    def __init__(
            self, moment_file=None, density_unit='mass', abs_z=True,
            population_density_name='DensityFromGrid',
            **kwargs
            ):
        dat = pd.read_csv(const.MOMENTS_DIR + '/' + moment_file,
            sep='\s+', comment='#')
        super().__init__(max_gc_dist=np.max(dat['r']), **kwargs)
        rho = dat.pivot(index='r', columns='z', values='rho')
        self.piv = rho
        self.interpolate_rho = RegularGridInterpolator((rho.index.to_numpy(),
            rho.columns.to_numpy()), rho.to_numpy(), 
            bounds_error=False, fill_value=0.0)
        self.density_unit = density_unit
        self.abs_z=abs_z
        self.population_density_name = population_density_name

    def density(self, r, theta, z):
        if self.abs_z:
            z = np.abs(z)

        return self.interpolate_rho(list(zip(r,z))) * 1e9
