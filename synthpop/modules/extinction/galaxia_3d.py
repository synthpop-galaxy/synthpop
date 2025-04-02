""" 
Extinction from Galaxia, based on Schlegel et al 1998 2-D map, with 3-D dust disk model

Extinction is given as E(B-V).

Publication DOI: 10.1088/0004-637X/730/1/3 (galaxia), 10.1086/305772 (2-d map)

Data available at: http://bhs.astro.berkeley.edu/GalaxiaData.tar.gz
"""

__all__ = ["Galaxia3D", ]
__author__ = "M.J. Huston"
__date__ = "2024-04-18"

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from .. import const
try:
    from ._extinction import ExtinctionMap
except ImportError:
    from _extinction import ExtinctionMap
import time
import ebf

class Galaxia3D(ExtinctionMap):
    """
    Extinction map used in Galaxia

    Methods
    -------
    extinction_in_map():
        gets extinction value from map for star positions
    """

    def __init__(self, **kwargs):
        # name of the extinction map used
        self.extinction_map_name = "Galaxia3D"
        self.ref_wavelength = 0.4361
        self.ref_wavelength2 = 0.5448
        self.A_or_E_type = 'E(B-V)' 

        # Set up 3D grid
        mapfile_3d = ebf.read(f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.ebf')
        map_grid_3d = mapfile_3d['exmap3d.xmms']
        map_data_3d = mapfile_3d['exmap3d.data']
        # Set up 3D interpolation
        l_grid_3d = np.append(np.arange(*map_grid_3d[0]),map_grid_3d[0][1])
        b_grid_3d = np.append(np.arange(*map_grid_3d[1]),map_grid_3d[1][1])
        self.r_grid = 10**np.append(np.arange(*map_grid_3d[2]),map_grid_3d[2][1])
        self.grid_interpolator_3d = RegularGridInterpolator((l_grid_3d,b_grid_3d,self.r_grid), 
            map_data_3d, bounds_error=False, fill_value=None, method='nearest')

        # Set up 3d grid
        mapfile_2d = ebf.read(f'{const.EXTINCTIONS_DIR}/Galaxia_Schlegel_4096.ebf')
        map_grid_2d = mapfile_2d['exmap2d.xmms']
        map_data_2d = mapfile_2d['exmap2d.data']
        # set up interpolator
        l_grid_2d = np.append(np.arange(*map_grid_2d[0]),map_grid_2d[0][1])
        b_grid_2d = np.append(np.arange(*map_grid_2d[1]),map_grid_2d[1][1])
        self.grid_interpolator_2d = RegularGridInterpolator((l_grid_2d,b_grid_2d), 
            map_data_2d, bounds_error=False, fill_value=None, method='nearest')

    def extinction_in_map(self, l_deg, b_deg, dist):
        """
        Estimates the extinction for a list of star positions.

        Parameters
        ----------
        l_deg: ndarray [degrees]
            galactic longitude
        b_deg: ndarray [degrees]
            galactic latitude
        dist: ndarray [kpc]
            radial distance from the Sun
        
        Returns
        -------
        extinction_value: ndarray [mag]
            extinction at each star position defined as self.A_or_E_type
        """
        use_l = l_deg + (l_deg<0)*360
        mapval_3d = self.grid_interpolator_3d((use_l, b_deg, dist))
        mapval_2d = self.grid_interpolator_2d((use_l, b_deg))

        # 2D scaling
        return mapval_2d*mapval_3d




