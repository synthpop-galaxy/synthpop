""" 
Extinction from Galaxia, based on Schlegel et al 1998 2-D map, with 3-D dust disk model

Extinction is given as E(B-V).

Publication DOI: 10.1088/0004-637X/730/1/3 (galaxia), 10.1086/305772 (2-d map)

Data available in .ebf form at: http://bhs.astro.berkeley.edu/GalaxiaData.tar.gz
"""

__all__ = ["Galaxia_3d", ]
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
import tarfile
import os
import requests
import h5py

class Galaxia_3d(ExtinctionMap):
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

        # Check for files and download if needed
        if not os.path.isdir(f'{const.EXTINCTIONS_DIR}'):
            os.mkdir(f'{const.EXTINCTIONS_DIR}')
        map_filename_3d = f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.h5'
        if not os.path.isfile(map_filename_3d):
            print('Fetching 3-d extinction map file.')
            map_url = 'https://lsu.box.com/shared/static/3hifnqy9u6lko3ebqdqwwt5ockim27t3'
            try:
                with open(map_filename_3d, "wb") as f:
                    r = requests.get(map_url, stream=True)        
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                    print('Map retrieved.')
            except:
                print(f'There was an error fetching the extinction map. Try downloading the file at {map_url} manually and placing it in {const.EXTINCTIONS_DIR}')
                raise
        map_filename_schlegel = f'{const.EXTINCTIONS_DIR}/Galaxia_Schlegel_4096.h5'
        if not os.path.isfile(map_filename_schlegel):
            print('Fetching Schlegel extinction map file.')
            map_url = 'https://lsu.box.com/shared/static/x3he8q1ybjrg4x98le551dmczj3ys1dx'
            try:
                with open(map_filename_schlegel, "wb") as f:
                    r = requests.get(map_url, stream=True)        
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                    print('Map retrieved.')
            except:
                print(f'There was an error fetching the extinction map. Try downloading the file at {map_url} manually and placing it in {const.EXTINCTIONS_DIR}')
                raise

        # Set up 3D grid interpolation for distance scaling
        mapfile_3d = h5py.File(map_filename_3d, 'r')
        map_grid_3d = np.array(mapfile_3d['xmms'])
        map_data_3d = np.array(mapfile_3d['data'])
        mapfile_3d.close()
        l_grid_3d = np.append(np.arange(*map_grid_3d[0]),map_grid_3d[0][1])
        b_grid_3d = np.append(np.arange(*map_grid_3d[1]),map_grid_3d[1][1])
        self.r_grid = 10**np.append(np.arange(*map_grid_3d[2]),map_grid_3d[2][1])
        self.grid_interpolator_3d = RegularGridInterpolator((l_grid_3d,b_grid_3d,self.r_grid), 
            map_data_3d, bounds_error=False, fill_value=None, method='linear')

        # Set up 2D grid interpolation for Schlegel extinction map
        mapfile_2d = h5py.File(map_filename_schlegel, 'r')
        map_grid_2d = np.array(mapfile_2d['xmms'])
        map_data_2d = np.array(mapfile_2d['data'])
        mapfile_2d.close()
        l_grid_2d = np.append(np.arange(*map_grid_2d[0]),map_grid_2d[0][1])
        b_grid_2d = np.append(np.arange(*map_grid_2d[1]),map_grid_2d[1][1])
        self.grid_interpolator_2d = RegularGridInterpolator((l_grid_2d,b_grid_2d), 
            map_data_2d, bounds_error=False, fill_value=None, method='linear')

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