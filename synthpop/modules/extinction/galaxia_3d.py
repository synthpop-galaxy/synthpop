""" 
Extinction from Galaxia, based on Schlegel et al 1998 2-D map, with 3-D dust disk model

Extinction is given as E(B-V).

Publication DOI: 10.1088/0004-637X/730/1/3 (galaxia), 10.1086/305772 (2-d map)

Data available at: http://bhs.astro.berkeley.edu/GalaxiaData.tar.gz
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
import ebf
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

    def __init__(self, use_h5=False, **kwargs):
        # name of the extinction map used
        self.extinction_map_name = "Galaxia3D"
        self.ref_wavelength = 0.4361
        self.ref_wavelength2 = 0.5448
        self.A_or_E_type = 'E(B-V)'
        
        if use_h5:
            if (not os.path.isfile(f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.h5')) or (not os.path.isfile(f'{const.EXTINCTIONS_DIR}/Galaxia_Schlegel_4096.h5')):
                OSError('Missing Galaxia map data in h5 format. Try using ebf or getting the data files.')
        else:
            if (not os.path.isfile(f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.ebf')) or (not os.path.isfile(f'{const.EXTINCTIONS_DIR}/Galaxia_Schlegel_4096.ebf')):
                print("Missing Galaxia map data - download and arrangement will take a few minutes.")
                if not os.path.isdir(f'{const.EXTINCTIONS_DIR}'):
                    os.mkdir(f'{const.EXTINCTIONS_DIR}')
                if not os.path.isfile(f'{const.EXTINCTIONS_DIR}/GalaxiaData.tar.gz'):
                    with open(f'{const.EXTINCTIONS_DIR}/GalaxiaData.tar.gz', "wb") as f:
                        r = requests.get("http://bhs.astro.berkeley.edu/GalaxiaData.tar.gz")
                        f.write(r.content)
                        print('Data retrieved.')
                with tarfile.open(f'{const.EXTINCTIONS_DIR}/GalaxiaData.tar.gz', "r") as f:
                    f.extract('GalaxiaData/Extinction/ExMap3d_1024.ebf', f'{const.EXTINCTIONS_DIR}/')
                    f.extract('GalaxiaData/Extinction/Schlegel_4096.ebf', f'{const.EXTINCTIONS_DIR}/')
                    os.rename(f'{const.EXTINCTIONS_DIR}/GalaxiaData/Extinction/ExMap3d_1024.ebf',
                                f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.ebf')
                    os.rename(f'{const.EXTINCTIONS_DIR}/GalaxiaData/Extinction/Schlegel_4096.ebf',
                                f'{const.EXTINCTIONS_DIR}/Galaxia_Schlegel_4096.ebf')
                os.remove(f'{const.EXTINCTIONS_DIR}/GalaxiaData.tar.gz')
                os.rmdir(f'{const.EXTINCTIONS_DIR}/GalaxiaData/Extinction')
                os.rmdir(f'{const.EXTINCTIONS_DIR}/GalaxiaData')
                print('Extinction file setup complete.')

        # Set up 3D grid interpolation for distance scaling
        if use_h5:
            mapfile_3d = h5py.File(f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.h5', 'r')
            map_grid_3d = np.array(mapfile_3d['xmms'])
            map_data_3d = np.array(mapfile_3d['data'])
            mapfile_3d.close()
        else:
            mapfile_3d = ebf.read(f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.ebf')
            map_grid_3d = mapfile_3d['exmap3d.xmms']
            map_data_3d = mapfile_3d['exmap3d.data']
        l_grid_3d = np.append(np.arange(*map_grid_3d[0]),map_grid_3d[0][1])
        b_grid_3d = np.append(np.arange(*map_grid_3d[1]),map_grid_3d[1][1])
        self.r_grid = 10**np.append(np.arange(*map_grid_3d[2]),map_grid_3d[2][1])
        self.grid_interpolator_3d = RegularGridInterpolator((l_grid_3d,b_grid_3d,self.r_grid), 
            map_data_3d, bounds_error=False, fill_value=None, method='linear')

        # Set up 2D grid interpolation for Schlegel extinction map
        if use_h5:
            mapfile_2d = h5py.File(f'{const.EXTINCTIONS_DIR}/Galaxia_Schlegel_4096.h5', 'r')
            map_grid_2d = np.array(mapfile_2d['xmms'])
            map_data_2d = np.array(mapfile_2d['data'])
            mapfile_2d.close()
        else:
            mapfile_2d = ebf.read(f'{const.EXTINCTIONS_DIR}/Galaxia_Schlegel_4096.ebf')
            map_grid_2d = mapfile_2d['exmap2d.xmms']
            map_data_2d = mapfile_2d['exmap2d.data']
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




