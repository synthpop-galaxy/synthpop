"""
Extinction map from Surot et al 2020. This is a 2-d map which may be used
as a screen at some distance, or projected into 3-d following the scheme used
in Galaxia (see galaxia_3d module / Sharma et al. 2011).

Extinction is provided as total extinction A_Ks at 2.15 microns

Publication DOI: 10.1051/0004-6361/202038346

Data file FTP: https://cdsarc.cds.unistra.fr/ftp/J/A+A/644/A140/ejkmap.dat.gz
"""

__all__ = ["HybridSurotNoguerasLara", ]
__author__ = "M.J. Huston"
__date__ = "2026-09-21"

import numpy as np
import pandas as pd
from .. import const
from scipy.spatial import KDTree
from ._extinction import ExtinctionMap
import time
from scipy.interpolate import RegularGridInterpolator
import requests
import os
import tarfile
import h5py
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import pdb

current_map_name = None
current_map_data = None

class HybridSurotNoguerasLara(ExtinctionMap):
    """
    Extinction map from Surot et al. (2020), but replacing w/ Nogueras-Lara et al. (2018)
    for the central region their map covers.
    
    Attributes
    ----------
    project_3d=True : boolean
        whether to project the map into 3-d or keep as 2-d
    dist_2d=8.15 : float
        for 3-d: distance where the map value is the true extinction value;
        for 2-d: distance where extinction is applied as a screen

    Methods
    -------
    extinction_in_map(l_deg, b_deg, dist):
        equivalent to lallement_ext_func
    """

    def __init__(self, project_3d=True, dist_2d=8.15, **kwargs):
        super().__init__(**kwargs)
        # name of the extinction map used
        self.extinction_map_name = "HybridSurotNoguerasLara"
        # effective wavelength for VISTA K_s bandpass (http://svo2.cab.inta-csic.es/theory/fps/index.php?id=Paranal/VISTA.Ks&&mode=browse&gname=Paranal&gname2=VISTA#filter)
        self.ref_wavelength = 2.152152
        self.A_or_E_type = "A_Ks"
        self.project_3d = project_3d
        self.dist_2d = dist_2d
        map_url = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/644/A140/ejkmap.dat.gz'

        # Fetch extinction map data if needed
        surot_map_file = f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table.h5'
        if not os.path.isfile(surot_map_file):
            if not os.path.isdir(f'{const.EXTINCTIONS_DIR}'):
                os.mkdir(f'{const.EXTINCTIONS_DIR}')
            print('Retrieving Surot extinction map file. This may take several minutes the first time.')
            try:
                map_url = 'https://lsu.box.com/shared/static/cwitks7jtzne0w32nhc00yj9pcpcrjoe'
                with open(surot_map_file, "wb") as f:
                    r = requests.get(map_url, stream=True)        
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                    print('Map retrieved.')
            except:
                print(f'There was an error fetching the extinction map.')
                raise
                
        # Fetch extinction map data if needed
        gns_map_file = f'{const.EXTINCTIONS_DIR}/NoguerasLara18_AK.fits'
        if not os.path.isfile(gns_map_file):
            if not os.path.isdir(f'{const.EXTINCTIONS_DIR}'):
                os.mkdir(f'{const.EXTINCTIONS_DIR}')
            print('Retrieving Nogueras-Lara+18 extinction map file. This may take several minutes the first time.')
            try:
                map_url = 'http://cdsarc.u-strasbg.fr/ftp/J/A+A/610/A83/fits/AK_28.fit'
                with open(gns_map_file, "wb") as f:
                    r = requests.get(map_url, stream=True)
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                    print('Map retrieved.')
            except:
                print(f'There was an error fetching the Nogueras-Lara+18 extinction map.')
                raise

        if project_3d:
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

        # Grab saved data from last population, if same map used.
        global current_map_name, current_map_data
        if (current_map_name is not None) and (current_map_data is not None):
            if current_map_name==self.extinction_map_name:
                map_data = current_map_data
                self.coord_tree=map_data[0]
                self.A_Ks_list=map_data[1]
                self.gns_wcs=map_data[2]
                self.gns_table=map_data[3]
            else:
                current_map_name = self.extinction_map_name
                # Data from surot_A_Ks_table1.csv
                tmp = pd.read_hdf(f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table.h5', key='data')
                #pd.read_csv(f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table_2D.csv',
                #    usecols=[0, 1, 2], sep=',', names=['l','b','A_Ks'])
                self.coord_tree = KDTree(np.transpose(np.array([tmp['l'],tmp['b']])))
                self.A_Ks_list = np.array(tmp['A_Ks'])
                tmp = fits.open(f'{const.EXTINCTIONS_DIR}/NoguerasLara18_AK.fits')
                self.gns_wcs = WCS(tmp[0].header)
                self.gns_table = tmp[0].data
                self.gns_table[self.gns_table==0.0] = np.nan
                current_map_data = [self.coord_tree, self.A_Ks_list, self.gns_wcs, self.gns_table]
        else:
            current_map_name = self.extinction_map_name
            # Data from surot_A_Ks_table1.csv
            tmp = pd.read_hdf(surot_map_file, key='data')
            #pd.read_csv(f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table_2D.csv',
            #    usecols=[0, 1, 2], sep=',', names=['l','b','A_Ks'])
            self.coord_tree = KDTree(np.transpose(np.array([tmp['l'],tmp['b']])))
            self.A_Ks_list = np.array(tmp['A_Ks'])
            tmp = fits.open(f'{const.EXTINCTIONS_DIR}/NoguerasLara18_AK.fits')
            self.gns_wcs = WCS(tmp[0].header)
            self.gns_table = tmp[0].data
            self.gns_table[self.gns_table==0.0] = np.nan
            current_map_data = [self.coord_tree, self.A_Ks_list, self.gns_wcs, self.gns_table]
                
        if self.project_3d:
            # Set up 3D grid interpolation for distance scaling
            mapfile_3d = h5py.File(f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.h5', 'r')
            map_grid_3d = np.array(mapfile_3d['xmms'])
            map_data_3d = np.array(mapfile_3d['data'])
            mapfile_3d.close()
            l_grid_3d = np.append(np.arange(*map_grid_3d[0]),map_grid_3d[0][1])
            b_grid_3d = np.append(np.arange(*map_grid_3d[1]),map_grid_3d[1][1])
            self.r_grid = 10**np.append(np.arange(*map_grid_3d[2]),map_grid_3d[2][1])
            self.grid_interpolator_3d = RegularGridInterpolator((l_grid_3d,b_grid_3d,self.r_grid), 
                map_data_3d, bounds_error=False, fill_value=None, method='linear')

    def surot_nl_ext(self, l_deg, b_deg, dist):
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
        use_l = l_deg - (l_deg>180)*360
        # Get GNS map values
        coords = SkyCoord(l=use_l, b=b_deg, unit='degree', frame='galactic')
        in_gns_map = self.gns_wcs.footprint_contains(coords)
        if np.any(in_gns_map):
            ext_value = np.full(len(b_deg), np.nan)
            gns_coords = coords[in_gns_map]
            pix_coords = self.gns_wcs.world_to_pixel(coords[in_gns_map])
            x_pix,y_pix = (np.round(pix_coords)).astype(int)
            x_pix[x_pix==self.gns_table.shape[1]] -= 1
            y_pix[y_pix==self.gns_table.shape[0]] -= 1
            ext_value[in_gns_map] = self.gns_table[y_pix,x_pix]
            idx_invalid = np.isnan(ext_value)
            if np.any(idx_invalid):
                _, min_dist_arg = self.coord_tree.query(np.transpose([use_l[idx_invalid], b_deg[idx_invalid]]))
                ext_value[idx_invalid] = self.A_Ks_list[min_dist_arg]
        else:
            _, min_dist_arg = self.coord_tree.query(np.transpose([use_l,b_deg]))
            ext_value = self.A_Ks_list[min_dist_arg]

        # Calculate the scaling from the 3-d interpolator
        if self.project_3d:
            use_l = l_deg + (l_deg<0)*360
            scale_value = self.grid_interpolator_3d(np.transpose([use_l,b_deg, np.minimum(dist, self.r_grid[-1])]))
            scale_norm = self.grid_interpolator_3d(np.transpose([use_l,b_deg, self.dist_2d*np.ones(len(use_l))]))
            scale_factor = scale_value/scale_norm
        else:
            scale_factor = (dist>self.dist_2d)

        return ext_value * scale_factor

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
        return self.surot_nl_ext(l_deg, b_deg, dist)
