"""
Extinction map from Lallement et al 2019, for dust within 3 kpc based on
Gaia & 2MASS observations.

Extinction is given as 'A0', or total extinction at reference wavelength 5500 Angstroms.

The user may determine whether to integrate over the map individually for each star, or to
integrate for the average sightline and interpolate for distances. The user can also adjust
the dr used for extinction integration.

Publication DOI: 10.1051/0004-6361/201834695

Data file FTP: https://cdsarc.cds.unistra.fr/ftp/J/A+A/625/A135/
"""

__all__ = ["Lallement", ]
__author__ = "M.J. Huston"
__date__ = "2024-11-06"

import gzip
import h5py
import shutil
import numpy as np
try: 
    from ._extinction import ExtinctionMap
    from .. import const
except ImportError:
    from _extinction import ExtinctionMap
    import const
import time
import os
import requests

current_map_name = None
current_map_data = None

class Lallement(ExtinctionMap):
    """
    Extinction map from Lallement et al. 2019

    Attributes
    ----------
    dr=0.001 : float
        step size for extinction integration in kpc
    per_sightline=True : boolean
        calculate extinction integral once per sightline if True (faster)
        calculate extinciton integral individually per star if False (slower)

    Methods
    -------
    lallement_ext_func(l_deg, b_deg, dist):
        get extinction value in map for list of star locations
    extinction_in_map(l_deg, b_deg, dist):
        equivalent to lallement_ext_func
    """

    def __init__(self, dr=0.001, per_sightline=True, **kwargs):
        super().__init__(**kwargs)
        # name of the extinction map used
        self.extinction_map_name = "Lallement"
        # A0 = value at 5500 angstroms
        self.ref_wavelength = 0.55
        self.A_or_E_type = 'A0'
        self.per_sightline = per_sightline # Calculate the extinction integral once per sightline, otherwise, individually for each star
        self.dr = dr #: step size for extinction integration in kpc
        if not os.path.isfile(f'{const.EXTINCTIONS_DIR}/map3D_GAIAdr2_feb2019.h5'):
            print("Missing Lallement extinction table. Download and unpacking may take several minutes but only needs done once.")
            map_url = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/625/A135/map3D_GAIAdr2_feb2019.h5.gz'
            map_filename = f'{const.EXTINCTIONS_DIR}/'+map_url.split("/")[-1]
            with open(map_filename, "wb") as f:
                r = requests.get(map_url)
                f.write(r.content)
                print('Map retrieved.')
            with gzip.open(map_filename,'rb') as f_in, open(map_filename[:-3],'wb') as f_out:
                shutil.copyfileobj(f_in,f_out)
            os.remove(map_filename)
            print('File unzipped; ready to use.')
            
        # Check whether the map is already loaded from the prior population
        global current_map_name, current_map_data
        if (current_map_name is not None) and (current_map_data is not None):
            if current_map_name==self.extinction_map_name:
                self.map_data = current_map_data
            else:
                self.map_data = np.array(h5py.File(f'{const.EXTINCTIONS_DIR}/map3D_GAIAdr2_feb2019.h5')['stilism']['cube_datas'])
                current_map_name = self.extinction_map_name
                current_map_data = self.map_data
        else:
            self.map_data = np.array(h5py.File(f'{const.EXTINCTIONS_DIR}/map3D_GAIAdr2_feb2019.h5')['stilism']['cube_datas'])
            current_map_name = self.extinction_map_name
            current_map_data = self.map_data
        # 5 pc grid spacing: -3 to 3 kpc in x,y; -400 to 400 pc in z
        self.grid_dr = 0.005
        self.grid_x_mid = 600
        self.grid_y_mid = 600
        self.grid_z_mid = 80
        self.x_extent=3.0
        self.y_extent=3.0
        self.z_extent=0.4
        # Data units are dmag/dpc at A0, 5500 angstrom
        
    def lallement_ext_func(self,l_deg,b_deg,dist):
        '''
        Get extinction value from Lallement et al. (2019) map
        for an array of star positions.
        
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
        '''
        dist_max = np.max(dist)
        # Take mean sightline for all stars, or do each star individually
        if self.per_sightline:
            dist_pts = np.arange(0,dist_max,self.dr)
            # Convert to nearest neighbor map array indices
            l_rad, b_rad = np.mean(l_deg)*np.pi/180, np.mean(b_deg)*np.pi/180
        else:
            dist_pts = np.arange(0,dist_max,self.dr)[np.newaxis,:]
            # Convert to nearest neighbor map array indices
            l_rad, b_rad = l_deg[:,np.newaxis]*np.pi/180, b_deg[:,np.newaxis]*np.pi/180
        xm_dists = dist_pts*np.cos(b_rad)*np.cos(l_rad)
        ym_dists = dist_pts*np.cos(b_rad)*np.sin(l_rad)
        zm_dists = dist_pts*np.sin(b_rad)
        xm_pts = np.maximum(np.minimum(np.around(xm_dists/self.grid_dr).astype(int), self.grid_x_mid), -self.grid_x_mid) + self.grid_x_mid
        ym_pts = np.maximum(np.minimum(np.around(ym_dists/self.grid_dr).astype(int), self.grid_y_mid), -self.grid_y_mid) + self.grid_y_mid
        zm_pts = np.maximum(np.minimum(np.around(zm_dists/self.grid_dr).astype(int), self.grid_z_mid), -self.grid_z_mid) + self.grid_z_mid
        # Find where each line exits the grid
        within_grid = (np.abs(xm_dists)<self.x_extent) & \
                        (np.abs(ym_dists)<self.y_extent) & \
                        (np.abs(zm_dists)<self.z_extent)
        # Get differential extinction along sightline
        dm_dr_pts = self.map_data[xm_pts,ym_pts,zm_pts]
        # Sum up dA/dr * dr for each point in foreground
        # to get total extinction for sightline or each star
        if self.per_sightline:
            m_pts = np.cumsum(dm_dr_pts*self.dr*1e03*within_grid)
            result = np.interp(dist, dist_pts, m_pts)
        else:
            result = np.sum(dm_dr_pts*self.dr*1.0e3*(dist_pts<(dist[:,np.newaxis])) * within_grid, axis=1)
        return result

    def extinction_in_map(self,l_deg,b_deg,dist):
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
        return self.lallement_ext_func(l_deg, b_deg, dist)
