"""
Extinction maps from Surot et al 2020
"""

__all__ = ["Surot", ]
__author__ = "M.J. Huston"
__date__ = "2023-06-26"
__license__ = "GPLv3"
__version__ = "1.0.0"

import numpy as np
import pandas as pd
from .. import const
from scipy.spatial import KDTree
from ._extinction import ExtinctionMap
import time
import ebf
from scipy.interpolate import RegularGridInterpolator
import requests
import os

current_map_name = None
current_map_data = None

class Surot(ExtinctionMap):
    """
    Extinction map from Surot et al. 2020
    On first use, the map is downloaded and converted from E(J-Ks) to A_Ks values
    
    Attributes
    ----------
    extinction_map_name : str
        name of the Extinction Map
    ref_wavelength : float
        reference wavelength for the extinction
    A_or_E_type : str
        Output type from the extinction map.
        If it starts with "A", A_or_E is handled  as a total extinction.
        If it starts with "E": A_or_E is handled as a color excess.
    project_3d : booolean
        True = project the map into 3-D with the Galaxia scheme, with the 
            extinction scaled to the 2-D map value at dist_2d
        False = model extinction as a flat screen at dist_2d
    dist_2d : float
        distance where the 2-D map value is the true value

    Methods
    -------
    extinction_in_map():
        function that returns extinction at input positions
    get_map_properties():
        returns the basic parameters of the extinction map
        used for Communication between ExtinctionLaw and ExtinctionMap
    """

    def __init__(self, project_3d=True, dist_2d=8.15, **kwargs):
        super().__init__(**kwargs)
        # name of the extinction map used
        self.extinction_map_name = "Surot"
        # effective wavelength for VISTA K_s bandpass (http://svo2.cab.inta-csic.es/theory/fps/index.php?id=Paranal/VISTA.Ks&&mode=browse&gname=Paranal&gname2=VISTA#filter)
        self.ref_wavelength = 2.152152
        self.A_or_E_type = "A_Ks"
        self.project_3d = project_3d
        self.dist_2d = dist_2d
        map_url = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/644/A140/ejkmap.dat.gz'

        # Fetch extinction map data if needed
        if not os.path.isfile(f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table.h5'):
            if not os.path.isfile(f'{const.EXTINCTIONS_DIR}/surot_'+map_url.split("/")[-1]):
                print("Missing Surot table. Download and formatting may take several minutes.")
                print('Downloading map file from VizieR...')
                map_filename = f'{const.EXTINCTIONS_DIR}/surot_'+map_url.split("/")[-1]
                with open(map_filename, "wb") as f:
                    r = requests.get(map_url)
                    f.write(r.content)
                    print('Map retrieved.')
            else:
                map_filename = f'{const.EXTINCTIONS_DIR}/surot_ejkmap.dat'
            print('Reading table...')
            E_JKs_map_df = pd.read_fwf(map_filename,compression='gzip', header=None)
            print('Reformatting values...')
            E_JKs_map = E_JKs_map_df.to_numpy()
            A_Ks_vals = 0.422167 * E_JKs_map[:,2] #conversion from Surot2020 paper
            entries = E_JKs_map.shape[0]
            print('Saving hdf5 version')
            map_output = 'surot_A_Ks_table.h5'
            surot_2d = np.zeros((entries, 3))
            surot_2d[:,0] = E_JKs_map[:,0]
            surot_2d[:,1] = E_JKs_map[:,1]
            surot_2d[:,2] = A_Ks_vals
            surot_2d_df = pd.DataFrame(surot_2d, columns=['l','b','A_Ks'])
            surot_2d_df.to_hdf(f'{const.EXTINCTIONS_DIR}/'+map_output, key='data', index=False, mode='w')
            print('File 2D version saved as '+map_output)

        # Grab saved data from last population, if same map used.
        global current_map_name, current_map_data
        if (current_map_name is not None) and (current_map_data is not None):
            if current_map_name==self.extinction_map_name:
                map_data = current_map_data
                self.coord_tree=map_data[0]
                self.A_Ks_list=map_data[1]
            else:
                current_map_name = self.extinction_map_name
                # Data from surot_A_Ks_table1.csv
                tmp = pd.read_hdf(f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table.h5', key='data')
                #pd.read_csv(f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table_2D.csv',
                #    usecols=[0, 1, 2], sep=',', names=['l','b','A_Ks'])
                self.coord_tree = KDTree(np.transpose(np.array([tmp['l'],tmp['b']])))
                self.A_Ks_list = np.array(tmp['A_Ks'])
                current_map_data = [self.coord_tree, self.A_Ks_list]
        else:
            current_map_name = self.extinction_map_name
            # Data from surot_A_Ks_table1.csv
            tmp = pd.read_hdf(f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table.h5', key='data')
            #pd.read_csv(f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table_2D.csv',
            #    usecols=[0, 1, 2], sep=',', names=['l','b','A_Ks'])
            self.coord_tree = KDTree(np.transpose(np.array([tmp['l'],tmp['b']])))
            self.A_Ks_list = np.array(tmp['A_Ks'])
            current_map_data = [self.coord_tree, self.A_Ks_list]
        
        if self.project_3d:
        # Set up 3D grid
            mapfile_3d = ebf.read(f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.ebf')
            map_grid_3d = mapfile_3d['exmap3d.xmms']
            map_data_3d = mapfile_3d['exmap3d.data']
            # Set up 3D interpolation
            l_grid_3d = np.append(np.arange(*map_grid_3d[0]),map_grid_3d[0][1])
            b_grid_3d = np.append(np.arange(*map_grid_3d[1]),map_grid_3d[1][1])
            self.r_grid = 10**np.append(np.arange(*map_grid_3d[2]),map_grid_3d[2][1])
            self.grid_interpolator_3d = RegularGridInterpolator((l_grid_3d,b_grid_3d,self.r_grid),
                map_data_3d)

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
        use_l = l_deg - (l_deg>180)*360
        _, min_dist_arg = self.coord_tree.query(np.transpose([use_l,b_deg]))
        ext_value = self.A_Ks_list[min_dist_arg]

        # Calculate the scaling from the 3-d interpolator
        if self.project_3d:
            use_l = l_deg + (l_deg<0)*360
            scale_value = self.grid_interpolator_3d(np.transpose([use_l,b_deg, dist]))[0]
            scale_norm = self.grid_interpolator_3d(np.transpose([use_l,b_deg, self.dist_2d*np.ones(len(use_l))]))[0]
            scale_factor = scale_value/scale_norm
        else:
            scale_factor = (dist>self.dist_2d)

        return ext_value * scale_factor
