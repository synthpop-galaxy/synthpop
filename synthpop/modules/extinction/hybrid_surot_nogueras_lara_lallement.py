"""
NEED TO UPDATE HERE

Extinction is provided as total extinction A_Ks at 2.15 microns
"""

__all__ = ["HybridSurotNoguerasLaraLallement", ]
__author__ = "M.J. Huston"
__date__ = "2026-09-21"

import numpy as np
import pandas as pd
from .. import const
from scipy.spatial import KDTree
try:
    from ._extinction import ExtinctionMap
    from .lallement2022 import Lallement2022
    from .hybrid_surot_nogueras_lara import HybridSurotNoguerasLara
    from .. import const
except ImportError:
    from _extinction import ExtinctionMap
    from lallement2022 import Lallement2022
    from hybrid_surot_nogueras_lara import HybridSurotNoguerasLara
    import constants as const
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

class HybridSurotNoguerasLaraLallement(HybridSurotNoguerasLara, ExtinctionMap):
    """
    Hybrid extinction map in 3-d

    Methods
    -------
    extinction_in_map(l_deg, b_deg, dist):
        equivalent to lallement_ext_func
    """

    def __init__(self, project_3d=True, dist_2d=8.15, **kwargs):
        super().__init__(**kwargs)
        # name of the extinction map used
        self.extinction_map_name = "HybridSurotNoguerasLaraLallement"
        # effective wavelength for VISTA K_s bandpass (http://svo2.cab.inta-csic.es/theory/fps/index.php?id=Paranal/VISTA.Ks&&mode=browse&gname=Paranal&gname2=VISTA#filter)
        self.ref_wavelength = 2.152152
        self.A_or_E_type = "A_Ks"
        self.lallement_file = f'{const.EXTINCTIONS_DIR}/lallement2022_at_surot_sightlines.h5'
        if not os.path.isfile(self.lallement_file):
            generate_lallement_surot_map(self.lallement_file)
        with h5py.File(self.lallement_file, "r") as f:
            l_grid = f['l_grid'][:]
            b_grid = f['b_grid'][:]
            r_grid = f['r_grid'][:]
            ext_grid = f['ext_grid'][:]
        self.lallement_dist_lim = 3.0
        self.lallement_interp = RegularGridInterpolator((l_grid,b_grid,r_grid),
            ext_grid, bounds_error=False, fill_value=None, method='linear')

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
        # Get Lallement point (or nearest point if too distant)
        ext_vals = self.lallement_interp(np.transpose([use_l, b_deg, np.minimum(dist, self.lallement_dist_lim)]))
        
        # Get the Surot/NL value for each point
        coords = SkyCoord(l=use_l, b=b_deg, unit='degree', frame='galactic')
        in_gns_map = self.gns_wcs.footprint_contains(coords)
        if np.any(in_gns_map):
            surot_gns_vals = np.full(len(b_deg), np.nan)
            gns_coords = coords[in_gns_map]
            pix_coords = self.gns_wcs.world_to_pixel(coords[in_gns_map])
            x_pix,y_pix = (np.round(pix_coords)).astype(int)
            x_pix[x_pix==self.gns_table.shape[1]] -= 1
            y_pix[y_pix==self.gns_table.shape[0]] -= 1
            surot_gns_vals[in_gns_map] = self.gns_table[y_pix,x_pix]
            idx_invalid = np.isnan(surot_gns_vals)
            if np.any(idx_invalid):
                _, min_dist_arg = self.coord_tree.query(np.transpose([use_l[idx_invalid], b_deg[idx_invalid]]))
                surot_gns_vals[idx_invalid] = self.A_Ks_list[min_dist_arg]
        else:
            _, min_dist_arg = self.coord_tree.query(np.transpose([use_l,b_deg]))
            surot_gns_vals = self.A_Ks_list[min_dist_arg]
        
        # Distance scaling
        scale_enddist = self.grid_interpolator_3d(np.transpose([use_l,b_deg, np.full(len(use_l), self.lallement_dist_lim)]))
        scale_value = self.grid_interpolator_3d(np.transpose([use_l,b_deg, np.minimum(dist, self.r_grid[-1])]))
        scale_norm = self.grid_interpolator_3d(np.transpose([use_l,b_deg, self.dist_2d*np.ones(len(use_l))]))
        
        which = dist>self.lallement_dist_lim
        true_scale_to_cen = surot_gns_vals[which] / ext_vals[which]
        disk_scale_to_cen = scale_norm[which] / scale_enddist[which]
        disk_scale_to_star = scale_value[which] / scale_enddist[which]
        assert(np.all(true_scale_to_cen)<=1.0)
        ext_vals[which] *= true_scale_to_cen / disk_scale_to_cen * disk_scale_to_star

        return ext_vals



def generate_lallement_surot_map(filename):
    try:
        from .SODC import SODC
    except ImportError:
        from SODC import SODC
    l_grid = np.linspace(-11,11, 221)
    b_grid = np.linspace(-11,6, 171)
    r_grid = np.linspace(0,3.0, 301)
    
    full_map = Lallement2022()
    law = SODC(R_V=2.5)
    
    vals = []
    
    for l in l_grid:
        vals.append([])
        for b in b_grid:
            val = full_map.extinction_in_map(np.full(len(r_grid), l), np.full(len(r_grid), b), r_grid) * law.Alambda_Aref(2.152152) / law.Alambda_Aref(0.55)
            vals[-1].append(val)

    with h5py.File(filename, "a") as f:
        f.create_dataset("l_grid", data=l_grid)
        f.create_dataset("b_grid", data=b_grid)
        f.create_dataset("r_grid", data=r_grid)
        f.create_dataset("ext_grid", data=np.array(vals))
