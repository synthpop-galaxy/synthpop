"""
Postprocessing module to convert from galactic coordinates to equatorial.
"""

__all__ = ["EquatorialCoordinates", ]
__author__ = "M.J. Huston"
__date__ = "2025-06-30"

import pandas as pd
import numpy as np
from ._post_processing import PostProcessing
from synthpop.synthpop_utils.coordinates_transformation import lb_to_ad, uvw_to_vrmuad

class EquatorialCoordinates(PostProcessing):
    """
    Postprocessing module to convert from galactic to equatorial coordinates
    
    Attributes
    ----------
    keep_galactic : bool
        selection for whether to keep galactic coordinates in data frame or remove
    """

    def __init__(self, model, logger, keep_galactic=False, **kwargs):
        super().__init__(model,logger, **kwargs)
        self.keep_galactic = keep_galactic

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Perform the magnitude conversions and returns the modified DataFrame.
        """

        # Get ra, dec from l, b & add to system_df
        l_arr, b_arr = system_df['l'].to_numpy(), system_df['b'].to_numpy()
        ra_arr, dec_arr = lb_to_ad(l_arr,b_arr)
        l_idx = np.where(system_df.columns=='l')[0][0]
        system_df.insert(l_idx, 'ra', ra_arr)
        system_df.insert(l_idx+1, 'dec', dec_arr)

        # Get mura, mudec from coords and kinemaitcs & add to system_df
        dist_arr = system_df['Dist'].to_numpy()
        u_arr, v_arr, w_arr = system_df['U'].to_numpy(), system_df['V'].to_numpy(), system_df['W'].to_numpy()
        vr_arr, mura_arr, mudec_arr = uvw_to_vrmuad(l_arr,b_arr,dist_arr, u_arr,v_arr,w_arr)
        mul_idx = np.where(system_df.columns=='mul')[0][0]
        system_df.insert(mul_idx, 'mura', mura_arr)
        system_df.insert(mul_idx+1, 'mudec', mudec_arr)

        # Remove galactic coordinates if selected
        if not self.keep_galactic:
            system_df.drop(columns=['l','b'],inplace=True)
            system_df.drop(columns=['mul','mub'],inplace=True)
        
        return system_df, companion_df
