"""
Postprocessing module to convert from galactic coordinates to equatorial.
"""

__all__ = ["EquatorialCoordinates", ]
__author__ = "M.J. Huston"
__date__ = "2025-06-30"

import pandas
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

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        Perform the magnitude conversions and returns the modified DataFrame.
        """

        # Get ra, dec from l, b & add to dataframe
        l_arr, b_arr = dataframe['l'].to_numpy(), dataframe['b'].to_numpy()
        ra_arr, dec_arr = lb_to_ad(l_arr,b_arr)
        l_idx = np.where(dataframe.columns=='l')[0][0]
        dataframe.insert(l_idx, 'ra', ra_arr)
        dataframe.insert(l_idx+1, 'dec', dec_arr)

        # Get mura, mudec from coords and kinemaitcs & add to dataframe
        dist_arr = dataframe['Dist'].to_numpy()
        u_arr, v_arr, w_arr = dataframe['U'].to_numpy(), dataframe['V'].to_numpy(), dataframe['W'].to_numpy()
        vr_arr, mura_arr, mudec_arr = uvw_to_vrmuad(l_arr,b_arr,dist_arr, u_arr,v_arr,w_arr)
        mul_idx = np.where(dataframe.columns=='mul')[0][0]
        dataframe.insert(mul_idx, 'mura', mura_arr)
        dataframe.insert(mul_idx+1, 'mudec', mudec_arr)

        # Remove galactic coordinates if selected
        if not self.keep_galactic:
            dataframe.drop(columns=['l','b'],inplace=True)
            dataframe.drop(columns=['mul','mub'],inplace=True)
        
        return dataframe
