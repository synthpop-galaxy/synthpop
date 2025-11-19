"""
Postprocessing module to calculate blended magnitude of all stars within a given blend radius.
Note: consider using no regular magnitude cut, then running additional_cuts after this module
    to cut on blend magnitude.
"""

__all__ = ["Blending", ]
__author__ = "M.J. Huston"
__date__ = "2025-05-08"

import pandas as pd 
import numpy as np
from ._post_processing import PostProcessing
from scipy.spatial import KDTree
try:
    from synthpop_utils.utils_functions import add_magnitudes, combine_system_mags
except:
    from .synthpop_utils.utils_functions import add_magnitudes, combine_system_mags

class Blending(PostProcessing):
    """
    Postprocessing module to calculate blended magnitudes in selected filters
    for a given blend radius. To use different radii for different filters,
    include multiple instances of this postproc in your configuration.
    
    Attributes
    ----------
    blend_radius : float
        blending radius in arcseconds
    filters : list of str
        names of filters to compute blended magnitude in;
        must be included in your config's chosen_bands
    """

    def __init__(self, model, blend_radius, filters, logger, **kwargs):
        super().__init__(model,logger, **kwargs)
        self.blend_radius = blend_radius
        self.filters = filters

    def get_obs_mags(self, dataframe, filt):
        mags_noext = np.array(dataframe[filt]) + 5*np.log10(np.array(dataframe.Dist)*100)
        exts = np.nan*np.ones(len(mags_noext))
        for pop in np.unique(dataframe['pop']).astype(int):
            in_pop = np.where(np.array(dataframe['pop'])==pop)[0]
            ipop = int(pop)
            ext_from_map = np.array(dataframe[self.model.populations[ipop].extinction.A_or_E_type])
            ext_at_filt = self.model.populations[ipop].extinction.extinction_at_lambda(self.model.parms.eff_wavelengths[filt], ext_from_map)
            exts[in_pop] = ext_at_filt[in_pop]
        return mags_noext+exts

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        
        delta_l = np.array(system_df['l'] - self.model.l_deg) * np.cos(self.model.b_deg*np.pi/180.0) * 3600
        delta_b = np.array(system_df['b'] - self.model.b_deg) * 3600
        pts = np.transpose([delta_l, delta_b])
        kdt = KDTree(pts)
        res = kdt.query_ball_point(pts, self.blend_radius)
        if (not self.model.parms.combine_system_mags) and (companion_df is not None):
            system_df_new = combine_system_mags(system_df.copy(), companion_df, self.model.populations[0].bands)
        else:
            system_df_new = system_df
        for filt in self.filters:
            if self.model.parms.obsmag:
                mags = np.array(system_df_new[filt])
            else:
                mags = self.get_obs_mags(system_df_new, filt)

            mag_arr = np.array(list(map(lambda i: add_magnitudes(mags[i]), res)))
            system_df[filt+'_bl'] = mag_arr
         
        return system_df, companion_df
