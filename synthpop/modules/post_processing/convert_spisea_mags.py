"""
Postprocessing module to convert magnitude systems for any
filters provided by SPISEA.
"""

__all__ = ["ConvertSpiseaMags", ]
__author__ = "M.J. Huston"
__date__ = "2026-02-21"

import pandas as pd 
import numpy as np
from ._post_processing import PostProcessing
from spisea import synthetic

class ConvertSpiseaMags(PostProcessing):
    """
    Postprocessing module to convert magnitude systems from Vega for any
    filters provided by SPISEA.
    
    Attributes
    ----------
    system : str
        magnitude system (options: AB)
    """

    def __init__(self, model, logger, system='AB', **kwargs):
        super().__init__(model,logger, **kwargs)
        if system.upper() not in ['AB', 'ST']:
            raise NotImplementedError(f"System {system} not availabel in SPISEA.")
        self.system = system.upper()

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Perform the magnitude conversions and returns the modified DataFrames.
        """

        for f in self.model.populations[0].bands:
            f_str = synthetic.get_obs_str(f)
            if self.system=='AB':
                conv = synthetic.calc_ab_vega_filter_conversion(f_str)
            elif self.system=='ST':
                conv = synthetic.calc_st_vega_filter_conversion(f_str)
            system_df.loc[:,f] += conv
            if companion_df is not None:
                companion_df.loc[:,f] += conv 
        
        return system_df, companion_df
