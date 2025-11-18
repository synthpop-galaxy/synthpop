"""
Add some derived and additional parameters to the multi-star systems
"""

__all__ = ["CombineTables", ]
__author__ = "M. Newman, M.J. Huston"
__date__ = "2025-10-15"

import pandas as pd
import numpy as np
from ._post_processing import PostProcessing
from scipy.stats import maxwell, uniform_direction
from synthpop.synthpop_utils.coordinates_transformation import CoordTrans
import pdb

class CombineTables(PostProcessing):
    """
    Post-processing to add some derived columns for binaries.

    Attributes
    ----------
    
    """

    def __init__(self, model, logger, **kwargs):
        super().__init__(model, logger, **kwargs)

    def do_post_processing(self, system_df: pd.DataFrame
            companion_df: pd.DataFrame) -> (pd.DataFrame pd.DataFrame):
        """
        Perform the post-processing and return the modified DataFrame.
        """
        
        raise ValueError("THIS DOESNT WORK YET")

        # Catch case of zero companions
        # if not np.any(system_df['Is_Binary'].to_numpy()>0):
        #     return system_df
        in_system = [np.where(system_df['primary_ID']==pid) for pid in system_df['primary_ID']]
        system_df.loc[:,"system_mass"] = system_df['Mass']
        system_df.loc[:,"system_logL"] = system_df['log_L']
        for i in system_df[system_df['Is_Binary']>0].index:
            system_df.loc[i, "system_mass"] = np.sum(system_df['Mass'][in_system[i][0]])
            system_df.loc[i, "system_logL"] = np.log10(np.nansum(10**system_df['log_L'][in_system[i][0]]))
            if system_df['Is_Binary'][i]>1:
                system_df.loc[i, "q"] = system_df['Mass'][i]/system_df['Mass'][
                                            system_df['primary_ID'][i]]
            
        return system_df, None
