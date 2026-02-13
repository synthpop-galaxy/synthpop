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
    Post-processing to mimic old file format with binaries.

    Attributes
    ----------
    
    """

    def __init__(self, model, logger, **kwargs):
        super().__init__(model, logger, **kwargs)

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Perform the post-processing and return the modified DataFrame.
        """
        
        if np.any(system_df.n_companions>1):
            raise NotImplementedError("THIS DOESNT WORK YET for higher order than binary systems")
        # Do some renames
        system_df.rename(columns={'system_Mass':'total_mass', 'system_idx':'ID'}, inplace=True)
        companion_df.rename(columns={'system_Mass': 'total_mass', 'system_idx':'primary_ID'}, inplace=True)
        companion_df.loc[:,'logP'] = np.log10(companion_df['period'])
        companion_df.loc[:,'ID'] = np.max(system_df['ID']) + np.arange(len(companion_df))
        # Set up combo columns
        system_df.loc[:,'combined_logL'] = np.nan
        system_df.loc[:,'combined_logP'] = np.nan
        system_df.loc[:,'q'] = np.nan
        system_df.loc[:,'Is_Binary'] = 0
        # Match up companions and primaries for combined quantities
        system_df.set_index('ID', inplace=True)
        companion_system_idxs = companion_df['primary_ID'].to_numpy()
        companion_df.loc[:,'q'] = companion_df['Mass'] / system_df['Mass'][companion_system_idxs].to_numpy()
        combined_logL = np.log10(10**companion_df['log_L'] + 10**system_df['log_L'][companion_system_idxs].to_numpy())
        companion_df.loc[:,'combined_logL'] = combined_logL
        companion_df.loc[:, 'Is_Binary'] = 2
        # Update primaries combined values
        system_df.loc[companion_system_idxs, 'combined_logL'] = combined_logL.to_numpy()
        system_df.loc[companion_system_idxs, 'combined_logP'] = companion_df['logP'].to_numpy()
        system_df.loc[companion_system_idxs, 'Is_Binary'] = 1
        system_df.reset_index(inplace=True)
        # Combine table
        system_df = pd.concat([system_df, companion_df])
        system_df.sort_values(['primary_ID', 'Is_Binary'])
            
        return system_df, None
