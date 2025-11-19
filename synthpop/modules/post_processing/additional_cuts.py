"""
Post-processing module to make additional cuts to the catalog.
"""

__all__ = ["AdditionalCuts", ]
__author__ = "M.J. Huston"
__date__ = "2024-05-14"

import pandas as pd
import numpy as np
from ._post_processing import PostProcessing

class AdditionalCuts(PostProcessing):
    """
    Post-processing module to make additional cuts on which stars are included in the catalog.
    These may be based on a minimum or maximum value for any column(s),
    or a minimum or maximum value of the difference between any column pair(s).
    
    Attributes
    ----------
    standard_cuts : list
        list of cuts to make by data column name, cut type, and cut limit;
        e.g. [['param_name1', 'cut_type1', cut_limit1],
        ['param_name2', 'cut_type2', cut_limit2]],
        where 'param_nameN' can be any column in the DataFrame,
        'cut_typeN' must be 'min' or 'max', and
        cutlimitN should be the int or float value to cut at.
    difference_cuts : list
        list of cuts to make based on column value differences (e.g. color);
        e.g. [['param_name1a', 'param_name1b', 'cut_type1', cut_limit1],
        ['param_name2a', 'param_name2b', 'cut_type2', cut_limit2]],
        where 'param_nameNm' can be any column in the DataFrame,
        'cut_typeN' must be 'min' or 'max', and
        cutlimitN should be the int or float value to cut at.
    """

    def __init__(self, model, logger, standard_cuts=None, difference_cuts=None, **kwargs):
        super().__init__(model,logger, **kwargs)
        self.standard_cuts = standard_cuts
        self.difference_cuts = difference_cuts

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Make the cuts based on the system table parameters
        and applies to both tables.
        """
        
        if not self.standard_cuts==None:
            for cut in self.standard_cuts:
                if cut[1]=='min':
                    system_df = system_df[system_df[cut[0]]>cut[2]]
                elif cut[1]=='max':
                    system_df = system_df[system_df[cut[0]]<cut[2]]
                else:
                    raise ValueError('Invalid cut type: '+cut[1]+'. Valid types are min and max')
        if not self.difference_cuts==None:
            for cut in self.difference_cuts:
                if cut[2]=='min':
                    system_df = system_df[(system_df[cut[0]]-system_df[cut[1]])>cut[3]]
                elif cut[2]=='max':
                    system_df = system_df[(system_df[cut[0]]-system_df[cut[1]])<cut[3]]
                else:
                    raise ValueError('Invalid cut type: '+cut[2]+'. Valid types are min and max')
        if companion_df is not None:
            companion_df = companion_df[np.isin(companion_df['system_idx'].to_numpy(),
                                        system_df['system_idx'].to_numpy)]
        return system_df, companion_df
