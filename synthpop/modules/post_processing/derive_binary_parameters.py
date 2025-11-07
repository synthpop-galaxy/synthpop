"""
Add some derived and additional parameters to the multi-star systems
"""

__all__ = ["DeriveBinaryParameters", ]
__author__ = "M. Newman, M.J. Huston"
__date__ = "2025-10-15"

import pandas
import numpy as np
from ._post_processing import PostProcessing
from scipy.stats import maxwell, uniform_direction
from synthpop.synthpop_utils.coordinates_transformation import CoordTrans
import pdb

class DeriveBinaryParameters(PostProcessing):
    """
    Post-processing to add some derived columns for binaries.

    Attributes
    ----------
    
    """

    def __init__(self, model, logger, **kwargs):
        super().__init__(model, logger, **kwargs)

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        Perform the post-processing and return the modified DataFrame.
        """

        # Catch case of zero companions
        # if not np.any(dataframe['Is_Binary'].to_numpy()>0):
        #     return dataframe
        in_system = [np.where(dataframe['primary_ID']==pid) for pid in dataframe['primary_ID']]
        dataframe.loc[:,"system_mass"] = dataframe['Mass']
        dataframe.loc[:,"system_logL"] = dataframe['log_L']
        for i in dataframe[dataframe['Is_Binary']>0].index:
            dataframe.loc[i, "system_mass"] = np.sum(dataframe['Mass'][in_system[i][0]])
            dataframe.loc[i, "system_logL"] = np.log10(np.nansum(10**dataframe['log_L'][in_system[i][0]]))
            if dataframe['Is_Binary'][i]>1:
                dataframe.loc[i, "q"] = dataframe['Mass'][i]/dataframe['Mass'][
                                            dataframe['primary_ID'][i]]
            
        return dataframe
