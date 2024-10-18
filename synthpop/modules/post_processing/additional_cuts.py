"""
This file contains post-processing to make additional cuts to the catalog
"""
__all__ = ["AdditionalCuts", ]
__author__ = "M.J. Huston"
__date__ = "2024-05-14"
__license__ = "GPLv3"
__version__ = "1.0.0"

import pandas
import numpy as np
from ._post_processing import PostProcessing

class AdditionalCuts(PostProcessing):
    def __init__(self, model, logger, standard_cuts=None, difference_cuts=None, **kwargs):
        """
        Parameters:
            model:
                SynthPop main model object
            standard_cuts:
                list of additional cuts to make to the catalog on individual columns
                expected format:
                [['param_name1', 'cut_type1', cut_limit1],
                ['param_name2', 'cut_type2', cut_limit2]]
                'param_nameN' can be any column in the DataFrame
                'cut_typeN' must be 'min' or 'max'
                cutlimitN should be the int or float value to cut at
            difference_cuts:
                list of additional cuts to make to the catalog on the difference between
                columns (e.g. colors)
                expected format:
                [['param_name1a', 'param_name1b', 'cut_type1', cut_limit1],
                ['param_name2a', 'param_name2b', 'cut_type2', cut_limit2]]
                'param_nameNm' can be any column in the DataFrame
                    note: the second will be subtracted from the first
                'cut_typeN' must be 'min' or 'max'
                cutlimitN should be the int or float value to cut at
        """
        self.model = model
        self.standard_cuts = standard_cuts
        self.difference_cuts = difference_cuts

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        It must accept the pandas data frame and must return a pandas data frame

        Parameters
        ----------
        dataframe : dataframe
            original SynthPop output as pandas data frame

        Returns
        -------
        dataframe : dataframe
            modified pandas data frame
        """
        dataframe = dataframe.reset_index()
        if not self.standard_cuts==None:
            for cut in self.standard_cuts:
                if cut[1]=='min':
                    idxs = np.where(dataframe[cut[0]]<cut[2])[0]
                elif cut[1]=='max':
                    idxs = np.where(dataframe[cut[0]]>cut[2])[0]
                else:
                    print("INVALID CUT TYPE")
                if len(idxs)>0:
                    dataframe = dataframe.drop(idxs)
        if not self.difference_cuts==None:
            for cut in self.difference_cuts:
                if cut[2]=='min':
                    idxs = np.where(dataframe[cut[0]]-dataframe[cut[1]]<cut[3])[0]
                elif cut[2]=='max':
                    idxs = np.where(dataframe[cut[0]]-dataframe[cut[1]]>cut[3])[0]
                else:
                    print("INVALID CUT TYPE")
                if len(idxs)>0:
                    dataframe = dataframe.drop(idxs)
        return dataframe
