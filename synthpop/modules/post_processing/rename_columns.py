"""
Post-processing module to replace certain output column names. This replaces the "col_names"
configuration line from <v1.1.0.
"""

__all__ = ["RenameColumns", ]
__author__ = "M.J. Huston"
__date__ = "2025-08-13"

import pandas
import numpy as np
from ._post_processing import PostProcessing

class RenameColumns(PostProcessing):
    """
    Post-processing module to rename columns
    
    Attributes
    ----------
    old_names : list
        list of columns to rename
    new_names : list
        list of new names in same order as old_names
    """

    def __init__(self, model, logger, old_names=None, new_names=None, **kwargs):
        super().__init__(model,logger, **kwargs)
        self.old_names = old_names
        self.new_names = new_names

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        Make the column name changes.
        """
        
        dataframe.rename(columns=dict(zip(self.old_names, self.new_names)), inplace=True)
        return dataframe
