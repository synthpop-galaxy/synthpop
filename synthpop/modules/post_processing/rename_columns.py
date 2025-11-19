"""
Post-processing module to replace certain output column names. This replaces the "col_names"
configuration line from <v1.1.0.
"""

__all__ = ["RenameColumns", ]
__author__ = "M.J. Huston"
__date__ = "2025-08-13"

import pandas as pd
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

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Make the column name changes.
        """
        system_df.rename(columns=dict(zip(self.old_names, self.new_names)), inplace=True)
        if companion_df is not None:
            for i, column in enumerate(self.old_names):
                if column in companion_df.columns:
                    companion_df.rename(columns={column: self.new_names[i]}, inplace=True)
        return system_df, companion_df
