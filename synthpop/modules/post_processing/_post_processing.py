"""
Base class for post-processing.
"""
__all__ = ["PostProcessing", ]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-11-16"

import pandas as pd
from types import ModuleType
from .. import default_sun
from typing import Tuple

class PostProcessing:
    """
    Base class for post-processing, for catalog modifications after the
    regular SynthPop generation process has completed.
    """
    
    def __init__(self,
            model: ModuleType = None,
            logger: ModuleType = None,
            sun: ModuleType = None,
            **kwargs):

        #: SynthPop model object
        #:        to access properties from the model
        #:        e.g: self.model.const, for all constants
        #:             self.model.params, for all input parameters
        #:             self.model.populations, for each population class and properties within
        #:             self.model.l_deg, self.model.b_deg,
        #:             self.model.filename_base for the directory
        #:                and filename of the output file (without extension)
        #:             ...
        self.model = model
        #: SynthPop logger
        self.logger = logger
        #: Solar and LSR parameters
        self.sun = sun if sun is not None else default_sun

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        This is a placeholder for the postprocessing.
        It must accept the pandas data frame and must return a pandas data frame.

        Parameters
        ----------
        system_df : pandas dataframe
            original SynthPop output for star systems
            (may be all individual stars if no multiplicity)
        companion_df : pandas dataframe
            original SynthPop output for companions
            (may be None for no multiplicity)

        Returns
        -------
        system_df : pandas dataframe
            modified star systems table
        companion_df : pandas dataframe
            modified companions table
        """
        return system_df, companion_df
