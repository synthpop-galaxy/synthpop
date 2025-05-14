"""
Base class for post-processing.
"""
__all__ = ["PostProcessing", ]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-11-16"

import pandas
from types import ModuleType
from .. import default_sun

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

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        This is a placeholder for the postprocessing.
        It must accept the pandas data frame and must return a pandas data frame.

        Parameters
        ----------
        dataframe : dataframe
            original SynthPop output as pandas data frame

        Returns
        -------
        dataframe : dataframe
            modified pandas data frame
        """
        return dataframe
