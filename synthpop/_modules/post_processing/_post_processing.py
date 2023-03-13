"""
This file contains the base class for the Post processing.
"""
__all__ = ["PostProcessing", ]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__license__ = "GPLv3"
__date__ = "2022-11-16"
__version__ = '1.0.0'

import pandas


class PostProcessing:
    def __init__(self, model, logger,  **kwargs):
        """
        Parameters:
            model: SynthPop
                to access properties from the model.
                e.g: self.model.const, for all constants
                    self.model.params, for all input parameters
                    self.model.populations, for each population class and properties within
                    self.model.l_deg, self.model.b_deg,
                    self.model.filename_base for the directory
                            and filename of the output file (without extension)
                    ...
            logger: SynthPopLogger
                 can be used to add messages to the logging.

        kwargs:
            keyword arguments from the configuration
                ("advanced_post_processing":{"name":"name_of_the_subclass", kwargs}
        """
        self.model = model
        self.logger = logger

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        This is a placeholder for the postprocessing
        replace it whit what ever you think is useful
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
        return dataframe
