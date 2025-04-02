"""
Post-processing subclass that saves the result from all locations in a combined csv file.
"""

__all__ = ["CombinedCsv", ]
__author__ = "J. Klüter"
__date__ = "2023-01-23"

import os
import pandas
from ._post_processing import PostProcessing

class CombinedCsv(PostProcessing):
    """
    Post-processing subclass that saves the result from all locations in a combined csv file.
    This is in addition to the single location files.

    NOTE that the file is overwritten every time the class is initialized.
    """

    def __init__(self, model, logger, combined_filename=None, **kwargs):
        super().__init__(model, logger, **kwargs)
        if combined_filename is None:
            self.combined_filename = os.path.join(model.parms.output_location,
                                f"{model.parms.model_name}.combined.csv")
        else:
            #: File name for combined output
            self.combined_filename = combined_filename
        if os.path.isfile(self.combined_filename):
            os.remove(self.combined_filename)

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        Combine all catalogs into one output csv file, and returns the unchanged DataFrame.
        """
        # check if the file exist, if so it will not write a new header
        file_exist = os.path.isfile(self.combined_filename)
        # convert dataframe to csv sting
        csv_data = dataframe.to_csv(index=False, header=file_exist)

        # open the file in append mode
        self.logger.info(f"attach {dataframe.shape[0]} to {self.combined_filename}")
        with open(self.combined_filename, "a") as f:
            # write the data to the file
            f.write(csv_data)

        return dataframe
