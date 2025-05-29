"""
Evolution module to store information for the SpiseaGenerator. Not valid for a standard StarGenerator.
"""

__all__ = ["Spisea", ]
__author__ = "M.J. Huston"
__date__ = "2024-11-20"
__license__ = "GPLv3"
__version__ = "1.0.0"

from ._evolution import CombineEvolution

class SpiseaCluster(CombineEvolution):
    """
    Placeholder object to store settings for use by the SpiseaGenerator, which 
    will generate and evolve stars as binned SPISEA clusters.
    """
    def __init__(self, evo_model_name="MISTv1.2", atm_func_name="get_merged_atmosphere", 
                    wd_atm_func_name="get_wd_atmosphere", ifmr_name=None, 
                    multiplicity=None, **kwargs):
        self.evo_model_name = evo_model_name
        self.atm_func_name = atm_func_name
        self.wd_atm_func_name = wd_atm_func_name
        self.ifmr_name=ifmr_name
        self.multiplicity = multiplicity