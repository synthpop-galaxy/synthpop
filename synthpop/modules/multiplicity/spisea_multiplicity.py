"""
Spisea multiplicity class holder
"""

__all__ = ["SpiseaMultiplicity"]
__author__ = "M.J. Huston"
__credits__ = ["M. Newman, M.J. Huston"]
__date__ = "2025-10-15"

import pandas as pd
import numpy as np
from ._multiplicity import Multiplicity
import pdb
import synthpop.constants as const
from spisea.imf import multiplicity as spisea_multiplicity

try:
    from constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR, DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)
except (ImportError, ValueError):
    from synthpop.constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR, DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)

class SpiseaMultiplicity(Multiplicity):
    def __init__(self, spisea_multiplicity_name="MultiplicityResolvedDK", spisea_multiplicity_kwargs={}, **kwargs):
        """
        Hi
        """
        super().__init__(**kwargs)
        self.name='SpiseaMultiplicity'
        self.spisea_multiplicity = getattr(spisea_multiplicity, spisea_multiplicity_name)(**spisea_multiplicity_kwargs)

    def generate_companions(self, pri_masses):
        """
        not valid here
        """
        raise ValueError("SpiseaMultiplicity can only be used by SpiseaGenerator") 
        
