"""
Filler class for zero multiplicity (single stars only)
"""

__all__ = ["NoMultiplicity"]
__author__ = "M.J. Huston"
__credits__ = ["M.J. Huston"]
__date__ = "2025-10-15"

import pandas as pd
import numpy as np
from ._multiplicity import Multiplicity
import pdb
import synthpop.constants as const

try:
	from constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR, DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)
except (ImportError, ValueError):
	from synthpop.constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR, DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)

class NoMultiplicity(Mutiplicity):
	def __init__(self, **kwargs):
		"""
		Hi
		"""
		super().__init__(**kwargs)
		self.name='NoMultiplicity'

	def generate_companions(self, pri_masses):
		"""
		Generates companion stars

		Parameters
		----------
		pri_masses : ndarray
			primary star masses

		Returns
		-------
		companions_table : dataframe
			companion properties
		"""
					
		return np.array([]), np.array([]), np.array([])
