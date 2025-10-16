"""
Binary companion generator, based on Raghavan et al. 2010

NOT WORKING YET
"""

__all__ = ["Multiplicity"]
__author__ = "M. Newman, M.J. Huston"
__credits__ = ["M. Newman, M.J. Huston"]
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

class Raghavan(Mutiplicity):
	def __init__(self, **kwargs):
		"""
		Hi
		"""
		super().__init__(**kwargs)

	@staticmethod
	def check_is_binary(pri_masses):
		"""
		Probabilistic determination of binary star status based on temperature bins and probabilities.

		Parameters
		----------
		pri_masses
			masses of the primaries

		Returns
		-------
		is_binary
			array of boolean values indicating whether each star is a binary
		"""
		
		# Function from Duchene 2013
		binary_frac = (pri_masses<=0.1)*0.2 + \
				   (pri_masses>0.1)*0.3836*pri_masses**0.27

		is_binary = np.random.rand(len(pri_masses))>binary_frac

		return (is_binary, binary_frac)

	@staticmethod
	def draw_companion_m_ratio(n):
		"""
		Mass ratio (M2/M1) calculation from toy probability function based on Figure 16 of Raghavan et. al 2010

		Returns
		-------
		random_mass_ratio
			float value for the mass ratio (M2/M1) of the binary system
		"""
		# Normalization factor (maximum value of N on Figure 16 of Raghavan 2010)
		#normalization_factor = 13

		# Separate systems into 3 different sections based on Raghavan 2010, Figure 16
		probability_bins = [0.103743, 0.778075+0.103743, 1.]	# Add previous bins to get *cumulative* values
		random_numbers = np.random.rand(n)
		bin_nos = np.search_sorted(probability_bins, random_numbers)

		random_mass_ratio = np.zeros(n)
		random_mass_ratio[bin_nos==0] = np.sqrt(2*np.random.rand(len(np.where(bin_nos==0)[0]))/50)
		random_mass_ratio[bin_nos==1] = np.random.uniform(0.2, 0.95, len(np.where(bin_nos==1)[0]))
		random_mass_ratio[bin_nos==2] = np.random.uniform(0.95, 1.0, len(np.where(bin_nos==2)[0]))

		return random_mass_ratio

	def draw_period(n):
		"""
		Find binary orbital period from Gaussian distribution in Figure 13 of Raghavan et al. 2010

		Returns
		-------
		logP
			float value of the logarithm of period
		"""
		# Draw a log(Period) from the Raghavan Figure 13 Gaussian
		logP = np.repeat(-1, n)
		mu = 5.03	# Value from Raghavan Figure 13
		sigma = 2.28	# Value from Raghavan Figure 13
		logP = np.random.normal(loc = mu, scale = sigma, size=n)

		# Keep drawing if we get a period less than 1 day (logP < 0)
		while np.any(logP <= 0) or np.any(logP >= 9):
			redraw_idx = ((logP <= 0) | (logP >= 9))
			logP[redraw_idx] = np.random.normal(loc=mu, scale=sigma, size=len(np.where(redraw_idx)[0]))

		return logP

	def generate_companions(self, pri_masses, opt_iso_props, bands):
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
				
		num_samples = len(dataframe)
		
		# Initialize ID column
		# Add a column for multiplicity fraction
		dataframe['ID'] = np.nan

		# Make a new column for identifying binaries [0=single, 1=primary, 2=secondary]
		dataframe['Is_Binary'] = 0
		
		# Make a new column that will point secondaries to the primary id (0 by default since initial dataframe stars are single)
		dataframe['primary_ID'] = np.nan

		# Identify primary stars and flag them in the new 'Is_Binary' column
		(binary_flags, binary_frac) = self.check_is_binary(dataframe['iMass'])
		
		# Add a column for multiplicity fraction
		dataframe['binary_frac'] = binary_frac

		for i, is_binary in enumerate(binary_flags):
			if is_binary:
				# Set the Is_Binary column to True for the identified binary star
				dataframe.at[i, 'Is_Binary'] = 1
		
		num_binaries = sum(binary_flags)
		
		# Initialize kinematic paramaters for the whole system (used at the very end)
		proper_motions = np.full((num_binaries, 3), np.nan)
		velocities = np.full((num_binaries, 3), np.nan)
		vr_lsr = np.repeat(np.nan, num_binaries)
			
		# Make a copy of the dataframe for secondaries
		secondary_df = dataframe.copy(deep=True)
		
		mass_ratio_df = pd.DataFrame()
		periods_list = []
			
		(population_density, imf, age, metallicity, kinematics, evolution, extinction) = self.model.populations[popid].assign_subclasses()

		# Requested properties
		props = set(const.REQ_ISO_PROPS + opt_iso_props + bands)			
		
		# Mask out only one population
		population_mask = (dataframe['Is_Binary'] == 1) & (dataframe["pop"] == popid)
		population_df = dataframe[population_mask]
		original_indices = population_df.index
				
		# Draw an initial mass	
		mass_ratios = pd.Series([self.draw_companion_m_ratio() for j in range(len(population_df))], index = original_indices)
		mass_ratio_df = pd.concat([mass_ratio_df, mass_ratios])#, ignore_index=True)
		pri_masses = population_df["iMass"]
		mass = pd.Series(mass_ratios.to_numpy() * pri_masses.to_numpy(), index=pri_masses.index)	# Array of companion initial masses
		
		
		# Run get_evolved_props
		ref_mag, s_props, final_phase_flag, inside_grid, not_evolved = self.model.populations[popid].generator.get_evolved_props(mass.to_numpy(), population_df['Fe/H_initial'].to_numpy(), population_df['age'].to_numpy(), props)
		
		# Switch back to original indices
		ref_mag = pd.Series(ref_mag, index=original_indices)
		s_props = pd.DataFrame(s_props, index=original_indices)
		inside_grid = pd.Series(inside_grid, index=original_indices)
		not_evolved = pd.Series(not_evolved, index=original_indices)

		# Run extract_properties
		m_evolved, props, user_props = self.model.populations[popid].extract_properties(mass, s_props, const.REQ_ISO_PROPS, opt_iso_props, inside_grid, not_evolved)
		
		m_evolved = pd.Series(m_evolved, index = original_indices)
		props = pd.DataFrame(props, index = original_indices)
		user_props = pd.DataFrame(user_props, index = original_indices)
		
		position = population_df[['x', 'y', 'z', 'Dist', 'l', 'b']].to_numpy()
		
		# Get maximum distance and step_size
		max_distance = getattr(self.model.populations[popid].pop_params, 'max_distance', glbl_params.max_distance)
		radii = np.linspace(0, max_distance, len(population_df)+1)

		r_inner = radii[:-1]
				
		# Run extract_magnitudes
		mag, extinction_in_map = self.model.populations[popid].extract_magnitudes(r_inner, position[:, 3:6], ref_mag, s_props)
		mag = pd.DataFrame(mag, index=original_indices, columns=bands)
		extinction_in_map = pd.DataFrame(extinction_in_map, index=original_indices)
		
		# Update secondary_df
		secondary_df.loc[population_df.index, 'Is_Binary'] = 2
		secondary_df.loc[population_df.index, 'iMass'] = mass
		secondary_df.loc[population_df.index, 'Mass'] = m_evolved
		secondary_df.loc[population_df.index, 'A_Ks'] = extinction_in_map
		secondary_df.loc[population_df.index, 'log_L'] = s_props["log_L"]
		secondary_df.loc[population_df.index, 'log_g'] = s_props["log_g"]
		secondary_df.loc[population_df.index, 'log_Teff'] = s_props["log_Teff"]
		secondary_df.loc[population_df.index, 'log_R'] = s_props["log_R"]

		for col in mag.columns:
			secondary_df.loc[population_df.index, col] = mag[col]
		
		# Draw periods for binary stars
		periods = [self.draw_period() for j in range(len(population_df))]
		periods_list.extend(periods)
	
		############# Add new df information to final catalog dataframe #############
		
		# Insert secondaries beneath their corresponding primaries in the full dataframe
		combined_list = []
		for i in range(len(dataframe)):
			combined_list.append(dataframe.iloc[i])
			if binary_flags[i]:
				combined_list.append(secondary_df.iloc[i])
		combined_df = pd.DataFrame(combined_list).reset_index(drop=True)
		#pdb.set_trace()
		# Add a new column with a six-digit identification number
		combined_df['ID'] = np.array([f"{i:09d}" for i in range(1, len(combined_df) + 1)])

		# Add the primary IDs to "primary_ID" columns of the secondary
		for i in range(1, len(combined_df)):
			if combined_df.loc[i, "Is_Binary"] == 2.0:
				combined_df.loc[i, "primary_ID"] = combined_df.loc[i - 1, "ID"]

		# Initialize new combined property columns with NaN values
		combined_df['combined_logL'] = np.nan
		combined_df['total_mass'] = np.nan
		combined_df['q'] = np.nan
		combined_df['combined_logP'] = np.nan

		
		# Find the indices of primary stars in combined_df
		prim_indices = combined_df[combined_df['Is_Binary'] == 1].index
				
		for index in prim_indices:
			
			# Combine luminosities
			logL1 = combined_df.loc[index, 'log_L']
			logL2 = combined_df.loc[index+1, 'log_L']
			combined_lum = np.log10(10**logL1 + 10**logL2)
			
			# Combine mass
			mass1 = combined_df.loc[index, 'Mass']
			mass2 = combined_df.loc[index+1, 'Mass']
			combined_mass = mass1 + mass2
			
			# Add combined properties to the primary star and secondary star rows
			combined_df.loc[index, 'combined_logL'] = combined_lum
			combined_df.loc[index, 'total_mass'] = combined_mass
			combined_df.loc[index+1, 'combined_logL'] = combined_lum
			combined_df.loc[index+1, 'total_mass'] = combined_mass
			
		
		mass_ratio_df = mass_ratio_df.reset_index()
		
		for reset_index, row in mass_ratio_df.iterrows():
			original_index = prim_indices[reset_index]
			ratio = row[0]
			combined_df.loc[original_index, 'q'] = ratio
			combined_df.loc[original_index, 'combined_logP'] = periods_list[reset_index]
			combined_df.loc[original_index+1, 'q'] = ratio
			combined_df.loc[original_index+1, 'combined_logP'] = periods_list[reset_index]
		
					
		return combined_df