"""
This file contains post-processing to account for binary populations, based on Raghavan et al. 2010

"""

__author__ = "Marz Newman"

import pandas as pd
import numpy as np
import math
import glob
import os
from ._post_processing import PostProcessing
import sys
import synthpop.constants as const
from synthpop.star_generator import StarGenerator
from synthpop.synthpop_utils.synthpop_logging import logger
from synthpop.synthpop_utils import Parameters
import synthpop.synthpop_utils as sp_utils
from synthpop.population import Population
from synthpop.modules.initial_mass_function import InitialMassFunction
from synthpop.modules.age import Age
from synthpop.modules.metallicity import Metallicity
try:
	from constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR, DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)
except (ImportError, ValueError):
	from synthpop.constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR, DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)

class Binary(PostProcessing):
	def __init__(self, **kwargs):
		"""
		Parameters:
			model:
				SynthPop main model object
		"""

		super().__init__(**kwargs)

	@staticmethod
	def check_is_binary(prim_masses):
		"""
		Probabilistic determination of binary star status based on temperature bins and probabilities.

		Parameters
		----------
		temperatures
			array of float values for star massas.

		Returns
		-------
		is_binary
			array of boolean values indicating whether each star is a binary star.
		"""
		
		is_binary = []
		num_isnan = 0

		# Function from Duchene 2013
		def multiplicity_fraction(prim_mass):
			# Keep fraction at 20% for <0.1 M_Sun
			if mass < 0.1:
				return 0.2
			# Everything else
			else:
				return 38.36*(prim_mass)**(0.27) / 100	# Divide by 100 to go from percent to probability
		
		binary_frac = []
		num_fraction_gt_1 = 0
		print("\n+++++++++++++++++++")
		for mass in prim_masses:
			rand_num = np.random.rand()
			fraction = multiplicity_fraction(mass)
			binary_status = rand_num < fraction
			if fraction >= 1:
				#print("Multiplicity fraction is greater than 1")
				num_fraction_gt_1 += 1
				print(fraction)
				print(binary_status)
			is_binary.append(binary_status)
			binary_frac.append(fraction)
		
		print(f"Number of multiplicity fractions > 1: {num_fraction_gt_1}")

		binary_count = sum(is_binary)

		print(f"Number of identified binary stars: {binary_count}")
		print(f"Length of binary_frac list: {len(binary_frac)}")

		return (is_binary, binary_frac)

	def draw_companion_m_ratio(self):
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

		# Plug in random number and pull a mass ratio from the toy probability function inspired by Figure 16 of D. Raghavan 2010
		def pull_from_distribution(random_number):
			# Probability function for the first, linear section
			if random_number < probability_bins[0]:
				P = np.random.rand()  # Uniform random number for inversion method
				#mass_ratio = np.sqrt((2 * P) / 27.5)
				mass_ratio = np.sqrt((2 * P) / 50)
			elif random_number < probability_bins[1]:
				mass_ratio = np.random.uniform(0.2, 0.95)
			else:
				mass_ratio = np.random.uniform(0.95, 1.)
			return mass_ratio

		random_number = np.random.rand()

		# Pull random mass ratio
		random_mass_ratio = pull_from_distribution(random_number)

		return random_mass_ratio

	def draw_period(self):
		"""
		Find binary orbital period from Gaussian distribution in Figure 13 of Raghavan et al. 2010

		Returns
		-------
		logP
			float value of the logarithm of period
		"""
		# Draw a log(Period) from the Raghavan Figure 13 Gaussian
		mu = 5.03	# Value from Raghavan Figure 13
		sigma = 2.28	# Value from Raghavan Figure 13
		logP = np.random.normal(loc = mu, scale = sigma)

		# Keep drawing if we get a period less than 1 day (logP < 0)
		while logP <= 0 or logP >= 9:
			logP = np.random.normal(loc=mu, scale=sigma)

		return logP

	def do_post_processing(self, dataframe: pd.DataFrame) -> pd.DataFrame:
		"""
		Evolves binary companion to obtain properties

		Parameters
		----------
		dataframe : dataframe
		original SynthPop output as pandas data frame

		Returns
		-------
		dataframe : dataframe
		modified pandas data frame
		"""
		print('\nRunning binary post processing')
				
		num_samples = len(dataframe)
		
		# Initialize ID column
		# Add a column for multiplicity fraction
		dataframe['ID'] = np.nan

		# Make a new column for identifying binaries [0=single, 1=primary, 2=secondary]
		dataframe['Is_Binary'] = 0
		
		# Make a new column that will point secondaries to the primary id (0 by default since initial dataframe stars are single)
		dataframe['primary_ID'] = np.nan

		# Identify primary stars and flag them in the new 'Is_Binary' column
		(binary_flags, binary_frac) = Binary.check_is_binary(dataframe['iMass'])
		
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
		
		# Which populations do the primaries belong to?
		used_pop_ids = dataframe[dataframe['Is_Binary'] == 1]['pop'].unique()
		
		# Load config file parameters
		default_config = DEFAULT_CONFIG_FILE
		model_dir = DEFAULT_MODEL_DIR+"/"+self.model.parms.model_name
		
		args = sp_utils.parser()
		specific_config = os.path.join(DEFAULT_CONFIG_DIR, args.specific_config)

		# Parameters for making an instance of Populations, taken from the config files and model directory
		glbl_params = sp_utils.Parameters(specific_config, default_config, model_dir)
		max_mass = self.model.parms.mass_lims['max_mass']
		
		# Collect all the dataframe column names - all of these are the same for every population, so I just used popid = 0
		headers = const.COL_NAMES + self.model.populations[0].glbl_params.opt_iso_props + self.model.populations[0].bands
		extinction_index = headers.index("ExtinctionInMap")
		
		# Replace "ExtinctionInMap" with the output of the extinction map - all of these are the same for every population, so I just used popid = 0
		if 'ExtinctionInMap' in headers:
			headers[extinction_index] = self.model.populations[0].extinction.A_or_E_type
			
		# Make a copy of the dataframe for secondaries
		secondary_df = dataframe.copy(deep=True)
		
		mass_ratio_df = pd.DataFrame()
		periods_list = []
				
		# Loop over populations 
		for i in used_pop_ids:
			i = int(i)
			popid = self.model.populations[i].popid
			(population_density, imf, age, metallicity, kinematics, evolution, extinction) = self.model.populations[popid].assign_subclasses()

			# Requested properties
			bands = list(evolution.bands)
			props = set(const.REQ_ISO_PROPS + glbl_params.opt_iso_props + bands)			
			
			# Mask out only one population
			population_mask = (dataframe['Is_Binary'] == 1) & (dataframe["pop"] == popid)
			population_df = dataframe[population_mask]
			original_indices = population_df.index
					
			# Create an instance of the Binary subclass for companions
			companions = [Binary() for j in range(len(population_df))]
			
			# Draw an initial mass	
			mass_ratios = pd.Series([companion.draw_companion_m_ratio() for companion in companions], index = original_indices)
			mass_ratio_df = pd.concat([mass_ratio_df, mass_ratios])#, ignore_index=True)
			prim_masses = population_df["iMass"]
			mass = pd.Series(mass_ratios.to_numpy() * prim_masses.to_numpy(), index=prim_masses.index)	# Array of companion initial masses
			
			#print(mass)
			
			# Run get_evolved_props
			ref_mag, s_props, final_phase_flag, inside_grid, not_evolved = self.model.populations[popid].generator.get_evolved_props(mass.to_numpy(), population_df['Fe/H_initial'].to_numpy(), population_df['age'].to_numpy(), props)
			
			# Switch back to original indices
			ref_mag = pd.Series(ref_mag, index=original_indices)
			s_props = pd.DataFrame(s_props, index=original_indices)
			inside_grid = pd.Series(inside_grid, index=original_indices)
			not_evolved = pd.Series(not_evolved, index=original_indices)

			# Run extract_properties
			m_evolved, props, user_props = self.model.populations[popid].extract_properties(mass, s_props, const.REQ_ISO_PROPS, glbl_params.opt_iso_props, inside_grid, not_evolved)
			
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
			periods = [companion.draw_period() for companion in companions]
			periods_list.extend(periods)
		
		############# Add new df information to final catalog dataframe #############
		
		# Insert secondaries beneath their corresponding primaries in the full dataframe
		combined_list = []
		for i in range(len(dataframe)):
			combined_list.append(dataframe.iloc[i])
			if binary_flags[i]:
				combined_list.append(secondary_df.iloc[i])
		combined_df = pd.DataFrame(combined_list).reset_index(drop=True)
			
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
