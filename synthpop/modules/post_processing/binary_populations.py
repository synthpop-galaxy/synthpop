"""
This file contains post-processing to account for binary populations, based on Raghavan et al. 2010

Not ready for use
"""

__author__ = "Marz Newman"

import pandas as pd
pd.set_option('display.max_rows', None)	# Debugging
import numpy as np
import math
import glob
import os
np.random.seed(1234)
from ._post_processing import PostProcessing
import sys
#sys.path.append('/home/marznewman/.local/lib/python3.8/site-packages/')
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


#print('Running post processing')

# Make a subclass of StarGenerator
class BinaryGenerator(StarGenerator):
	def __init__(self, imf_module, age_module, met_module, evolution_module, glbl_params, logger, binary_property):
		super().__init__(imf_module, age_module, met_module, evolution_module, glbl_params, logger)
		self.binary_property = binary_property

	def generate_binary_star(self):
		print("BinaryGenerator class works")

class Binary(PostProcessing):
	def __init__(self, **kwargs):
		"""
		Parameters:
			model:
				SynthPop main model object
		"""
		#self.imf=None
		#self.age=None
		#self.metallicity=None
		#self.evolution=None
		#(self.population_density, self.imf, self.age, self.metallicity, self.kinematics, self.evolution, self.extinction) = self.assign_subclasses()
		#self.generator = StarGenerator(self.imf, self.age, self.metallicity, self.evolution, self.glbl_params, logger)

		super().__init__(**kwargs)

	@staticmethod
	def temp_is_binary(temperatures):
		"""
		Probabilistic determination of binary star status based on temperature bins and probabilities.

		Parameters
		----------
		temperatures
			array of float values for star temperatures.

		Returns
		-------
		is_binary
			array of boolean values indicating whether each star is a binary star.
		"""
		#print("Length of temperatures: ",len(temperatures))
		# Temperature bins
		#T_bins = np.array([3850, 5300, 5920, 7240, 9500, 31000, 41000])
#		T_bins = np.array([3000, 3500, 4000, 4500, 5000, 5500, 6000])
		# probability for each type for each temperature bin
#		T_binary = np.array([0.05, 0.15, 0.3, 0.35, 0.1, 0.05, 0.0])
		# determine if each temperature is a binary star
		#print(temperatures)

		# Mamajek 2013 Table 5 temperature bins
		#T_bins = np.array([44900, 31400, 7220, 6350, 5770, 4830, 3850, 2570, 450])
		T_bins = np.array([2570, 3850, 4830, 5770, 6350, 7220, 31400, 44900])
#		spTyp_bins = np.array([0, 1, 2, 3, 4, 5, 6, 7])	# Spectral type indices: 0 for O-B, 1 for B-F, 2 for F0-F6, 3 for F6-G2, 4 for G2-K3, 5 for K3-M0, 6 for M0-M8, and 7 for M8-Y
		spTyp_bins = np.array([0, 1, 2, 3, 4, 5, 6, 7])	# Spectral type indices: 0 for M8-Y, 1 for M0-M8, 2 for K3-M0, 3 for G2-K3, 4 for F6-G2, 5 for F0-F6, 6 for B-F, and 7 for O-B
		percentage_bins = [20, 35, 37.5, 40, 50, 60, 70, 75]
		probability_bins = percentage_bins / np.sum(percentage_bins)
		#print(probability_bins)

		is_binary = []

		# Initialize a count for identified binary stars
		binary_count = 0

		num_isnan = 0


#		for temp in temperatures:
			# Sort star temperature into the correct temperature bin
			#print("Temperature: "+str(temp))
			#			if math.isnan(temp):
			#print(str(temp)+" Temperature is nan")
#				bin_index = 0
#				num_isnan = num_isnan+1
#			else:
#				bin_index = np.searchsorted(T_bins, temp)
			#print("Bin index: "+str(bin_index))
			# Get the probability for the star to be binary, given the temperature bin
			#probability = T_binary[bin_index - 1] if bin_index > 0 else 0
#			probability = T_binary[bin_index] if bin_index < len(T_binary) else 0
			#print("Probability: "+str(probability))
#			rand_num=np.random.rand()
			#print("Random number draw: "+str(rand_num))
#			binary_status = rand_num < probability
#			is_binary.append(binary_status)
			#print(rand_num < probability)

#		binary_count = sum(is_binary)

		for temp in temperatures:
			if math.isnan(temp):
				#print(str(temp)+" Temperature is nan")
				bin_index = 0
				num_isnan = num_isnan+1
				#print("NaN")
			else:
				bin_index = np.searchsorted(T_bins, temp)
				#print(bin_index)
			#print(bin_index, len(spTyp_bins))
			#probability = spTyp_bins[bin_index] if bin_index < len(spTyp_bins) else 0
			#print(temp, probability)
			#print(temp, bin_index)
			#print(bin_index, probability)
			probability = probability_bins[bin_index]
			#print(temp, bin_index, probability)
			rand_num=np.random.rand()
			binary_status = rand_num < probability
			is_binary.append(binary_status)

		binary_count = sum(is_binary)

		print(f"Number of identified binary stars: {binary_count}")

		return is_binary

	def draw_companion_m_ratio(self):
		"""
		Mass ratio (M2/M1) calculation from toy probability function based on Figure 16 of Raghavan et. al 2010

		Returns
		-------
		random_mass_ratio
			float value for the mass ratio (M2/M1) of the binary system
		"""
		# Normalization factor (maximum value of N on Figure 16 of Raghavan 2010)
		normalization_factor = 13

		# Randomly pull a number between zero and 1
		#random_number = np.random.rand()

		# Separate systems into 3 different sections based on Raghavan 2010, Figure 16
		#sectors = [1, 2, 3]
		#probability_bins = np.array([0.0909, 0.7909, 0.1182])*normalization_factor
		probability_bins = [0.0909, 0.7909, 0.1182]
		#sector_from_probability = {0.0909:1, 0.7909:2, 0.1182:3}

		# Plug in random number and pull a mass ratio from the toy probability function inspired by Figure 16 of D. Raghavan 2010
		def pull_from_distribution(random_number):
#			random_probability = random_number * normalization_factor
#			if random_probability < 5.5:
				#mass_ratio = random_probability / 22
				#mass_ratio = random_probability * 22
				#mass_ratio = random_probability * 27.5
#				mass_ratio = random_probability / 27.5
#			elif random_probability < 7.0:
				#mass_ratio = np.random.uniform(0.25, 1)
				#mass_ratio = np.random.uniform(0.25, 0.95)
#				mass_ratio = np.random.uniform(0.2, 0.95)
#			else:
				#mass_ratio = 1.0
				#mass_ratio = 0.95
#				mass_ratio = np.random.uniform(0.95, 1)
			####if random_number < (5.5 / normalization_factor):
			####	mass_ratio = random_number / 27.5
			####elif random_number < (7.0 / normalization_factor):
			####	mass_ratio = np.random.uniform(0.2, 0.95)
			####else:
			####	mass_ratio = np.random.uniform(0.95, 1)

			# Probability function for the first, linear section
			if random_number < probability_bins[0]:
				#print("Sector 1")
				#mass_ratio = random_number / 27.5
				#mass_ratio = random_number / 22
				P = np.random.rand()  # Uniform random number for inversion method
				mass_ratio = np.sqrt((2 * P) / 27.5)
			elif random_number < probability_bins[1]:
				#print("Sector 2")
				mass_ratio = np.random.uniform(0.2, 0.95)
				#mass_ratio = np.random.uniform(0.25, 0.95)
			else: 
				#print("Sector 3")
				mass_ratio = np.random.uniform(0.95, 1)
			#print(random_number, random_probability, mass_ratio)
			#print(random_number, mass_ratio)
			return mass_ratio

		# Normalize mass ratio distribution
		#distribution_max = mass_ratio_distribution(1)
		#def normalize_distribution(x):
		#	return mass_ratio_distribution(x) / distribution_max

		# Invert mass ratio distribution to solve for mass in terms of N
		#def inverse_mass_distribution(y):
		#	if y <= 0.25 * 24 * 0.25**2:
		#		return (2 * y / 24)**0.5
		#	elif y <= 0.25 * 24 * 0.25**2 + (0.65 * 5.5) * (0.9 - 0.25):
		#		return 0.25 + ((y - 0.25 * 24 * 0.25**2) / (0.65 * 5.5))**0.5
		#	else:
		#		return 0.9 + (y - 0.25 * 24 * 0.25**2 - (0.65 * 5.5) * (0.9 - 0.25)) / 13

		# Generate random numbers from the distribution
		#def generate_random():
		#	return inverse_mass_distribution(random.uniform(0, 1))
		#random_y = random.uniform(0, max(24 * 0.25**2 * 0.25, 0.25 * 24 * 0.25**2 + (0.65 * 5.5) * (0.9 - 0.25), 1))

		# Get the corresponding x value using the inverse CDF
		#random_x = inverse_cdf(random_y)

		# Randomly pull a number between zero and 1
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
			float value of the logarithn of period
		"""
		# Draw a log(Period) from the Raghavan Figure 13 Gaussian
		mu = 5.03	# Value from Raghavan Figure 13
		sigma = 2.28	# Value from Raghavan Figure 13
		logP = np.random.normal(loc = mu, scale = sigma)

		# Keep drawing if we get a period less than 1 day (logP < 0)
		while logP <= 0 or logP >= 9:
			logP = np.random.normal(loc=mu, scale=sigma)

		#return 10**logP
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
		
		#print(dataframe.columns)

		# Add a new column with random 6-digit identification numbers
		num_samples = len(dataframe)
		#print(len(dataframe))
		#print(dataframe)
		###identification_numbers = np.random.randint(100000, 999999, size=num_samples)
		###dataframe['ID'] = identification_numbers

		# Make a new column for identifying binaries [0=single, 1=primary, 2=secondary]
		dataframe['Is_Binary'] = 0
		#dataframe['Is_Binary'] = dataframe['Is_Binary'].astype(int)

		# Make a new column that will point secondaries to the primary id (0 by default since initial dataframe stars are single)
		dataframe['primary_ID'] = 0

		print("Range of temperatures: "+str(min(10**dataframe['logTeff']))+", "+str(max(10**dataframe['logTeff'])))

		# Identify binary stars and add new rows
		binary_flags = Binary.temp_is_binary(10**dataframe['logTeff'])	# 10** because Teff is log

		for i, is_binary in enumerate(binary_flags):
			#print(i, is_binary)
			if is_binary:
				# Set the Is_Binary column to True for the identified binary star
				dataframe.at[i, 'Is_Binary'] = 1
				#print(dataframe.at[i, 'Is_Binary'])

		# Add a new column with a nine-digit identification number
		####identification_numbers = np.random.randint(100000000, 999999999, size=len(dataframe))
		#dataframe['ID'] = np.array([f"{i:09d}" for i in range(1, len(dataframe) + 1)])
		#print(dataframe['ID'])
		# Multiply by 10 so that there is a zero at the end of all IDs
		####dataframe['ID'] = identification_numbers * 10#+ added_companions_df['Is_Binary']


		# Insert rows below binary stars
		new_rows = []
		####prim_sec_data = []
		#print("##################################")
		for index, row in dataframe.iterrows():
			#print(index)
			#new_rows.append(row)
			new_rows.append(row.to_frame().transpose())  # Convert the row to a DataFrame and append
			#print(new_rows)
			if row['Is_Binary'] == 1:
				#print("Is_Binary=1")
				# Insert a row with zeros below binary stars
				new_row = pd.DataFrame([[0] * dataframe.shape[1]], columns=dataframe.columns)
				new_row["Is_Binary"] = 2
				#print(new_row)
				new_rows.append(new_row)

				# Set identification numbers (add 1 for primaries and 2 to secondaries)
				#primary_id = row['ID'] + 1
				#secondary_id = row['ID'] + 2

				#row['ID'] = primary_id
				#new_row['ID'] = secondary_id

				# Make a second dataframe with primary and secondary IDs
				#prim_sec_data.append({'primary_ID': primary_id, 'secondary_ID': secondary_id})
			added_companions_df = pd.concat(new_rows, ignore_index=True)
			#print(added_companions_df)
			#print(index)
		#print("complete")
		####prim_sec_df = pd.DataFrame(prim_sec_data)
		#prim_sec_df.to_csv('/home/marznewman/spdev-3/synthpop-dev/synthpop/outputfiles/binary_populations/prim_sec.csv', index=False)
		#print(prim_sec_df)
		####prim_sec_df.to_csv(self.model.parms.output_location+'/prim_sec.csv', index=False)
		#print("complete")


		# Add a new column with a six-digit identification number
#		identification_numbers = np.random.randint(100000, 999999, size=len(added_companions_df))
		# Multiply by 10 so that there is a zero at the end of all IDs
#		added_companions_df['System_ID'] = identification_numbers * 10#+ added_companions_df['Is_Binary']
		added_companions_df['ID'] = np.array([f"{i:09d}" for i in range(1, len(added_companions_df) + 1)])
		#print(added_companions_df['ID'])
		#print("complete")

		# Add the primary IDs to "primary_ID" columns of the secondary
		for i in range(1, len(added_companions_df)):
			#print(i)
			if added_companions_df.loc[i, "Is_Binary"] == 2:
				added_companions_df.loc[i, "primary_ID"] = added_companions_df.loc[i - 1, "ID"]

		# Hard-coded files/directories for making an instance of Parameters
		####specific_config = "/home/marznewman/spwork/config_files/binary_pop_test.synthpop_conf"
		#default_config = "/home/marznewman/spwork/config_files/_default.synthpop_conf"
		default_config = DEFAULT_CONFIG_FILE
		#model_dir = "/home/marznewman/spwork/models/besancon_Robin2003"
		model_dir = DEFAULT_MODEL_DIR+"/"+self.model.parms.model_name

		# Get a list of population json files
		population_files = sorted(glob.glob(os.path.join(model_dir, "*pop.json")) + glob.glob(os.path.join(model_dir, "*popjson")))

		# Make a dictionary with the population id and corresponding population json file
		population_params = {
		pop_id: sp_utils.PopParams.parse_jsonfile(population_file)
		for pop_id, population_file in enumerate(population_files)
		}

		# Parameters for making an instance of Populations, taken from the config files and model directory
		####glbl_params = sp_utils.Parameters(specific_config, default_config, model_dir)

		# Get glbl_params by iterating over population numbers 0 through 9 (may cause bugs in the future)
		#print(self.model.populations[0].glbl_params)
		for i in range(0,10):
			#print(i)
			glbl_params = self.model.populations[i].glbl_params	

		# Instance of Populations
		populations = [
			Population(pop_params, pop_id, glbl_params)
			for pop_id, pop_params in population_params.items()
			]

		primaries = added_companions_df[added_companions_df['Is_Binary']==1].reset_index(drop=False)
		position = primaries[['x', 'y', 'z', 'Dist', 'l', 'b']].to_numpy()

		# Trying to figure out how to pull max_mass from the configuration file
		#print(default_config)
		#print(default_config.keys() if isinstance(default_config, dict) else dir(default_config))
		#max_mass = default_config.get('max_mass', None)
		max_mass = self.model.parms.mass_lims['max_mass']

		######## Try to use the assign_subclass method ########
		for pop_id, population in enumerate(populations):
			(population_density, imf, age, metallicity, kinematics, evolution, extinction) = populations[pop_id].assign_subclasses()
			#populations[pop_id].max_mass = populations[pop_id].glbl_params.mass_lims['max_mass']
			populations[pop_id].generator = StarGenerator(imf, age, metallicity, evolution, glbl_params, position, max_mass, logger)
			#print("AAAAAAAA")
			#print(dir(populations))
			#print(populations[pop_id].extinction)
			#print(extinction.get_extinctions())
			#print(position)
			#print(position[:, 3:6])
			#####extinction_in_map, extinction_dict = populations[pop_id].extinction.get_extinctions(position[:, 4], position[:, 5], position[:, 3])	# turns out I was actually already getting extinction_in_map from a line below
			
		
		#print(extinction_in_map)
			
		######## Try to make an instance of StarGenerator ########
		#generator = StarGenerator(
		#	imf, age, metallicity, evolution,
		#	glbl_params, logger
		#	)

		####### Try to use the get_evolved_props method ########
		bands = list(evolution.bands)
		# requested properties
		props = set(const.REQ_ISO_PROPS + glbl_params.opt_iso_props + bands)

		# Define a dataframe with just primaries
#		primaries = added_companions_df[added_companions_df['Is_Binary']==1].reset_index(drop=False)
		#print(primaries)

		# Define array for primary metallicities and ages to use for secondaries
		metallicities = primaries['Fe/H_initial']
		#print(metallicities)
		ages = primaries['age']
		#print(ages)
		prim_masses = primaries['iMass']
		#print(prim_masses)
		#print(min(prim_masses), max(prim_masses))

		# Create an instance of the Binary subclass for each companion
		companions = [Binary() for i in range(len(primaries))]
		#print(companions)

		'''
		# Make an array of initial masses, pulled from a mass ratio distribution
		mass_ratios = []
		for i in range(len(primaries)):
			mass_ratios.append(companions[i].draw_companion_m_ratio())
		masses = np.array(mass_ratios * prim_masses)	# Array of initial masses
		'''

		mass_ratios = [companion.draw_companion_m_ratio() for companion in companions]
		#print("\n\nMass Ratios")
		#print(mass_ratios)
		masses = np.array(mass_ratios * prim_masses)	# Array of companion initial masses
		#print(masses)
		
		# Combine primary/single stars with the companions
		

		# Debug check
		#for i in range(0, len(prim_masses)):
		#	print(prim_masses[i], masses[i])

		#print(mass_ratios)
		
		#print("AAAAAAAAAAAAAAA")
		#print(props)

		# Run get_evolved_props
		#ref_mag, s_props, final_phase_flag, inside_grid, not_evolved = generator.get_evolved_props(masses, metallicities, ages, props)
		for pop_id, population in enumerate(populations):
			ref_mag, s_props, final_phase_flag, inside_grid, not_evolved = populations[pop_id].generator.get_evolved_props(masses, metallicities, ages, props)

		# Run extract_properties
		for pop_id, population in enumerate(populations):
			m_evolved, props, user_props = populations[pop_id].extract_properties(masses, s_props, const.REQ_ISO_PROPS, glbl_params.opt_iso_props, inside_grid, not_evolved)

		###### Define variables needed for running extract magnitudes ######
		total_stars = len(companions)
		#total_stars = len(added_companions_df)
		#print(total_stars)

		# get maximum distance and step_size
		for pop_id, population in enumerate(populations):
			max_distance = getattr(populations[pop_id].pop_params, 'max_distance', glbl_params.max_distance)
			#step_size = getattr(populations[pop_id].pop_params, 'distance_step_size', glbl_params.distance_step_size)

		#print(max_distance)
		#print(step_size)

		# Radii we will step over
		#radii = np.arange(0, max_distance + step_size, step_size)
		radii = np.linspace(0, max_distance, total_stars+1)
		#print(radii.shape)

		# reduce number of stars by the scale factor and fract_above_min_mass
		#total_stars = np.random.poisson(n_star_expected * (1 - frac_lowmass[1]) / self.scale_factor)

		# Hard-coded variable because I have no idea what it is and it requires me to go down a rabbit hole
		#missing_stars = [100000, 100001]

		# Inner radius of the slice
		#print(radii[:-1])
		#r_inner = np.repeat(radii[:-1], total_stars)
		r_inner = radii[:-1]
		#print(r_inner.shape)


#		for pop_id, population in enumerate(populations):
		#position = np.vstack([np.column_stack(populations[pop_id].position.draw_random_point_in_slice(r_inner, r_outer, n_stars)) for r_inner, r_outer, n_stars in zip(radii, radii[1:], 1)])
#			positions_list = []
#			for r_inner, r_outer in zip(radii, radii[1:]):
#				x, y, z, d_kpc, star_l_deg, star_b_deg = populations[pop_id].position.draw_random_point_in_slice(r_inner, r_outer, total_stars)
#				positions_list.append(np.column_stack((x, y, z, d_kpc, star_l_deg, star_b_deg)))

		# Define the position array from existing primaries
#		position = primaries[['x', 'y', 'z', 'Dist', 'l', 'b']].to_numpy()
		#print(position[:, 3:6])


		# Extract magnitudes and convert to observed magnitudes if needed
#		for pop_id, population in enumerate(populations):
#			#print(pop_id)
#			#print(population)
#			#print(populations[pop_id])
#			mags, extinction_in_map = populations[pop_id].extract_magnitudes(r_inner, position[:, 3:6], ref_mag, s_props)
#			print(len(mags))
#			print(mags)
#			print('\n\n\n\n\n')
		#mags, extinction_in_map = populations[pop_id].extract_magnitudes(r_inner, [3, 4, 5], ref_mag, s_props)
		
		#print(extinction_in_map)
		#print("population_params")
		#print(population_params)
		
		#print("AAAAA")
		#print(len(r_inner))
		#print(len(s_props))
		#print(position[:, 3:6])
		#print(s_props)
		
		#print(len(populations))
		#print(len(radii))
		#print(radii)
		
		#for pop_id, population in enumerate(populations):
		mags, extinction_in_map = populations[0].extract_magnitudes(r_inner, position[:, 3:6], ref_mag, s_props)	# Why do I only do this for pop_id = 0? I have no idea
		
		#print(mags)
		#print(extinction_in_map, len(extinction_in_map))
		#print(mags)

		############ Add new properties to companion stars in dataframe ############
		# Define a list of initial parameters
		initial_parameters = np.column_stack([masses, ages, metallicities])

		# collect all the column names
		# required_properties + optional_properties + magnitudes
		#headers = const.COL_NAMES + populations[pop_id].glbl_params.col_names + populations[pop_id].bands
		headers = const.COL_NAMES + populations[0].glbl_params.col_names + populations[0].bands	# Once again not iterating over populations, just using zero which worries me
		#print(headers)
		
		# replace "ExtinctionInMap" with the output of the extinction map
		if 'ExtinctionInMap' in headers:
			extinction_index = headers.index("ExtinctionInMap")
			#print(extinction_index)
			headers[extinction_index] = populations[0].extinction.A_or_E_type	# The fact that I only use the 0 population worries me, but it looked like it was the same for all populations
			
		#print(headers)

		# Fill a couple of arrays filled with nans for kinematic parameters
		#proper_motions = np.zeros(len(companions))
		#velocities = np.zeros(len(companions))
		#vr_lsr = np.zeros(len(companions))
		proper_motions = np.full((len(companions), 3), np.nan)
		velocities = np.full((len(companions), 3), np.nan)
		vr_lsr = np.repeat(np.nan, len(companions))

		# Convert Table to pd.DataFrame
		#companion_properties_df = populations[pop_id].convert_to_dataframe(populations[pop_id].popid, initial_parameters, m_evolved, final_phase_flag, position[:, 3:6], proper_motions, position[:, 0:3], velocities, vr_lsr, extinction_in_map, props, user_props, mags, headers)
		# I think I had the position arguments backwards before
		companion_properties_df = []
		for pop_id, population in enumerate(populations):
			df = populations[pop_id].convert_to_dataframe(populations[pop_id].popid, initial_parameters, m_evolved, final_phase_flag, position[:, 0:3], proper_motions, position[:, 3:6], velocities, vr_lsr, extinction_in_map, props, user_props, mags, headers)
			companion_properties_df.append(df)
		
		# Convert to pandas dataframe
		companion_properties_df = pd.concat(companion_properties_df)

		############# Add new df information to final catalog dataframe #############
		companion_indices = added_companions_df[added_companions_df['Is_Binary'] == 2].index
		#print(companion_indices)

		for i, idx in enumerate(companion_indices):
			added_companions_df.loc[idx, companion_properties_df.columns] = companion_properties_df.iloc[i]

		############# Make a new dataframe of combined properties #############
		# Separate primary and secondary stars
		###primary_stars = added_companions_df[added_companions_df['Is_Binary'] == 1].reset_index(drop=True)
		###secondary_stars = added_companions_df[added_companions_df['Is_Binary'] == 2].reset_index(drop=True)
		primary_stars = added_companions_df[added_companions_df['Is_Binary'] == 1]	# Dont drop the index because we'll add it to the main dataframe
		secondary_stars = added_companions_df[added_companions_df['Is_Binary'] == 2]

		#print(primary_stars)
		#print("----------------------------")
		#print(secondary_stars)
		#print("\n\n\n")

		# Initialize new combined property columns with NaN values
		added_companions_df['combined_logL'] = np.nan
		added_companions_df['total_mass'] = np.nan
		added_companions_df['q'] = np.nan
		added_companions_df['combined_logP'] = np.nan

		# Draw periods for binary stars
		periods = [companion.draw_period() for companion in companions]

		# Calculate combined luminosities using a loop
		#combined_magnitudes = []	# I don't think I need this
		###combined_luminosities = []
		###combined_masses = []	
		#print('\n\n\n\nLuminosities')
		#for i in range(len(primary_stars)):
		for prim_index, row in primary_stars.iterrows():
			#print(prim_index)
			# Combine luminosities
			logL1 = primary_stars.loc[prim_index, 'logL']
			logL2 = secondary_stars.loc[prim_index+1, 'logL']
			combined_lum = np.log10(10**logL1 + 10**logL2)
			###combined_luminosities.append(combined_lum)
			#print(logL1, logL2)

			# Combine mass
			mass1 = primary_stars.loc[prim_index, 'Mass']
			mass2 = secondary_stars.loc[prim_index+1, 'Mass']
			combined_mass = mass1 + mass2
			###combined_masses.append(combined_mass)

			# Add combined properties to the primary star and secondary star rows
			added_companions_df.loc[prim_index, 'combined_logL'] = combined_lum
			added_companions_df.loc[prim_index, 'total_mass'] = combined_mass
			#added_companions_df.loc[prim_index, 'logP'] = 
			added_companions_df.loc[prim_index+1, 'combined_logL'] = combined_lum
			added_companions_df.loc[prim_index+1, 'total_mass'] = combined_mass


		#print(combined_luminosities)
		'''
		# Make a list of periods
		periods = []
		for i in range(len(primary_stars)):
			periods.append(companions[i].draw_period())
		'''

		for i, companion_index in enumerate(companion_indices):
			#print(i, companion_index)
			added_companions_df.loc[companion_index-1, 'q'] = mass_ratios[i]
			added_companions_df.loc[companion_index-1, 'combined_logP'] = periods[i]
			added_companions_df.loc[companion_index, 'q'] = mass_ratios[i]
			added_companions_df.loc[companion_index, 'combined_logP'] = periods[i]


		'''	
		# Create the combined_props_df DataFrame
		combined_props = {
			#'system_ID': primary_stars['ID']-1,
			'primary_ID': primary_stars['ID'],
			#'mag': combined_magnitudes,
			'combined_logL': combined_luminosities,
			'total_mass' : combined_masses,
			'q': mass_ratios,
			'logP': periods
		}

		combined_props_df = pd.DataFrame(combined_props)
		'''

		# Output
		#combined_props_df.to_csv('/home/marznewman/spdev-3/synthpop-dev/synthpop/outputfiles/binary_populations/combined_props.csv', index=False)
		####combined_props_df.to_csv(self.model.parms.output_location+'/combined_props.csv', index=False)	# This was when I created a separate output file for combined properties

		### Check if I even have the luminosities of companions
		#print(added_companions_df['logL'])
		#added_companions_df[['logL']].to_csv('logL_column.csv', index=False)
		
		# Check if extinctions are not nans/zeros
		#print(added_companions_df["ExtinctionInMap"])
		#print(added_companions_df["A_Ks"])
		#print(companion_indices)

		return added_companions_df

