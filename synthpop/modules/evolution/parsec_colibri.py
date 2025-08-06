"""
This module provides the loading and handling of PARSEC-Colibri isochrones.

Generated via CMD 3.8 (http://stev.oapd.inaf.it/cmd) on Tue Jul 22 19:14:39 UTC 2025.
Isochrones based on PARSEC release v1.2S +  COLIBRI S_37 + S_35 + PR16.
"""
__all__ = ["PARSEC",]
__author__ = "M.J. Huston"
__date__ = "2025-08-04"

import os
import warnings
import json
import tqdm
import sys
import time
import pdb

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from io import StringIO

import numpy as np
import pandas as pd
import requests

from ._evolution import EvolutionIsochrones, ISOCHRONES_DIR, EVOLUTION_DIR, FILTERS_DIR
# import a "standard" interpolator
from .charon_interpolator import CharonInterpolator
#from .lagrange_interpolator import LagrangeInterpolator

# global variable to store the isochrones
pc_isochrones = None

class ParsecColibri(EvolutionIsochrones, CharonInterpolator):
    """
    Parsec-Colibri Isochrone class

    Attributes
    ----------
    FOLDER : location to store the isochrones

    isochrones : dict

    min_mass :

    max_mass :

    met_to_file_iso :

    non_mag_cols :

    file_met :

    bands :

    iso_ages :

    mass_range :


    Methods
    -------
    __init__(self, chosen_bands, iso_props, use_global=False) : None

    get_mass_ranges(isochrones) :  dict
        estimate the mass range for each metallicity age

    convert_to_met(file_isochrones) : ndarray

    load_isochrones(self, magsys) : dict

    convert_to_cmd_to_h5(filename) : pandas.DataFrame
        convert the ascii file to an .h5 file
        Returns the loaded table

    download_isochrones(self, magsys_name) : None
        Download the isochrone system and unpack them
    """

    # folder where isochrone files can be found
    FOLDER = f"{ISOCHRONES_DIR}/PARSEC_COLIBRI"

    pc_non_mag_cols = ['Zini', '[Fe/H]_init', 'log10_isochrone_age_yr', 'initial_mass', 'int_IMF', 'star_mass', 'log_L', 'log_Teff', 
                    'log_g', 'phase', 'McoreTP', 
                    'C_O', 'period0', 'period1', 'period2', 'period3', 'period4', 'pmode', 'Mloss', 'tau1m', 
                    'X', 'Y', 'Xc', 'Xn', 'Xo', 'Cexcess', 'Z', 'mbol']
    pc_req_cols_to_mist = {'MH':'[Fe/H]_init', 'logAge':'log10_isochrone_age_yr', 'Mini':'initial_mass', 
                            'Mass':'star_mass', 'label':'phase', 'mbolmag':'mbol', 'logL':'log_L',
                            'logTe':'log_Teff', 'logg':'log_g'}
    pc_phase_to_mist = {0:-1, 1:0, 2:2, 3:2, 4:3, 5:3, 6:3, 7:4, 8:5, 9:6}

    # lowest and highest mass where these isochrones should be used
    # (can be outside the covered range, in such cases the closets grid_points are used)
    # used area can be overwritten in the config.json file
    max_mass = 350
    min_mass = 0.1
    isochrones_name = 'PARSEC_COLIBRI'

    def __init__(self, columns, use_global=True, **kwargs):
        """
        pull out all the information from the isochrone file
        and puts it into tracks under the index of the lowest
        non-magnitude properties in MIST isochrones

        Parameters
        ----------
        columns : list
            list of columns
        use_global : Bool
            store or use isochrones as global variable
        """

        with open(f"{FILTERS_DIR}/effective_wavelengths_parsec.json") as f:
            self.all_filter_eff_wavs = json.load(f)
        with open(f"{FILTERS_DIR}/parsec_columns.json") as f:
            self.parsec_columns = json.load(f)

        self.magsys, self.non_mag_cols, self.bands = self.get_mag_systems(columns)

        # Check for isochrone directory and create if needed
        os.makedirs(self.FOLDER, exist_ok=True)

        self.file_met = np.array([-2.2, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5])
        self.met_to_file_iso = [self.convert_met_to_string(met) for met in self.file_met]

        # Check for isochrones in requested magnitude systems 
        # download from MIST if needed
        for msys in self.magsys:
            if not os.path.isdir(f"{self.FOLDER}/{msys}"):
                self.download_isochrones(msys)
                OSError("Isochrone files not found.")

        # load the isochrones
        if use_global:
            global pc_isochrones
            if pc_isochrones is not None:
                if np.array_equal(np.sort(np.unique(self.bands+self.non_mag_cols)).astype(str), np.sort(np.unique(pc_isochrones.keys())).astype(str)):
                    self.isochrones = pc_isochrones
                else:
                    pc_isochrones = self.isochrones = self.load_isochrones(self.magsys)
            else:
                pc_isochrones = self.isochrones = self.load_isochrones(self.magsys)
        else:
            self.isochrones = self.load_isochrones(self.magsys)
        self.isochrones_grouped = self.isochrones.groupby(["[Fe/H]_init", "log10_isochrone_age_yr"])

        # Get ages
        self.iso_ages = 10 ** self.isochrones['log10_isochrone_age_yr'].unique()

        # Get mass range for each metallicity and age
        self.mass_range = self.get_mass_ranges(self.isochrones_grouped)

        # call super after loading the isochrones so the interpolator has access to the data
        super().__init__(**kwargs)

    def get_mag_systems(self, columns):
        """
        Parameters
        ----------
        columns : list
            list specifying columns and filter systems from the MIST Isochrones
        Returns
        -------
        collected_magsys : dict
            Contains the columns for each needed filter system
        non_mag_cols : list
            a list of all non-magnitude columns
        all_bands : list
            a list of all magnitude columns
        """
        # global mist_columns
        # if len(mist_columns) == 0:
        #     with open(f"{EVOLUTION_DIR}/mist_columns.json") as f:
        #         mist_columns.update(json.load(f))

        # placeholder for output
        collected_magsys = {}
        all_bands = []
        non_mag_cols = []

        for col in columns:
            if isinstance(col, list):
                raise TypeError("'chosen_bands' must be a list of photometric systems or a " + \
                                "dictionary of systems and select filters for PARSEC_COLIBRI isochrones")
            elif isinstance(col, tuple):
                # magsys filters pairs are provided
                magsys, bands = col
                # check if parsec has this filter system
                if magsys not in self.parsec_columns['phot_short_names']:
                    raise ValueError(f"Can not assign {magsys}")

                if isinstance(bands, str):
                    bands = bands,
                magsys_all_bands = self.get_iso_phot_filt_list(magsys)
                if bands[0] == "all":
                    # get all filter for this filter system
                    bands = self.get_iso_phot_filt_list(magsys)
                else:
                    not_found = []
                    # check if all bands belongs to magsys
                    for i in bands:
                        if i not in magsys_all_bands:
                            not_found.append(i)
                    if len(not_found) != 0:
                        raise ValueError(f"Can not assign {not_found} to {magsys}")

            elif col in self.pc_non_mag_cols:
                magsys = 'any'
                bands = [col, ]

            elif col in self.parsec_columns['phot_short_names']:
                # filter systems are passed
                magsys = col
                bands = self.get_iso_phot_filt_list(magsys)
            else:
                raise ValueError(f"Can not assign {col} to a magnitude system")

            if magsys=='any':
                non_mag_cols.extend(bands)
            else:
                bands_long= (magsys+'_'+b for b in bands)
                all_bands.extend(bands_long)

            # create entry in collected_magsys
            # use set to avoid duplicates
            if magsys not in collected_magsys:
                collected_magsys[magsys] = set()

            collected_magsys[magsys].update(bands)

        # attach non-mag columns to first phot system data
        if "any" in collected_magsys:
            val = collected_magsys.pop("any")
            first_key = list(collected_magsys.keys())[0]
            collected_magsys[first_key].update(val)

        return collected_magsys, non_mag_cols, all_bands

    @staticmethod
    def get_mass_ranges(isochrones_grouped):
        """
        Estimate the mass range for each metallicity and age
        Parameters
        ----------
        isochrones_grouped : pandas.DataFrameGroupeBy

        Returns
        -------
        mass_ranges : pandas.DataFrame
            mass range for each age and metallicity
        """

        max_values = isochrones_grouped['initial_mass'].max()
        min_values = isochrones_grouped['initial_mass'].max()
        min_values.name='min_mass'
        max_values.name='max_mass'
        mass_range = pd.concat([min_values,max_values], axis=1)
        return mass_range

    def convert_met_to_string(self, met):
        return ('m' if met<0 else 'p')+f'{np.abs(met):1.2f}'

    def get_filenames(self, magsys, met_string):
        filename_pattern = f"{self.FOLDER}/{magsys}/" \
                           f"PARSECv1.2_COLIBRI_mh_{met_string}_{magsys}.dat"

    @staticmethod
    def get_columns(filename):
        with open(filename) as f:
            for line in f:
                if line.startswith('#'):
                    prev = line
                else:
                    break
        return prev[2:-1].split()

    def load_isochrones(self, magsys):
        """
        load the isochrones for each of the given magnitude systems and all metallicities
        if needed it will download the isochrones from the server

        Parameters
        ----------
        magsys : list
            list of magnitude systems.

        Returns
        -------
        isochrones : dict
            dictionary of separate PandasDataFrames for each metallicity.
        """

        isochrones = {}
        # Load columns
        for i, (magsys_name, chosen_columns) in enumerate(magsys.items()):
            # Get file header from one file for column names
            use_columns = [(magsys_name+'_'+c if (c not in self.non_mag_cols) else c) for c in chosen_columns]

            # Load in columns for each metallicity
            j=0
            for met_string, file_met in zip(self.met_to_file_iso, self.file_met):
                filename = self.get_filenames(magsys_name, met_string)
                if not os.path.isfile(f'{filename}.dat'):
                    self.download_isochrones(magsys_name, select_met_idx=j)

                if os.path.isfile(f'{filename}.h5'):
                    # file exist in a Hierarchical Data Format (.h5)
                    df = pd.read_hdf(f'{filename}.h5', 'data')
                else:
                    # convert ascii to Hierarchical Data Format
                    df = self.convert_dat_to_h5(filename)
                # Adjust column names as needed
                adjust_columns = self.pc_req_cols_to_mist.copy()
                adjust_columns.update({b+'mag':magsys_name+'_'+b for b in chosen_columns})
                df.rename(columns=adjust_columns, inplace=True)
                if file_met==-2.2:
                    df['[Fe/H]_init'] = -2.2
                df.sort_values(by=['[Fe/H]_init','log10_isochrone_age_yr', 'initial_mass'], inplace=True)
                for n_phase in self.pc_phase_to_mist:
                    df.loc[df['phase']==n_phase, 'phase'] = self.pc_phase_to_mist[n_phase]

                if i == 0:
                    isochrones[file_met] = df[use_columns].copy()
                else:
                    isochrones[file_met][use_columns] = df[use_columns]
                j+=1

        return pd.concat(isochrones.values())

    def convert_dat_to_h5(self, filename):
        """
        Convert the isochrones from an ascii table into a Hierarchical Data Format file (.h5)
        Parameters
        ----------
        filename : str
            filename of the ascii file

        Returns
        -------
        table : pandas.DataFrame
            isochrone table
        """
        print(f"convert {filename} to hdf5", end="\r")
        df = pd.read_csv(filename+'.dat', comment='#', sep='\s+', names=self.get_columns(filename+'.dat'))
        # save table as h5 file
        df.to_hdf(filename+'.h5', key='data')
        return df

    def get_filenames(self, magsys, met_string):
        return f"{self.FOLDER}/{magsys}/PARSECv1.2_COLIBRI_mh_{met_string}_{magsys}"

    def get_iso_phot_file_name(self, phot_sys_short):
        idx = np.where(phot_sys_short==np.array(self.parsec_columns['phot_short_names']))[0][0]
        return self.parsec_columns['phot_file_strs'][idx]

    def get_iso_phot_filt_list(self, phot_sys_short):
        idx = np.where(phot_sys_short==np.array(self.parsec_columns['phot_short_names']))[0][0]
        return self.parsec_columns['phot_filt_lists'][idx]

    def download_isochrones(self, magsys, select_met_idx=None):
        """
        Generate and download the PARSEC-COLIBRI isochrones

        Parameters
        ----------
        magsys_name : string
            magnitude system
        Returns
        -------

        """

        # Open CMD 3.8 web form
        driver = webdriver.Chrome() # Or Firefox, Edge, etc.
        driver.get("https://stev.oapd.inaf.it/cgi-bin/cmd")
        driver.implicitly_wait(0.2)
        os.makedirs(f"{self.FOLDER}/{magsys}",exist_ok=True)

        if select_met_idx is None:
            print(f'Beginning isochrone generation for {magsys}. This may take a couple minutes the first time.')
            loop_met_file, loop_met = self.met_to_file_iso, self.file_met
        else:
            loop_met_file, loop_met = [self.met_to_file_iso[select_met_idx]], [self.file_met[select_met_idx]]
        #Loop over metallicity range
        for met_str, met in zip(loop_met_file, loop_met):
            print(f'Beginning isochrone generation for {magsys}, [M/H]={met}.')

            # Select evolutionary tracks
            parsec_button = driver.find_element(By.CSS_SELECTOR, "input[type='radio'][value='parsec_CAF09_v1.2S']")
            parsec_button.click()
            colibri_button = driver.find_element(By.CSS_SELECTOR, "input[type='radio'][value='parsec_CAF09_v1.2S_S_LMC_08_web']")
            colibri_button.click()

            #Set additional evolution parameters
            n_inTPC_box = driver.find_element(By.NAME, "n_inTPC")
            n_inTPC_box.clear()
            n_inTPC_box.send_keys("10")
            eta_reimers_box = driver.find_element(By.NAME, "eta_reimers")
            eta_reimers_box.clear()
            eta_reimers_box.send_keys("0.2")

            # Select photometric system and bolometric corrections
            phot_dropdown = driver.find_element(by=By.NAME, value="photsys_file")
            phot_select = Select(phot_dropdown)
            phot_select.select_by_value(value=self.get_iso_phot_file_name(magsys))

            ybc_button = driver.find_element(By.CSS_SELECTOR, "input[type='radio'][value='YBCnewVega']")
            ybc_button.click()

            # Circumstellar dust settings
            dust_m_button = driver.find_element(By.CSS_SELECTOR, "input[type='radio'][value='dpmod60alox40']")
            dust_m_button.click()
            dust_c_button = driver.find_element(By.CSS_SELECTOR, "input[type='radio'][value='AMCSIC15']")
            dust_c_button.click()

            # Apply no extinction - we'll do this per star in SynthPop
            Av_box = driver.find_element(By.NAME, "extinction_av")
            Av_box.clear()
            Av_box.send_keys("0.0")

            # Long period variability
            lpv_button = driver.find_element(By.CSS_SELECTOR, "input[type='radio'][name='kind_LPV'][value='4']")
            lpv_button.click()

            # Initial mass function
            imf_dropdown = driver.find_element(by=By.NAME, value="imf_file")
            imf_select = Select(imf_dropdown)
            imf_select.select_by_value(value="tab_imf/imf_kroupa_orig.dat")

            # Set isochrone age sampling
            logage_button = driver.find_element(By.CSS_SELECTOR, "input[type='radio'][name='isoc_isagelog'][value='1']")
            logage_button.click()
            logage_min_box = driver.find_element(By.NAME, "isoc_lagelow")
            logage_min_box.clear()
            logage_min_box.send_keys("5.0")
            logage_max_box = driver.find_element(By.NAME, "isoc_lageupp")
            logage_max_box.clear()
            logage_max_box.send_keys("10.13")
            logage_d_box = driver.find_element(By.NAME, "isoc_dlage")
            logage_d_box.clear()
            logage_d_box.send_keys("0.05")

            # Set single metallicity value (can loop over entire thing with different metallicities to get full isochrone set)
            mh_button = driver.find_element(By.CSS_SELECTOR, "input[type='radio'][name='isoc_ismetlog'][value='1']")
            mh_button.click()
            mh_min_box = driver.find_element(By.NAME, "isoc_metlow")
            mh_min_box.clear()
            mh_min_box.send_keys(f"{met}")
            mh_d_box = driver.find_element(By.NAME, "isoc_dmet")
            mh_d_box.clear()
            mh_d_box.send_keys("0.0")

            # Select isochrone output and submit form
            out_button = driver.find_element(By.CSS_SELECTOR, "input[type='radio'][name='output_kind'][value='0']")
            out_button.click()
            out_button = driver.find_element(By.CSS_SELECTOR, "input[type='submit'][name='submit_form']")
            out_button.click()

            # Get the link for the output file, waiting if it needs more time to load
            time.sleep(5)
            while True:
                try:
                    elems = driver.find_elements(By.XPATH, "//a[@href]")
                    links = [elem.get_attribute('href') for elem in elems]
                    output_link = links[np.where([('output' in l) for l in links])[0][0]]
                    print(f'Beginning download.')
                    break
                except:
                    print('Download link not ready, wait 5 seconds and try again.')
                    time.sleep(5)

            # Save the data to a file
            with open(self.get_filenames(magsys, met_str)+'.dat', "wb") as f:
                r = requests.get(output_link, stream=True)        
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
                print('Isochrone saved.')

            # Return to web form
            reset_button = driver.find_element(By.CSS_SELECTOR, "input[type='submit'][name='reset_form']")
            reset_button.click()
            time.sleep(1)

        # Close the browser session
        driver.close()

        print('Isochrones downloaded: %s' % magsys)

        return
