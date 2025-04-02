"""
This module provides the loading and handling of MIST isochrones.
It automatically downloads the isochrone form the online MIST service.

MIST documentation: https://waps.cfa.harvard.edu/MIST/

MIST papers: Dotter (2016), Choi et al. (2016);
DOIs: 10.3847/0067-0049/222/1/8, 10.3847/0004-637X/823/2/102

The mist_columns.json file provides a guide to which isochrone files include which columns/filters.
"""
__all__ = ["MIST",]
__author__ = "M.J. Huston"
__date__ = "2023-03-01"

import glob
import importlib
import linecache
import os
import tarfile
import warnings
import json
import tqdm
import sys

import numpy as np
import pandas
import requests

from ._evolution import EvolutionIsochrones, ISOCHRONES_DIR, EVOLUTION_DIR
# import a "standard" interpolator
from .charon_interpolator import CharonInterpolator
#from .lagrange_interpolator import LagrangeInterpolator

# global variable to store the isochrones
mist_isochrones = None
mist_columns = {}

class MIST(EvolutionIsochrones, CharonInterpolator):
    """
    MIST Isochrone class

    Attributes
    ----------
    FOLDER : location to store the isochrones

    isochrones : dict

    min_mass :

    max_mass :

    met_to_file_iso :

    none_mag_cols :

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
    FOLDER = f"{ISOCHRONES_DIR}/mist"

    # lowest and highest mass where MIST isochrones should be used
    # (can be outside the covered range, in such cases the closets grid_points are used)
    # used area can be overwritten in the config.json file
    max_mass = 250
    min_mass = 0.1
    isochrones_name = 'MIST'

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

        self.magsys, self.none_mag_cols, self.bands = self.get_mag_systems(columns)

        # Check for isochrone directory and create if needed
        os.makedirs(self.FOLDER, exist_ok=True)

        # Check for isochrones in requested magnitude systems 
        # download from MIST if needed
        for msys in self.magsys:
            if not os.path.isdir(f"{self.FOLDER}/{msys}"):
                self.download_isochrones(msys)

        # Get metallicities
        feh_pos = sorted([x[x.find('_feh_') + 5:x.find('_feh_') + 10]
            for x in glob.glob(
                f'{self.FOLDER}/{msys}/MIST_*feh_p*.iso.cmd')])

        feh_neg = sorted([x[x.find('_feh_') + 5:x.find('_feh_') + 10]
            for x in glob.glob(
                self.FOLDER + '/' + msys + '/MIST_*feh_m*.iso.cmd')],
            reverse=True)

        self.met_to_file_iso = np.array(feh_neg + feh_pos)
        self.file_met = self.convert_to_met(self.met_to_file_iso)

        # load the isochrones
        if use_global:
            global mist_isochrones
            if mist_isochrones is not None:
                if np.array_equal(np.sort(np.unique(self.bands+self.none_mag_cols)).astype(str), np.sort(np.unique(mist_isochrones.keys())).astype(str)):
                    self.isochrones = mist_isochrones
                else:
                    mist_isochrones = self.isochrones = self.load_isochrones(self.magsys)
            else:
                mist_isochrones = self.isochrones = self.load_isochrones(self.magsys)
        else:
            self.isochrones = self.load_isochrones(self.magsys)
        self.isochrones_grouped = self.isochrones.groupby(["[Fe/H]_init", "log10_isochrone_age_yr"])

        # Get ages
        self.iso_ages = 10 ** self.isochrones['log10_isochrone_age_yr'].unique()

        # Get mass range for each metallicity and age
        self.mass_range = self.get_mass_ranges(self.isochrones_grouped)

        # call super after loading the isochrones so the interpolator has axcess to the data
        super().__init__(**kwargs)


    @staticmethod
    def get_mag_systems(columns):
        """
        Parameters
        ----------
        columns : list
            list specifying columns and filter systems from the MIST Isochrones
        Returns
        -------
        collected_magsys : dict
            Contains the columns for each needed filter system
        none_mag_cols : list
            a list of all non-magnitude columns
        all_bands : list
            a list of all magnitude columns
        """
        global mist_columns
        if len(mist_columns) == 0:
            with open(f"{EVOLUTION_DIR}/mist_columns.json") as f:
                mist_columns.update(json.load(f))

        # placeholder for output
        collected_magsys = {}
        all_bands = []
        none_mag_cols = []

        for col in columns:
            if isinstance(col, list):
                raise TypeError("'chosen_bands' must be a flat array or a dictionary, "
                                "see 'config_files/_default.synthpop_conf'")
            elif isinstance(col, tuple):
                # magsys filters pairs are provided
                magsys, bands = col
                if isinstance(bands, str):
                    bands = bands,
                # check if mist has this filter system
                if magsys not in mist_columns.values():
                    raise ValueError(f"Can not assign {magsys}")

                if bands[0] == "all":
                    # get all filter for this filter system
                    bands = tuple([filt for filt, ms in mist_columns.items() if ms == magsys])

                else:
                    not_found = []
                    # check if all bands belongs to magsys
                    for i in bands:
                        if i not in mist_columns:
                            not_found.append(i)
                        elif mist_columns[i] != magsys:
                            not_found.append(i)
                    if len(not_found) != 0:
                        raise ValueError(f"Can not assign {not_found} to {magsys}")

            elif col in mist_columns:
                # filters are passed as list  or columns

                # get filter system for filter
                magsys = mist_columns.get(col)
                bands = [col, ]

            elif col in mist_columns.values():
                # filter systems are passed
                magsys = col

                # get all filter for this filter system
                bands = tuple([filt for filt, ms in mist_columns.items() if ms == magsys])
            else:
                raise ValueError(f"Can not assign {col} to a magnitude system")

            if magsys in ['cmd', 'basic', 'full']:
                none_mag_cols.extend(bands)
            else:
                all_bands.extend(bands)

            # create entry in collected_magsys
            # use set to avoid duplicates
            if magsys not in collected_magsys:
                collected_magsys[magsys] = set()

            collected_magsys[magsys].update(bands)

        # clean up
        # attach basic to full if full needs to be loaded anyhow
        if "basic" in collected_magsys and "full" in collected_magsys:
            val = collected_magsys.pop("basic")
            collected_magsys["full"].update(val)

            # attach columns to first cmd data
        if "cmd" in collected_magsys:
            val = collected_magsys.pop("cmd")
            first_key = list(collected_magsys.keys())[:2]
            i = 0
            while first_key[i] in ['basic', 'full']:
                i += 1
            collected_magsys[first_key[i]].update(val)

        return collected_magsys, none_mag_cols, all_bands

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
        mass_range = pandas.concat([min_values,max_values], axis=1)
        return mass_range

    @staticmethod
    def convert_to_met(file_isochrones):
        return np.array(
            [float(x.replace('m', '-').replace('p', '')) for x in file_isochrones]
            )

    def get_filenames(self, magsys, met_string):
        filename_pattern = f"{self.FOLDER}/{magsys}/" \
                           f"MIST_v1.2_feh_{met_string}_afe_p0.0_vvcrit0.4_{magsys}.iso"

        if magsys in ["basic", "full"]:
            return filename_pattern
        else:
            return filename_pattern + ".cmd"

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
            use_columns = list(chosen_columns)

            # Load in columns for each metallicity
            for met_string, file_met in zip(self.met_to_file_iso, self.file_met):
                filename = self.get_filenames(magsys_name, met_string)

                # check if PyTables is installed
                if importlib.util.find_spec("tables") is None:
                    warnings.warn("PyTables should be installed "
                                  "to significantly speed up the loading process")
                    df = self.read_csv(filename)

                else:
                    if os.path.isfile(f'{filename}.h5'):
                        # file exist in a Hierarchical Data Format (.h5)
                        df = pandas.read_hdf(f'{filename}.h5', 'data')

                    else:
                        # convert ascii to Hierarchical Data Format
                        df = self.convert_to_cmd_to_h5(filename)
                # Correct Phase for WD
                df.loc[df['EEP'] > 808, 'phase'] = 5
                df.loc[df['EEP'] > 1409, 'phase'] = 6
                df.loc[(df['EEP'] > 707) & (df['phase'] == 0) , 'phase'] = 4
                df.loc[(df['EEP'] > 605) & (df['phase'] == 0) , 'phase'] = 3
                df.loc[(df['EEP'] > 453) & (df['phase'] == 0) , 'phase'] = 2

                if i == 0:
                    isochrones[file_met] = df[use_columns].copy()
                else:
                    isochrones[file_met][use_columns] = df[use_columns]

        return pandas.concat(isochrones.values())

    @staticmethod
    def get_columns(filename):
        with open(filename) as f:
            for line in f:
                if line.startswith('#'):
                    prev = line
                else:
                    break
        return prev[2:-1].split()

    @classmethod
    def read_csv(cls, filename):
        # get column names
        cols = cls.get_columns(filename)
        # load table from ascii file
        df = pandas.read_csv(filename, sep='\s+', comment='#',
            skip_blank_lines=True, low_memory=False, header=None, names=cols)
        return df

    @classmethod
    def convert_to_cmd_to_h5(cls, filename):
        """
        convert the MIST isochrones from an ascii table into a Hierarchical Data Format file (.h5)
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
        df = cls.read_csv(filename)
        # save table as h5 file
        df.to_hdf(f'{filename}.h5', key='data', mode='w')
        return df

    def download_isochrones(self, magsys_name):
        """
        Download the MIST isochrones from the Server

        Parameters
        ----------
        magsys_name : string
            magnitude system
        Returns
        -------

        """
        filename = f'{self.FOLDER}/MIST_v1.2_vvcrit0.4_{magsys_name}.txz'
        suffix = "_isos" if magsys_name in ["basic", "full"] else ""
        url = (f'http://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/'
               f'MIST_v1.2_vvcrit0.4_{magsys_name}{suffix}.txz')
        filesize = {"basic": 221e6, "full": 630e6}.get(magsys_name, 0)
        chunk_size = 1024
        with requests.get(url, stream=True, allow_redirects=True) as r, \
                open(filename, "wb") as f, \
                tqdm.tqdm(
                    unit="B",  # unit string to be displayed.
                    unit_scale=True,  # let tqdm to determine the scale in kilo, mega..etc.
                    unit_divisor=1024,  # is used when unit_scale is true
                    total=filesize,  # the total iteration.
                    file=sys.stdout,  # default goes to stderr, this is the display on console.
                    desc=filename  # prefix to be displayed on progress bar.
                    ) as progress:
            for chunk in r.iter_content(chunk_size=chunk_size):
                # download the file chunk by chunk
                datasize = f.write(chunk)
                progress.update(datasize)

        with tarfile.open(f'{self.FOLDER}/MIST_v1.2_vvcrit0.4_{magsys_name}.txz') as tf:
            tf.extractall(self.FOLDER)
        os.rename(f'{self.FOLDER}/MIST_v1.2_vvcrit0.4_{magsys_name}{suffix}',
            f'{self.FOLDER}/{magsys_name}')
        os.remove(f'{self.FOLDER}/MIST_v1.2_vvcrit0.4_{magsys_name}.txz')

        print('Isochrones downloaded: %s' % magsys_name)
