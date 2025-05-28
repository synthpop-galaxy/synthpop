"""
This module provides the loading and handling of SPISEA isochrones.
It automatically generates isochrone from the SPISEA package.
"""
__all__ = ["Spisea",]
__author__ = "M.J. Huston"
__date__ = "2025-05-14"

import glob
import importlib
import linecache
import os
import tarfile
import warnings
import json
import tqdm
import sys
from spisea import synthetic, reddening, atmospheres, evolution
from .. import const
import pdb

import numpy as np
import pandas
import requests

from ._evolution import EvolutionIsochrones, ISOCHRONES_DIR, EVOLUTION_DIR
# import a "standard" interpolator
from .charon_interpolator import CharonInterpolator
#from .lagrange_interpolator import LagrangeInterpolator

# global variable to store the isochrones
spisea_isochrones = None
spisea_columns = {}

class SpiseaIso(EvolutionIsochrones, CharonInterpolator):
    """
    SPISEA Isochrone class
    """

    # folder where isochrone files can be found
    FOLDER = f"{ISOCHRONES_DIR}/spisea"

    # lowest and highest mass where MIST isochrones should be used
    # (can be outside the covered range, in such cases the closets grid_points are used)
    # used area can be overwritten in the config.json file
    max_mass = 300
    min_mass = 0.1
    isochrones_name = 'SPISEA'

    def __init__(self, columns, use_global=True, evolution_model_name='MISTv1.2', **kwargs):
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

        self.allowed_non_mag_cols = ["[Fe/H]_init", "log10_isochrone_age_yr", 'phase',
            'star_mass', 'initial_mass', 'log_L', 'log_R', 'log_Teff', 'log_g']
        with open(f"{EVOLUTION_DIR}/spisea_filters.json") as f:
            self.spisea_filters = json.load(f)

        self.magsys, self.non_mag_cols, self.magsys_bands = self.get_cols(columns)
        self.evolution_model_name = evolution_model_name

        if (evolution_model_name=='MISTv1.2') or (evolution_model_name=='MISTv1.0'):
            self.evolution_model = evolution.MISTv1(version=float(evolution_model_name[-3:]))
            self.feh_grid = np.array([-4.0,-3.5,-3.0,-2.5,-2.0,-1.75,-1.5,-1.25,
                                      -1.0,-0.75,-0.5,-0.25,0,0.25,0.5])
            self.log_age_list = np.linspace(5.0,10.3,107)
            self.log_age_list[0] = 5.01
            assert len(self.feh_grid)==len(self.evolution_model.z_list)
        else:
            raise ValueError("Invalid SPISEA evolution_model. Only MISTv1.0 and MISTv1.2 are available at this time.")
        self.atm_func = atmospheres.get_merged_atmosphere
        self.red_law = reddening.RedLawHosek18b() # Doesn't matter - we do 0 extinction here and apply it later

        # load the isochrones
        if use_global:
            global spisea_isochrones
            if spisea_isochrones is not None:
                if np.array_equal(np.sort(np.unique(columns)).astype(str), np.sort(np.unique(spisea_isochrones.keys())).astype(str)):
                    self.isochrones = spisea_isochrones
                else:
                    spisea_isochrones = self.isochrones = self.load_isochrones()
            else:
                spisea_isochrones = self.isochrones = self.load_isochrones()
        else:
            spisea_isochrones = self.isochrones = self.load_isochrones()
        self.isochrones_grouped = self.isochrones.groupby(["[Fe/H]_init", "log10_isochrone_age_yr"])

        # Get mass range for each metallicity and age
        self.mass_range = self.get_mass_ranges(self.isochrones_grouped)

        #pdb.set_trace()

        # call super after loading the isochrones so the interpolator has axcess to the data
        super().__init__(**kwargs)

    def get_filename(self, magsys, met_string):
        filename_pattern = f"{self.FOLDER}/{magsys}/" \
                           f"{self.evolution_model_name}_feh_{met_string}.iso"
        return filename_pattern 

    def get_cols(self, columns):
        iso_files = {}
        non_mag_cols = []

        for column in columns:
            if column in self.allowed_non_mag_cols:
                non_mag_cols.append(column)
            else:
                try:
                    column_split = column.split('-')
                    if len(column_split)==2:
                        msys = column_split[0]
                    elif len(column_split)==3:
                        msys = column_split[0]+'-'+column_split[1]
                    filt = column_split[-1]
                    if (msys not in iso_files):
                        iso_files[msys] = []
                    iso_files[msys].append(column)
                except:
                    raise ValueError('Invalid column '+column+' for SPISEA isochrones.')

        return list(iso_files.keys()), list(np.unique(non_mag_cols)), iso_files

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

    def generate_isochrones(self, magsyses, feh):
        """
        generate new isochrone from spisea, convert into our format, and save
        """
        print(f'  Generating isochrone for requested filter systems in feh={feh}')
        i=np.where(self.feh_grid==feh)[0][0]
        m_h = np.log10(self.evolution_model.z_list[i] / self.evolution_model.z_solar)
        isos_tmp = []

        # Go ahead and generate isochrone for full filter system now, 
        #  so we don't have to recompute them again in the future.
        calc_filters = []
        for magsys in magsyses:
            msys_split = magsys.split('-')
            if len(msys_split)==1:
                filts = self.spisea_filters[msys_split[0]]
                calc_filters += [msys_split[0]+','+filt for filt in filts]
            else:
                filts = self.spisea_filters[msys_split[0]][msys_split[1]]
                calc_filters += [msys_split[0]+','+msys_split[1]+','+filt for filt in filts]

        for log_age in self.log_age_list:
            iso_tmp = synthetic.IsochronePhot(log_age, 0.0, 10.0, metallicity=m_h,
                            evo_model=self.evolution_model,
                            atm_func=self.atm_func,
                            red_law=self.red_law,
                            filters=calc_filters,
                            iso_dir=self.FOLDER+'/tmp/',
                            recomp=True).points
            iso_tmp['log10_isochrone_age_yr']=log_age
            isos_tmp.append(iso_tmp.to_pandas())
        isos_cat_tmp = pandas.concat(isos_tmp)
        iso = isos_cat_tmp[['log10_isochrone_age_yr', 'phase']]
        iso.insert(2, '[Fe/H]_init', feh)
        iso.insert(3, 'star_mass', isos_cat_tmp['mass_current'])
        iso.insert(4, 'initial_mass', isos_cat_tmp['mass'])
        iso.insert(5, 'log_L', np.log10(isos_cat_tmp['L']/const.Lsun_w))
        iso.insert(6, 'log_R', np.log10(isos_cat_tmp['R']/const.Rsun_m))
        iso.insert(7, 'log_Teff', np.log10(isos_cat_tmp['Teff']))
        iso.insert(8, 'log_g', isos_cat_tmp['logg'])
        ccount = 9
        for col in calc_filters:
            iso.insert(ccount, col.replace(',','-'), isos_cat_tmp['m_'+synthetic.get_filter_col_name(col)])
            ccount+=1
        #pdb.set_trace()
        return iso


    def load_isochrones(self):
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

        # isochrones = {}
        # # Load columns
        # for i, magsys_name in enumerate(self.magsys):
        #     use_columns = list(self.magsys_bands[magsys_name])
        #     # Get basic properties from first filter set isochrone
        #     if i==0:
        #         use_columns += self.non_mag_cols
        #     if not os.path.isdir(f'{self.FOLDER}/{magsys_name}'):
        #         print(f"Missing SPISEA {self.evolution_model_name} isochrones for {self.magsys_name}.")
        #         print(f"Generating these may take a few hours the first time, then they will be saved for quick loading in the future.")
        #         os.makedirs(f'{self.FOLDER}/{magsys_name}')

        #     # Load in columns for each metallicity
        #     for file_met in self.feh_grid:
        #         met_string = ('m' if file_met<0.0 else 'p') + f'{np.abs(file_met):1.2f}'
        #         filename = self.get_filename(magsys_name, met_string)

        #         if os.path.isfile(f'{filename}.h5'):
        #             df = pandas.read_hdf(f'{filename}.h5', key='data')
        #         else:
        #             df = self.generate_isochrones(magsys_name, file_met)
        #             df.to_hdf(f'{filename}.h5', key='data')

        #         if i == 0:
        #             isochrones[file_met] = df[use_columns].copy()
        #         else:
        #             isochrones[file_met][use_columns] = df[use_columns]
        #pdb.set_trace()

        need_to_generate = {}
        isochrones = {}

        # Check which isochrones need generated
        for file_met in self.feh_grid:
            met_string = ('m' if file_met<0.0 else 'p') + f'{np.abs(file_met):1.2f}'
            for magsys_name in self.magsys:
            # Check whether we have a file saved for this filter set
                filename = self.get_filename(magsys_name, met_string)
                if not os.path.isfile(f'{filename}.h5'):
                    if file_met not in need_to_generate.keys():
                        need_to_generate[file_met] = [magsys_name]
                    else:
                        need_to_generate[file_met].append(magsys_name)
            # Generate needed isochrones for given metallicity
            if file_met in need_to_generate.keys():
                df = self.generate_isochrones(need_to_generate[file_met], file_met)
                for magsys_name in need_to_generate[file_met]:
                    os.makedirs(f'{self.FOLDER}/{magsys_name}', exist_ok=True)
                    filename = self.get_filename(magsys_name, met_string)
                    msys_split = magsys_name.split('-')
                    if len(msys_split)==1:
                        filts = self.spisea_filters[msys_split[0]]
                        use_filters = [msys_split[0]+'-'+filt for filt in filts]
                    else:
                        filts = self.spisea_filters[msys_split[0]][msys_split[1]]
                        use_filters = [msys_split[0]+'-'+msys_split[1]+'-'+filt for filt in filts]
                    use_columns = list(self.non_mag_cols) + list(use_filters)
                    df[use_columns].to_hdf(f'{filename}.h5', key='data')

        # Load or generate columns
        for file_met in self.feh_grid:
            met_string = ('m' if file_met<0.0 else 'p') + f'{np.abs(file_met):1.2f}'

            for i, magsys_name in enumerate(self.magsys):
                use_columns = list(self.magsys_bands[magsys_name])
                # Get basic properties from first filter set isochrone
                if i==0:
                    use_columns += self.non_mag_cols
                df = pandas.read_hdf(filename, key='data')

                if i == 0:
                    isochrones[file_met] = df[use_columns].copy()
                else:
                    isochrones[file_met][use_columns] = df[use_columns]

        return pandas.concat(isochrones.values())


