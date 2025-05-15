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
        #with open(f"{EVOLUTION_DIR}/spisea_filters.json") as f:
        #    self.spisea_filters = json.load(f)

        self.non_mag_cols, self.mag_cols = self.get_cols(columns)
        self.evolution_model_name = evolution_model_name

        if (evolution_model_name=='MISTv1.2') or (evolution_model_name=='MISTv1.0'):
            self.evolution_model = evolution.MISTv1(version=float(evolution_model_name[-3:]))
            self.feh_grid = np.array([-4.0,-3.5,-3.0,-2.5,-2.0,-1.75,-1.5,-1.25,
                                      -1.0,-0.75,-0.5,-0.25,0,0.25,0.5])
            self.log_age_list = np.linspace(5.0,10.3,107)
            self.log_age_list[0] = 5.01
            print(np.log10(np.array(self.evolution_model.z_list) / self.evolution_model.z_solar))
            assert len(self.feh_grid)==len(self.evolution_model.z_list)
        else:
            raise ValueError("Invalid SPISEA evolution_model. Only MISTv1.0 and MISTv1.2 are available at this time.")
        self.atm_func = atmospheres.get_merged_atmosphere
        self.red_law = reddening.RedLawHosek18b() # Doesn't matter - we do 0 extinction here and apply it later

        # load the isochrones
        if use_global:
            global spisea_isochrones
            if spisea_isochrones is not None:
                self.isochrones = spisea_isochrones
            else:
                self.isochrones = self.load_isochrones()
        else:
            self.isochrones = self.load_isochrones()
        self.isochrones_grouped = self.isochrones.groupby(["[Fe/H]_init", "log10_isochrone_age_yr"])

        # Get mass range for each metallicity and age
        self.mass_range = self.get_mass_ranges(self.isochrones_grouped)

        # call super after loading the isochrones so the interpolator has axcess to the data
        super().__init__(**kwargs)

    def get_cols(self, columns):
        #iso_files = {'basic':[]}

        # for column in columns:
        #     if column in self.allowed_non_mag_cols:
        #         iso_files['basic'].append(column)
        #     else:
        #         try:
        #             column_split = column.split(',')
        #             if len(column_split)==2:
        #                 filt = column_split[1]
        #             elif len(column_split)==3:
        #                 filt = column_split[1]+'_'+column_split[2]
        #             msys = column_split[0]
        #             if ~(msys in iso_files):
        #                 iso_files[msys] = []
        #             iso_files[msys].append(filt)
        #         except:
        #             raise ValueError('Invalid column '+column+' for SPISEA isochrones.')

        non_mag_cols = []
        mag_cols = []

        for column in columns:
             if column in self.allowed_non_mag_cols:
                 non_mag_cols.append(column)
             else:
                mag_cols.append(column)
        #         try:
        #             column_split = column.split(',')
        #             if len(column_split)==2:
        #                 filt = column_split[1]
        #             elif len(column_split)==3:
        #                 filt = column_split[1]+'_'+column_split[2]
        #             msys = column_split[0]
        #             if ~(msys in iso_files):
        #                 iso_files[msys] = []
        #             iso_files[msys].append(filt)
        #         except:
        #             raise ValueError('Invalid column '+column+' for SPISEA isochrones.')

        return non_mag_cols, mag_cols

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

    def generate_isochrones(self, feh, m_h):
        """
        generate new isochrone from spisea, convert into our format, and save
        """

        #feh_str = ('m' if feh<0.0 else 'p') + f'{np.abs(feh):1.2f}'
        isos_tmp = []
        calc_filters = []
        for col in self.mag_cols:
            calc_filters.append(col.replace('-',','))
        for log_age in self.log_age_list:
            print(m_h)
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
        iso = isos_cat_tmp['log10_isochrone_age_yr', 'phase']
        iso['[Fe/H]_init'] = feh
        iso['star_mass'] = isos_cat_tmp['current_mass']
        iso['initial_mass'] = isos_cat_tmp['mass']
        iso['log_L'] = np.log10(isos_cat_tmp['L']/const.Lsun_w)
        iso['log_R'] = np.log10(isos_cat_tmp['R']/const.Rsun_m)
        iso['log_Teff'] = np.log10(isos_cat_tmp['Teff'])
        iso['log_g'] = isos_cat_tmp['logg']
        for col in self.mag_cols:
            iso[col] = isos_cat_tmp[col.replace('-',',')]
        return iso

    def load_isochrones(self):
        """
        load the isochrones for each of the given magnitude systems and all metallicities

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
        #new_cols_needed = []
        
        for i,feh in enumerate(self.feh_grid):
            #feh_str = ('m' if nn<0.0 else 'p') + f'{np.abs(feh):1.2f}'
            #for iso_file in iso_files:
                # os.mkdirs(self.FOLDER+'/'+iso_file, exist_ok=True)
                # fname = self.FOLDER+'/'+iso_file+f'/'+self.evolution_model_name+'_iso_{feh_str}.h5'
                # if os.isfile(fname):
                #     iso_tmp = pandas.read_hdf(fname, key='iso')
                #     for col in iso_files[iso_file]:
                #         if ~(col in isochrones[feh]):
                #             new_cols_needed.append(col)
                # else: 
            m_h = np.log10(self.evolution_model.z_list[i] / self.evolution_model.z_solar)
            isochrones[feh] = self.generate_isochrones(feh, m_h)
                #     iso_tmp.to_hdf(fname, key='iso')


        return pandas.concat(isochrones.values())


