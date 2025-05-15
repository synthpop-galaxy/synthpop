"""
This module provides the loading and handling of MIST isochrones.
It automatically downloads the isochrone form the online MIST service.

MIST documentation: https://waps.cfa.harvard.edu/MIST/

MIST papers: Dotter (2016), Choi et al. (2016);
DOIs: 10.3847/0067-0049/222/1/8, 10.3847/0004-637X/823/2/102

The mist_columns.json file provides a guide to which isochrone files include which columns/filters.
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
import spisea
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

    def __init__(self, columns, use_global=True, evolution_model_name='MISTv1', **kwargs):
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

        self.non_mag_cols, self.filters = self.get_cols(columns)

        if evolution_model=='MISTv1':
            self.evolution_model = spisea.evolution.MISTv1()
            self.feh_grid = np.array([-4.0,-3.5,-3.0,-2.5,-2.0,-1.75,-1.5,-1.25,
                                      -1.0,-0.75,-0.5,-0.25,0,0.25,0.5])
            self.log_age_list = np.linspace(5.0,10.3,107)
            assert len(feh_grid)==len(self.evolution_model.z_list)
        else:
            raise ValueError("Invalid SPISEA evolution_model. Only MISTv1 is available at this time.")
        self.atm_func = spisea.atmospheres.get_merged_atmosphere
        self.red_law = spisea.reddening.RedLawHosek18b()

        self.isochrones = self.load_isochrones()
        self.isochrones_grouped = self.isochrones.groupby(["[Fe/H]_init", "log10_isochrone_age_yr"])

        # Get ages
        self.iso_ages = 10 ** self.age_grid

        # Get mass range for each metallicity and age
        self.mass_range = self.get_mass_ranges(self.isochrones_grouped)

        # call super after loading the isochrones so the interpolator has axcess to the data
        super().__init__(**kwargs)

    @staticmethod
    def get_cols(columns):


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

    def load_isochrones(self, magsys):
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
        # Load columns
        for i,feh in enumerate(self.feh_grid):
            isos_tmp = []
            m_h = np.log10(self.evolution_model.z_list[i] / self.evolution_model.z_solar)
            for j,log_age in enumerate(self.log_age_list):
                iso_tmp = spisea.synthetic.IsochronePhot(log_age, 0.0, 10.0, metallicity=m_h,
                                evo_model=self.evolution_model,
                                atm_func=self.atm_func,
                                red_law=self.red_law,
                                filters=self.filters,
                                iso_dir=self.FOLDER+'/tmp/').points
                isos_tmp.append(iso_tmp.to_pandas())
                isos_tmp[-1]['log10_isochrone_age_yr']=log_age
            isos_cat_tmp = pandas.concat(isos_tmp)
            isochrones[feh] = isos_cat_tmp['log10_isochrone_age_yr', 'phase']
            isochrones[feh]['[Fe/H]_init'] = feh
            isochrones[feh]['star_mass'] = isos_cat_tmp['current_mass']
            isochrones[feh]['initial_mass'] = isos_cat_tmp['mass']
            isochrones[feh]['log_L'] = np.log10(isos_cat_tmp['L']/const.Lsun_w)
            isochrones[feh]['log_R'] = np.log10(isos_cat_tmp['R']/const.Rsun_m)
            isochrones[feh]['log_Teff'] = np.log10(isos_cat_tmp['Teff'])
            isochrones[feh]['log_g'] = isos_cat_tmp['logg']

        return pandas.concat(isochrones.values())


