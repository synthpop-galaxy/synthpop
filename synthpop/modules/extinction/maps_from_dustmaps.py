""" Extinction maps collected from the dustmap module

Using dustmaps allows to evaluate the dustmap for each star. 
which allows accurate results for larger cones.

Before using this module you need to fatch the data once.
See https://dustmaps.readthedocs.io/en/latest/installation.html

"""

__all__ = ["MapsFromDustmaps",]
__author__ = ["J. Klüter", "M.J. Huston"]
__date__ = "2022-11-05"
__license__ = "GPLv3"
__version__ = "1.0.0"

import os
import astropy.units as u

import dustmaps.bayestar
import dustmaps.bh
import dustmaps.chen2014
import dustmaps.csfd
import dustmaps.gaia_tge
import dustmaps.iphas
import dustmaps.leike_ensslin_2019
import dustmaps.leike2020
import dustmaps.edenhofer2023
import dustmaps.lenz2017
import dustmaps.marshall
import dustmaps.pg2010
import dustmaps.planck
import dustmaps.sfd

from astropy.coordinates import SkyCoord
try:
    from ._extinction import ExtinctionMap
except ImportError:
    from _extinction import ExtinctionMap
    
MAPS_INFO = {
    "bayestar": {
        'dim': 3,
        'returns': 'E(B-V)',
        "Query": dustmaps.bayestar.BayestarQuery,
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "options": {"max_samples": 0},
        "kwargs": {"mode":'best'}},
    "bh": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.bh.BHQuery,
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551},
    "chen2014": {
        'dim': 3,
        'returns': 'A_r',
        "Query": dustmaps.chen2014.Chen2014Query,
        "lambda_eff": 0.622},
    "csfd": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.csfd.CSFDQuery,
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551},
    "gaia_tge":{
        'dim': 2,
        'returns': 'A0',
        "Query": dustmaps.gaia_tge.GaiaTGEQuery,
        "lambda_eff": 0.5414},
    "iphas": {
        'dim': 3,
        'returns': 'A0',
        "Query": dustmaps.iphas.IPHASQuery,
        "lambda_eff": 0.5495},
    "leike_ensslin_2019": {
        'dim': 3,
        'returns': 'e-foldings_GaiaG',
        "Query": dustmaps.leike_ensslin_2019.LeikeEnsslin2019Query,
        "lambda_eff": 0.673},
    "leike_2020": {
        'dim': 3,
        'returns': 'e-foldings_GaiaG',
        "Query": dustmaps.leike2020.Leike2020Query,
        "lambda_eff": 0.673},
    "edenhofer2023": {
        'dim': 3,
        'returns': 'E(B-V)',
        "Query": dustmaps.edenhofer2023.Edenhofer2023Query,
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "options": {"integrated":True}},
    "lenz2017": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.lenz2017.Lenz2017Query,
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551},
    "marshall": {
        'dim': 3,
        'returns': 'A_Ks',
        "Query": dustmaps.marshall.MarshallQuery,
        "lambda_eff": 2.152},
    "pg2010": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.pg2010.PG2010Query,
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551},
    "planck_gnlc": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.planck.PlanckGNILCQuery,
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551},
    "planck": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.planck.PlanckQuery,
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551},
    "sfd": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.sfd.SFDQuery,
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551}
    }

# empty dictionary to store dustmaps query
_query_dict = {}

class MapsFromDustmaps(ExtinctionMap):
    def __init__(self, dustmap_name=None, dist_2d = 0.0, **kwargs):
        super().__init__(**kwargs)
        #pre loaded query functions to share between populations.
        global _query_dict
        self.extinction_map_name = f"dustmaps.{dustmap_name}"

        if dustmap_name is None:
            raise ValueError("dustmap_name needs to be defined")

        if dustmap_name not in MAPS_INFO.keys():
            raise ValueError(f"{dustmap_name} does not exist in dustmap_name\n If it exists in the dustmaps module, please submit an issue or pull request to enable it in SynthPop")
        if dustmap_name.startswith("leike"):
            raise NotImplementedError('leike2020 & leike_ensslin_2019 are not implemented yet')
        map_props = MAPS_INFO[dustmap_name]

        if not os.path.isdir(os.path.join(dustmaps.std_paths.data_dir(), dustmap_name)):
            url = 'https://dustmaps.readthedocs.io/en/latest/installation.html'
            module = map_props["Query"].__module__
            
            try:
                print("Downloading",self.extinction_map_name)
                getattr(dustmaps, dustmap_name).fetch()
            except:
                raise FileNotFoundError(
                    f"Data for '{dustmap_name} dustmap are not fetched\n"
                    f"when a dustmaps data directory is specified, data can be fetched using:\n\n"
                    f">>> import {module}\n"
                    f">>> {module}.fetch()\n\n"
                    f"please see '{url}' for further details")

        self.is_3D = map_props['dim'] == 3
        if not self.is_3D:
            self.dist_2d = dist_2d

        # select the query object for the given dustmap
        if dustmap_name not in _query_dict:
            _query_dict[dustmap_name] = map_props['Query'](**map_props.get('options', {}))

        self.kwargs = map_props.get('kwargs', {})
        self.query = _query_dict[dustmap_name]
        self.ref_wavelength = map_props['lambda_eff']
        self.A_or_E_type = map_props['returns']
        if self.A_or_E_type.startswith("E"):
            self.ref_wavelength2 = map_props['lambda_eff2']

        self.bands = []  # list of filters
        self.eff_wavelengths = {}  # effective wavelength for each band

    def extinction_in_map(self, l_deg, b_deg, dist):
        """
        Estimates the extinction for a list of star positions.

        Parameters
        ----------
        l_deg: ndarray [degrees]
            galactic longitude
        b_deg: ndarray [degrees]
            galactic latitude
        dist: ndarray [kpc]
            radial distance from the Sun
        
        Returns
        -------
        extinction_value: ndarray [mag]
            extinction at each star position defined as self.A_or_E_type
        """
        if self.is_3D:
            coords = SkyCoord(l_deg * u.deg, b_deg * u.deg, distance=dist * u.kpc, frame='galactic')
            dist_factor = 1.0
        else:
            coords = SkyCoord(l_deg * u.deg, b_deg * u.deg, frame='galactic')
            dist_factor = (dist>self.dist_2d)
        return self.query(coords, **self.kwargs)*dist_factor
