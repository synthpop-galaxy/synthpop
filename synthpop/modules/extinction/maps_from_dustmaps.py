"""
Extinction maps collected from the dustmap module (Green, 2018)

Extinction form (i.e. reference wavelength / color excess) depends on
map used.

Publication DOI: 10.21105/joss.00695

Before using this module you need to fatch the data once.
See https://dustmaps.readthedocs.io/en/latest/installation.html
"""

__all__ = ["MapsFromDustmaps",]
__author__ = ["J. Klüter", "M.J. Huston"]
__date__ = "2022-11-05"

import os
import astropy.units as u
import importlib

import dustmaps

from astropy.coordinates import SkyCoord
try:
    from ._extinction import ExtinctionMap
except ImportError:
    from _extinction import ExtinctionMap
    
MAPS_INFO = {
    "bayestar": {
        'dim': 3,
        'returns': 'E(B-V)',
        "Query": "BayestarQuery",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "options": {"max_samples": 0},
        "kwargs": {"mode":'best'},
        "need_fetch":True},
    "bh": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": "BHQuery",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "need_fetch":False},
    "chen2014": {
        'dim': 3,
        'returns': 'A_r',
        "Query": "Chen2014Query",
        "lambda_eff": 0.622,
        "need_fetch":True},
    "csfd": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": "CSFDQuery",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "need_fetch":True},
    "decaps": {
        'dim': 3,
        'returns': 'E(B-V)',
        "Query": "DECaPSQueryLite",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "need_fetch":True},
    "gaia_tge":{
        'dim': 2,
        'returns': 'A0',
        "Query": "GaiaTGEQuery",
        "lambda_eff": 0.5414,
        "need_fetch":True},
    "iphas": {
        'dim': 3,
        'returns': 'A0',
        "Query": "IPHASQuery",
        "lambda_eff": 0.5495,
        "need_fetch":True},
    "leike_ensslin_2019": {
        'dim': 3,
        'returns': 'e-foldings_GaiaG',
        "Query": "LeikeEnsslin2019Query",
        "lambda_eff": 0.673,
        "need_fetch":True},
    "leike_2020": {
        'dim': 3,
        'returns': 'e-foldings_GaiaG',
        "Query": "Leike2020Query",
        "lambda_eff": 0.673,
        "need_fetch":True},
    "edenhofer2023": {
        'dim': 3,
        'returns': 'E(B-V)',
        "Query": "Edenhofer2023Query",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "options": {"integrated":True},
        "need_fetch":True},
    "lenz2017": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": "Lenz2017Query",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "need_fetch":True},
    "marshall": {
        'dim': 3,
        'returns': 'A_Ks',
        "Query": "MarshallQuery",
        "lambda_eff": 2.152,
        "need_fetch":True},
    "pg2010": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": "PG2010Query",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "need_fetch":True},
    "planck_gnlc": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": "PlanckGNILCQuery",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "need_fetch":True},
    "planck": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": "PlanckQuery",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "need_fetch":True},
    "sfd": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": "SFDQuery",
        "lambda_eff": 0.493,
        "lambda_eff2": 0.551,
        "need_fetch":True}
    }

# empty dictionary to store dustmaps query
_query_dict = {}

class MapsFromDustmaps(ExtinctionMap):
    """
    Extinction maps from dustmaps

    Attributes
    ----------
    dustmap_name : string
        name of the map from dustmaps to use.
        options are "bayestar", "bh", "chen2014", "csfd", "decaps",
        "gaia_tge", "iphas", "leike_ensslin_2019", "leike_2020",
        "edenhofer2023", "lenz2017", "marshall", "pg2010",
        "planck_gnlc", "planck", "sfd"
    dist_2d=0.0 : float
        if a 2-d map is applied, this sets the distance [kpc] where
        the extinction is applied as a single infinitely thin screen

    Methods
    -------
    extinction_in_map(l_deg, b_deg, dist):
        estimates extinction for a list of star positions
    """

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
        
        map_module = importlib.import_module('dustmaps.'+dustmap_name)

        if map_props['need_fetch'] and (not os.path.isdir(os.path.join(dustmaps.std_paths.data_dir(), dustmap_name))):
            url = 'https://dustmaps.readthedocs.io/en/latest/installation.html'
            
            try:
                print("Downloading",self.extinction_map_name)
                map_module.fetch()
            except:
                raise FileNotFoundError(
                    f"Data for '{dustmap_name}' dustmap are not fetched\n"
                    f"when a dustmaps data directory is specified, data can be fetched using:\n\n"
                    f">>> import dustmaps.{dustmap_name}\n"
                    f">>> dustmaps.{dustmap_name}.fetch()\n\n"
                    f"please see '{url}' for further details")

        self.is_3D = map_props['dim'] == 3
        if not self.is_3D:
            self.dist_2d = dist_2d

        # select the query object for the given dustmap
        if dustmap_name not in _query_dict:
            _query_dict[dustmap_name] = getattr(map_module, map_props['Query'])(**map_props.get('options', {}))

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
