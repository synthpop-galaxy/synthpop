""" Extinction maps collected from the dustmap module
Currently only marshall is finished since the others need a transformation between eg E(B-V) and AV
to be comparable with the extinction laws.

Using dustmaps allows to evaluate the dustmap for each star. 
which allows accurate results for larger cones.

"""

__all__ = ["MapsFromDustmaps",]
__author__ = "J. Klüter"
__date__ = "2022-11-05"
__license__ = "GPLv3"
__version__ = "1.0.0"

import astropy.units as u

import dustmaps.sfd
import dustmaps.planck
import dustmaps.bayestar
import dustmaps.iphas
import dustmaps.marshall
import dustmaps.chen2014
import dustmaps.lenz2017
import dustmaps.pg2010
import dustmaps.leike_ensslin_2019
import dustmaps.leike2020

from astropy.coordinates import SkyCoord
from ._extinction import ExtinctionMap

MAPS_INFO = {
    "sfd": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.sfd.SFDQuery,
        "lambda_eff": 0.551},
    "planck": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.planck.PlanckQuery,
        "lambda_eff": 0.551},
    "bayestar": {
        'dim': 3,
        'returns': 'E(B-V)',
        "Query": dustmaps.bayestar.BayestarQuery,
        "lambda_eff": 0.551,
        "options": {"max_samples": 1}},
    "iphas": {
        'dim': 3,
        'returns': 'A0',
        "Query": dustmaps.iphas.IPHASQuery,
        "lambda_eff": 0.5495},
    "marshall": {
        'dim': 3,
        'returns': 'A_Ks',
        "Query": dustmaps.marshall.MarshallQuery,
        "lambda_eff": 2.152},
    "chen2014": {
        'dim': 3,
        'returns': 'A_r',
        "Query": dustmaps.chen2014.Chen2014Query,
        "lambda_eff": 0.622},
    "lenz2017": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.lenz2017.Lenz2017Query,
        "lambda_eff": 0.551},
    "pg2010": {
        'dim': 2,
        'returns': 'E(B-V)',
        "Query": dustmaps.pg2010.PG2010Query,
        "lambda_eff": 0.551},
    "leike_ensslin_2019": {
        'dim': 3,
        'returns': 'e-foldings_GaiaG',
        "Query": dustmaps.leike_ensslin_2019.LeikeEnsslin2019Query,
        "lambda_eff": 0.673},
    "leike2020": {
        'dim': 3,
        'returns': 'e-foldings_GaiaG',
        "Query": dustmaps.leike2020.Leike2020Query,
        "lambda_eff": 0.673}
    }

# empty dictionary to store dustmaps query
_query_dict = {}


class MapsFromDustmaps(ExtinctionMap):
    def __init__(self, dustmap_name=None, return_functions=True, **kwargs):
        global _query_dict
        self.extinction_map_name = f"dustmaps.{dustmap_name}"
        self.return_functions = return_functions

        if dustmap_name is None:
            raise ValueError("dustmap_name needs to be defined")

        if dustmap_name not in MAPS_INFO.keys():
            raise ValueError(f"{dustmap_name} does not exist in dustmap_name")
        if dustmap_name.startswith("leike"):
            raise NotImplementedError('leike2020 & leike_ensslin_2019 are not implemented yet')
        map_props = MAPS_INFO[dustmap_name]

        self.is_3D = map_props['dim'] == 3
        if not self.is_3D:
            raise NotImplementedError('2D-maps are not implemented yet')

        # select the query object for the given dustmap
        if dustmap_name not in _query_dict:
            _query_dict[dustmap_name] = map_props['Query'](**map_props.get('options', {}))

        self.query = _query_dict[dustmap_name]
        self.ref_wavelength = map_props['lambda_eff']
        self.A_or_E_type = map_props['returns']

        # placeholder for location, filter properties, etc.
        self.l_deg = None
        self.b_deg = None

        self.bands = []  # list of filters
        self.eff_wavelengths = {}  # effective wavelength for each band

    def get_map(self, l_deg, b_deg, dist):
        """
        read extinction from map


        Parameters
        ----------
        l_deg: float or nd_array [degree]
            galactic longitude
        b_deg: float or nd_array [degree]
            galactic latitude
        dist: float or nd_array [kpc]
            distance

        Returns
        -------
        A_or_E: float or nd_array
            Extinction ore color excess
        """
        if self.is_3D:
            coords = SkyCoord(l_deg * u.deg, b_deg * u.deg, distance=dist * u.kpc, frame='galactic')
        else:
            coords = SkyCoord(l_deg * u.deg, b_deg * u.deg, frame='galactic')
        return self.query(coords)

    def update_extinction_in_map(self, radius):
        """
        Estimates the extinction for the current sight-line and radial distance
        store the result into self.extinction_in_map.

        Parameters
        ----------
        radius: float [kpc]
            radial distance of the current slice
        """

        if self.return_functions:
            self.extinction_in_map = self.get_map
        else:
            self.extinction_in_map = self.get_map(self.l_deg, self.b_deg, radius)
