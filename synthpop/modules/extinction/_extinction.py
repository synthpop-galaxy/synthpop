"""
This module provides the parent classes for the Extinction:
An Extinction Class is a child of an ExtinctionMap and one or more ExtinctionLaws.
Those will be  combined with  the CombineExtinction factory functions.

ExtinctionLaw: Parent class for extinction laws
ExtinctionMap: Parent class for extinction maps

CombineExtinction: Factory that combines an ExtinctionLaw subclass
                   and an ExtinctionMap subclass into one Extinction Class
"""

__all__ = ["ExtinctionLaw", "ExtinctionMap", "CombineExtinction", "EXTINCTION_DIR"]
__author__ = "J. Klüter, M.J. Huston, S. Johnson"
__credits__ = ["J. Klüter", "M.J. Huston", "S. Johnson", "A. Aronica", "M. Penny"]
__date__ = "2022-07-09"

import os
import inspect
from typing import Callable, Tuple, List, Dict
from types import ModuleType
from .. import const
from abc import ABC, abstractmethod
import numpy as np

EXTINCTION_DIR = const.EXTINCTIONS_DIR

class ExtinctionLaw(ABC):
    """
    Parent class for an extinction law.
    It is enough to provide only a function for Alambda_Aref within the Subclasses.
    However, only extinction_at_lambda will be called by SynthPop.
    So a complete different handling of the extinction can be implemented

    Attributes
    ----------
    extinction_law_name: str
        name of extinction law
    law_ref_wavelength: float
        reference wavelength for the extinction law in microns
        lambda where Alambda_Aref = 1
    min_wavelength: float
        lower limit of valid wavelengths for the law given in microns
    max_wavelength: float
        upper limit of valid wavelengths for the law given in microns
        
    map_name: str
        name of the extinction map, when combined into an Extinction Object.
        will be set by the set_map_properties function
    ref_wavelength : float
        reference wavelength for the extinction map in microns
        for an A_lambda map, lambda
        for an E(lambda1-lambda2) map, lambda1
    ref_wavelength2 : float
        secondary reference wavelength for the extinction map in microns
        for an A_lambda map, none
        for an E(lambda1-lambda2) map, lambda2
    A_or_E_type : str
        specify if an extinction Alambda_ref, or a color excess is given by A_or_E
        if it starts with "A": A_or_E is a total extinction
        if it starts with "E": A_or_E is a color excess
        will be set by the set_map_properties function

    Methods
    ------- 
    Alambda_Aref(eff_wavelength)
        returns the absorption coefficient for a given effective wavelength,
        relative to the reference extinction.

    Alambda_Amap(eff_wavelength)
        wrapper such that the returned absorption is relative to the
        reference wavelength of the Extinction map

    extinction_at_lambda(eff_wavelength, A_or_E)
        estimates the total extinction at a given effective wavelength,
        and a given total extinction or color excess from the extinction map

    set_map_properties():
        method to set the important properties from the extinction map.
        This is the reference Wavelength, and the Extinction type
        (color excess or total extinction)
    """

    def __init__(self, logger: ModuleType = None,  **kwargs):
        self.logger =logger
        # Extinction law parameters
        self.extinction_law_name = "NONE"
        self.law_ref_wavelength = None
        self.min_wavelength = None 
        self.max_wavelength = None
        # Parameters to link this to the extinction map
        self.map_name = "NONE"
        self.ref_wavelength = None
        self.A_or_E_type = None
        self.ref_wavelength2 = None

    @abstractmethod
    def Alambda_Aref(self, eff_wavelength: float) -> float:
        """
        estimate A_lambda / A_ref,

        Parameters
        ----------
        eff_wavelength : float
            Effective Wavelength of the filter for which the extinction should be determined.
            in micrometer

        """
        raise NotImplementedError('No extinction law set')

    def Alambda_Amap(self, eff_wavelength: float) -> float:
        """
        Calculate the extinction relative to the specified filter
        for a given effective wavelength.
                 
        Arguments
        ---------
        eff_wavelength : float
            Effective Wavelength of the filter for which the extinction should be determined.
            in micrometer
        """
        Afilt_Aref = self.Alambda_Aref(self.ref_wavelength)
        Alambda_Aref = self.Alambda_Aref(eff_wavelength)
        return Alambda_Aref / Afilt_Aref

    def extinction_at_lambda(
            self, eff_wavelength: float,
            A_or_E: float or np.ndarray) -> float or np.ndarray:

        """
        Returns the extinction at a given effective wavelength
        Arguments
        ---------
        eff_wavelength : float
            Effective Wavelength of the filter for which the extinction should be determined.
            in micrometer
        A_or_E : float, numpy array
            total extinction or color excess provided by the Extinction Map
            the type is stored as A_or_E_type
        """
        Alambda_Aref = self.Alambda_Amap(eff_wavelength)
        if self.A_or_E_type.startswith("A"):
            # a total extinction was provided.
            A_ref = A_or_E
            return Alambda_Aref * A_ref

        if self.A_or_E_type.startswith("E"):
            # a color excess is provided
            E_ref_ref2 = A_or_E
            Aref2_Aref = self.Alambda_Amap(self.ref_wavelength2)
            A_ref = E_ref_ref2 / (1 - Aref2_Aref)
            return Alambda_Aref * A_ref

    def set_map_properties(self, map_name, ref_wavelength, A_or_E_type):
        # set the information of the output of the extinction map
        self.map_name = map_name
        self.ref_wavelength = ref_wavelength
        self.A_or_E_type = A_or_E_type
        if A_or_E_type.startswith("E"):
            self.ref_wavelength2 = ref_wavelength2

class ExtinctionMap(ABC):
    """
    Parent class to represent an Extinction Map:

    Attributes
    ----------
    extinction_map_name : str
        name of the Extinction Map
    l_deg : float
        galactic longitude in degree set by "update_line_of_sight"
    b_deg : float
        galactic latitude in degree set by "update_line_of_sight"

    ref_wavelength : float 
        reference wavelength for the extinction
        for an A_lambda map, lambda
        for an E(lambda1-lambda2) map, lambda1
        
    ref_wavelength2 : float
        secondary reference wavelength for the extinction
        for an A_lambda map, none
        for an E(lambda1-lambda2) map, lambda2

    A_or_E : float or function
        total extinction or color excess, from the extinction map.
        if it is a function it will be called

    A_or_E_type : str
        Output type from the extinction map.
        If it starts with "A", A_or_E is handled  as a total extinction.
        If it starts with "E": A_or_E is handled as a color excess.

    Methods
    -------
    update_line_of_sight(l_deg, b_deg) :
        specified the new galactic coordinates

    update_extinction():
        placeholder for function that updates the total extinction or color excess
        in self.extinction_map_name

    get_map_properties():
        returns the basic parameters of the extinction map
        used for Communication between ExtinctionLaw and ExtinctionMap

    """

    def __init__(self, logger: ModuleType = None, **kwargs):
        self.logger = logger
        # Basic properties of the extinction map.
        self.extinction_map_name = "NONE"
        self.ref_wavelength = None
        self.ref_wavelength2 = None
        self.A_or_E_type = 'for example A_V or E(B-V)'

    @abstractmethod
    def extinction_in_map(self, radius: float):
        """
        Estimates the extinction for the current sight-line and radial distance
        store the result into self.extinction_in_map.

        Parameters
        ----------
        radius: float [kpc]
            radial distance of the current slice

        """
        self.extinction_in_map = None
        raise NotImplementedError('No extinction map set')

    def get_map_properties(self) -> Tuple[str, float, str]:
        """
        return the extinction map properties,
        used to communicate between extinction map, and extinction law
        in case of multiple extinction laws.
        """
        return self.extinction_map_name, self.ref_wavelength, self.A_or_E_type


def CombineExtinction(ext_map=ExtinctionMap, ext_law=ExtinctionLaw) \
        -> "Extinction":
    """
    This function combines an Extinction Map and an ExtinctionLaw to a combined Class

    Parameters
    ----------
    ext_map : ExtinctionMap_class
        Class for the Extinction map. Must not be initialized
    
    ext_law : ExtinctionLaw_class or List[ExtinctionMap_class]
        Class or list of classes for the Extinction Law
        If it is a list, it uses multiple extinction laws
        Use the which_extinction_law to specify which Extinction Law should be adopted
    """
    if inspect.isclass(ext_law):
        ext_law = {ext_law}  # convert a single class into a set
    else:
        ext_law = set(ext_law)  # removing duplicates in a list

    class Extinction(ext_map, *ext_law):
        """
        Combination of an ExtinctionMap class, and one or more ExtinctionLaw classes.

        Attributes
        ----------

        Methods
        --------

        """

        def __init__(
                self, ext_map_kwargs: ModuleType,
                ext_law_kwargs: ModuleType,
                logger: ModuleType = None
                ):
            self.logger = logger
            self.bands = []
            self.eff_wavelengths = {}
            self.ext_law_index = {}
            self.ext_law_dict = {}
            self.extinction_dict = {}
            self.ext_law_kwargs = ext_law_kwargs

            # initializing the extinction law
            if len(ext_law) == 1:  # only one extinction law is set
                if isinstance(ext_law_kwargs, list):
                    ext_law_kwargs = ext_law_kwargs[0]
                ext_law.pop().__init__(self,logger=self.logger, **ext_law_kwargs.init_kwargs)
                self.multi_laws = False

            else:
                # multiple extinction laws are set.
                # needs to initialize each extinction law separately,
                # to avoid overwriting the functions
                self.multi_laws = True
                # store initialized extinction laws in a dictionary to prevent overwriting
                for i, law in enumerate(ext_law):
                    self.ext_law_dict[i] = law(logger=self.logger, **ext_law_kwargs[i].init_kwargs)

            # initializing the extinction map
            ext_map.__init__(self, logger=self.logger, **ext_map_kwargs.init_kwargs)

            if self.multi_laws:
                # pass properties from extinction map to extinction law
                for i, law in enumerate(ext_law):
                    self.ext_law_dict[i].set_map_properties(*self.get_map_properties())

        def get_extinctions(
                self, l_deg: np.ndarray, b_deg: np.ndarray, distance_kpc: np.ndarray
                ) -> Tuple[np.ndarray, Dict]:
            """

            Parameters
            ----------
            l_deg, b_deg
            distance_kpc

            Returns
            -------
            ext_in_map : extinctions at the given location

            extinction_dict : dict
                extinction for each band
            """

            ext_in_map = self.extinction_in_map(l_deg, b_deg, distance_kpc)
            
            if self.multi_laws:  # multiple extinction laws.
                extinction_dict = {
                    band: self.ext_law_dict[self.ext_law_index[band]].extinction_at_lambda(
                        self.eff_wavelengths[band], ext_in_map
                        )
                    for band in self.bands
                    }

            else:  # only one extinction law
                extinction_dict = {
                    band: self.extinction_at_lambda(
                        self.eff_wavelengths[band], ext_in_map)
                    for band in self.bands
                    }

            return ext_in_map, extinction_dict

        def get_ext_law_index(self, eff_wavelengths: float, band: str) -> int:
            """
            This function selects the extinction law
            When multiple extinction laws are provided:

            Arguments
            ---------
            eff_wavelengths : float
                Effective Wavelength of the filter for which the extinction should be determined.
                in micrometer
            band : string
                Name of the filter.
                used when an extinction law is defined for specific wavelength bands

            """
            for i, ext_law_kwarg in enumerate(self.ext_law_kwargs):
                # check if the band is specified 
                if band in ext_law_kwarg.bands:
                    return i

                if ext_law_kwarg.lambdas is not None \
                        and eff_wavelengths not in ext_law_kwarg.lambdas:
                    continue
                if ext_law_kwarg.min_lambda is not None \
                        and eff_wavelengths < ext_law_kwarg.min_lambda:
                    continue
                if ext_law_kwarg.max_lambda is not None \
                        and eff_wavelengths < ext_law_kwarg.max_lambda:
                    continue

                return i

        def set_bands(self, bands: List[str], eff_wavelengths: Dict[str, float]):
            """
            set the wavelength bands and effective wavelengths

            Parameters
            ----------
            bands: List[str]
                List of wavelength bands, for which the extinctions will be estimates
            eff_wavelengths: dict[str, float]
                Dictionary with all the effective wavelengths for each band.
            """
            self.bands = bands
            self.eff_wavelengths = eff_wavelengths

            if self.multi_laws:  # find the corresponding extinction law for each band.
                self.ext_law_index = {band: self.get_ext_law_index(
                    eff_wavelengths[band], band)
                    for band in bands}

    return Extinction
