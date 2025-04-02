"""
This is a module containing the base classes for the evolution process.

The Evolution class is split into a EvolutionIsochrones and EvolutionInterpolator.
This allows to vary both independently.
Those are combined into a common Evolution Class using the CombineEvolution function

A EvolutionIsochrones subclass should also have an EvolutionInterpolator as second Parent
to specify a default interpolator. This interpolator is used whenever no Interpolator is defined.
"""

__all__ = [
    "EvolutionIsochrones", "EvolutionInterpolator", "CombineEvolution",
    "ISOCHRONES_DIR", "EVOLUTION_DIR", "MUST_HAVE_COLUMNS"]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-06-29"

import os
from abc import ABC, abstractmethod
from types import ModuleType

import numpy as np

from .. import const

ISOCHRONES_DIR = const.ISOCHRONES_DIR
EVOLUTION_DIR = os.path.dirname(__file__)

# columns needed by the interpolator
MUST_HAVE_COLUMNS = ["log10_isochrone_age_yr", "initial_mass", "[Fe/H]_init", "phase"]


class EvolutionIsochrones(ABC):
    """
    The EvolutionIsochrones base class. The appropriate subclass is
    assigned based on Population[evolution_class][name] through the "get_subclass" factory.
    must specify at least this "log10_isochrone_age_yr", "initial_mass", "[Fe/H]_init", "star_mass"

    Current included subclasses:
    ----------------------------
    -MIST': Using the MIST isochrones systems,


    Attributes
    ----------
    isochrones :  Pandas DataFrame
        DataFrame of the loaded isochrones
    isochrones_grouped : Pandas DataFrameGroupBy
        Isochrones groupend in metallicity and age

    more Attributes can be specified in the subclasses

    Methods
    -------
    __init__(self,**kwargs) : None
        initialize the EvolutionIsochrones class from a Population dictionary
    get_isochrones(self) : None
        return the Isochrones Data Frame
    """
    isochrones_name = None
    isochrones_grouped = None
    def __init__(self, columns=None, logger: ModuleType = None,  **kwargs):
        """

        Parameters
        ----------
        columns: list of columns to be loaded
        kwargs : additional parameters specified in the evolution_class
        """
        super().__init__(**kwargs)
        self.columns = columns
        self.logger = logger

    def get_isochrones(self):
        """
        Returns
        -------
        isochrones : Pandas.DataFrame
            Isochrone Grid
        """
        return self.isochrones

    def get_mass_min(
            self, band: str, magnitude_limit: float, distances: np.ndarray or float,
            max_age: float = None, extinction_at_slice_front: np.ndarray or float = None
            ) -> np.ndarray:
        """
        gets the minimum mass of stars still passing the magnitude_limit

        Parameters
        ----------
        magnitude_limit : float
            maximum observed magnitude, (set distance to 0.01 if absolute magnitudes),
        distances : np.ndarray or float
            distances in kpc
        max_age: float [Gyr]
        Returns
        -------
        min_mass :  np.ndarray or float
        Returns minimum mass for stars brighter than magnitude_limit
        """

        abs_mag_lim = magnitude_limit - 5 * np.log10(np.maximum(distances, 1e-8) * 100)
        #print('abs_mag_lim_init',abs_mag_lim)
        if extinction_at_slice_front is not None:
            abs_mag_lim -= extinction_at_slice_front

        #print('abs_mag_lim',abs_mag_lim)
        #print('ages', max_age, self.iso_ages)

        if max_age is None:
            max_age = np.log10(self.iso_ages).max()
        else:
            next_age = self.iso_ages[np.searchsorted(self.iso_ages[:-1], max_age * 1e9)]
            max_age = np.log10(next_age)
        oldest = self.isochrones_grouped.get_group((min(self.file_met),max_age))

        no_mass_loss = oldest[1 - oldest['star_mass'] / oldest['initial_mass'] < 1e-4]
        yes_mass_loss = oldest[1 - oldest['star_mass'] / oldest['initial_mass'] >= 1e-4]
        #print(yes_mass_loss['initial_mass'],yes_mass_loss['Bessell_I'])
        #print(self.isochrones_grouped['initial_mass'])
        closest = np.searchsorted(-no_mass_loss[band].iloc[1:], -abs_mag_lim)
        masses = no_mass_loss.iloc[closest].initial_mass.values
        masses[closest < 10] = 0
        #print(masses)
        closest = np.searchsorted(-no_mass_loss[band].iloc[1:], -abs_mag_lim)
        masses = no_mass_loss.iloc[closest].initial_mass.values
        masses[closest < 10] = 0
        return masses


class EvolutionInterpolator(ABC):
    """
    The EvolutionInterpolator base class. The appropriate subclass is
    assigned based on Population[evolution_class][Interpolator] through the "get_subclass" factory.
    Each EvolutionIsochrone Class should have a specified default Interpolator

    Current included subclasses:
    ----------------------------
    -LagrangeInterpolator: Using Lagrange Polynomials to interpolate
            between Consecutive stars in the Isochrone grid


    Attributes
    ----------
    accept_np_arrays :  Bool
        should be set to True if the interpolator can work with multiple stars given as array.

    more Attributes can be specified in the subclasses

    Methods
    -------
    __init__(self,**kwargs) : None
        initialize the EvolutionIsochrones class from a Population dictionary
    get_evolved_props(self) : pandas.DataFrame
        return the Isochrones

    """
    interpolator_name = None
    # specify if the Interpolator can hande numpy arrays of multiple stars
    accept_np_arrays = False

    # AttributeError which needed to be provided by EvolutionIsochrones
    isochrones = None
    met_to_file_iso = None
    file_met = None
    iso_ages = None
    def __init__(self,  logger: ModuleType = None, **kwargs):
        self.logger = logger


    @abstractmethod
    def get_evolved_props(self, m_init, met, age, props, **kwargs):
        """
        Get the interpolated properties from the isochrones.

        Parameters
        ----------
        m_init : float, ndarray [Msun]
            initial mass for the stars
        met : float, ndarray [dex]
            initial metallicity for the stars
        age : float, ndarray [Gyr]
            initial age for the stars
        props : set
            a set of properties
        kwargs : dict, optional
            further ignored keywords


        Returns
        -------
        s_track : dict
            interpolator properties for all keywords specified in props
        in_grid : bool , ndarray
            indicator if the star is within the Isochrone Grid
        final_phase : bool, ndarray
            indicator if the star is in its final phase and might need special treatment.
            (can be the same as in_grid)

        Note: if accept_np_arrays is set to True.
            each entry in s_track has to be a ndarray of the same shape as the input parameters

        """
        raise NotImplementedError('No interpolator subclass defined')
        return s_track, in_grid, final_phase


def CombineEvolution(Isochrones=None, Interpolator=None):
    """
    This function combines the isochrones and Interpolator into a common Evolution Class
    Parameters
    ----------
    Isochrones : class
        Class of the Isochrone system
    Interpolator : class or None
        Class of the Interpolator
    Returns
    -------
    Evolution: class
    """
    if Isochrones is None:
        raise NotImplementedError('No Isochrones System is specified ')

    if (Interpolator is None) or (Isochrones in Interpolator.__subclasses__()):
        # The provided interpolator is the standard interpolator for the Isochrones class
        class Evolution(Isochrones):
            """
            A Combined Class of an Isochrones grid and an interpolator.
            """

            def __init__(self, int_kwargs=None, iso_kwargs=None,  logger: ModuleType = None):
                """
                Parameters
                ----------
                int_kwargs: dict, None optional
                    keyword arguments for the interpolator (Not used)
                iso_kwargs: dict
                     keyword arguments for the Isochrone
                """
                # initialize Isochrones
                if iso_kwargs is None:
                    iso_kwargs = {}
                super().__init__(logger=logger, **iso_kwargs)

                # add Isochrones docstring
                if Isochrones.__doc__ is not None:
                    self.__doc__ = ''.join(
                        (self.__doc__, 'Docstring from Isochrones:', Isochrones.__doc__)
                        )

        return Evolution

    else:
        class Evolution(Interpolator, Isochrones):
            """
            A Combined Class of an Isochrones grid and an Interpolator.
            """

            def __init__(self, int_kwargs=None, iso_kwargs=None, logger: ModuleType = None):
                """
                Parameters
                ----------
                int_kwargs: dict, None optional
                    keyword arguments for the interpolator
                iso_kwargs: dict
                     keyword arguments for the Isochrone
                """


                # initialize Isochrones
                if iso_kwargs is None:
                    iso_kwargs = {}
                Isochrones.__init__(self, **iso_kwargs)

                # initialize Interpolator
                if int_kwargs is None:
                    int_kwargs = {}
                Interpolator.__init__(self, **int_kwargs)
                self.logger = logger
                # add Isochrones docstring
                if Isochrones.__doc__ is not None:
                    self.__doc__ = ''.join(
                        (self.__doc__, 'Docstring from Isochrones:', Isochrones.__doc__)
                        )
                # add Interpolator docstring
                if Interpolator.__doc__ is not None:
                    self.__doc__ = ''.join(
                        (self.__doc__, 'Docstring from Interpolator:', Interpolator.__doc__)
                        )

        return Evolution
