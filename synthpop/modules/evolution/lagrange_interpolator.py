"""
An interpolator for the isochrone grid

Uses a 2nd order Lagrange polynomial over mass (linear at grid edges)
and linear over age and metallicity
"""
__all__ = ["LagrangeInterpolator"]

import numpy as np
from ._evolution import EvolutionInterpolator
import matplotlib.pyplot as plt


class LagrangeInterpolator(EvolutionInterpolator):
    """
    Performs an interpolation using lagrangian polynomials
    2 order for the mass (linear if it is at the edge of the isochrone grid)
    and a linear interpolation for age and  metallicity

    Methods
    -------
    lagrange_poly(x, grid, fgrid) : ndarray
        performs interpolation using lagrange polynomials
    grid_n4(grid, value) : ndarray, ndarray
        get the closest points from the isochrones mass grid
    combined_grid(met_key, mass, iso_ages, props) : ndarray, ndarray
        get the grid points for different metallicities and ages
    interp_props(self,track,mass,age,props) : dict
        interp the properties using the isochrones.
    get_evolved_props () : dict, ndarray
    """
    interpolator_name = 'Lagrange_Polynomials'

    accept_np_arrays = True

    def __init__(self, isochrones=None, **kwargs):
        super().__init__(**kwargs)
        if isochrones is not None:
            self.isochrones = isochrones
            self.isochrones_grouped = isochrones.groupby(["[Fe/H]_init", "log10_isochrone_age_yr"])
        self.log_ages = self.isochrones['log10_isochrone_age_yr'].unique()
        self.lin_ages = 10 ** (self.log_ages - 9)
        self.iso_met = self.isochrones['[Fe/H]_init'].unique()

    # @staticmethod
    def lagrange_poly(self, x, grid, fgrid):
        """
        cubic interpolation function using the lagrange polynomial
        Poly = \\sum_{i} f_{i} * \\Pi_{j, j!=i} (x-xj)/(xi-xj)
        Interpolation function
        Uses a lower order if the grid has fewer than 4 unique grid points
        (combined grid is defined in a way that it produces 4,2 or 1 unique grid points)

        Parameters
        ----------
        x: ndarray
            values at which the grid should be interpolated
        grid: ndarray
            grid points. must have shape [..., 1, j]
        fgrid: ndarray
            function values for each grid point
            must have shape [..., i, j] where

        Returns
        -------
        f: ndarray
            interpolated value for f(x)
        """

        x = x.reshape(grid.shape[:-1])
        f = np.zeros(fgrid.shape[:-1])

        for i in range(4):  # \\sum_{i} f_{i} ...
            f_i = fgrid[..., i].copy()
            for j in range(4):  # \\Pi_{j, j!=i} (x-xj)/(xi-xj)
                if i == j:
                    continue
                # distance between grid point i and j
                dij = grid[..., i] - grid[..., j]
                dn = dij[..., 0] != 0
                # check if the grid j appears twice:
                for k in range(j):
                    dn *= (grid[..., 0, k] != grid[..., 0, j])
                # multiply with (x-xj)/(xi-xj)
                f_i[dn] *= (x[dn] - grid[dn, :, j]) / dij[dn]
                # check if grid point exist twice, and only select first one
                # i.e. force polynomial to lower order
                if j < i:
                    f_i[dij[..., 0] == 0] *= 0  # else 1
            f += f_i
        return f

    @staticmethod
    def grid_n4(grid, value):
        """
        Gets the nearest four points in grid to value.
        If the value is inside the grid it will return the
        two values on either side (y) and their indices (x)
        If the value is outside the grid it returns 4 times the closest value ant its indices
        If the value is below/above the second lowest/or highest value it returns 2 times the
        next-lowest/second-highest value and two time the second-lowest/next-highest value
        and their indices
        This triggers the interpolation to switch to a linear interpolation between both values

        Parameters
        ----------
        grid: pandas.DataFrame
        value: ndarray

        Returns
        -------
        indices : ndarray
            indices for the 4 closest neighbours (two upper + two lower)
        values  : ndarray
            grid values for the 4 closest neighbours (two upper + two lower)
        """

        if not isinstance(value, np.ndarray):
            value = np.array(value)
            if len(value.shape) == 0:
                value = value.reshape(-1)
        edge = len(grid)

        # make sure the grid is sorted
        if not np.array_equal(grid, np.sort(grid)):
            print(value, grid)
            raise ValueError("Grid is not sorted")

        # arrays to store the results
        ind = np.zeros((*value.shape, 4), int)
        val = np.zeros((*value.shape, 4))

        # back tracking revels that value is never a list, so the output will not be a list
        # find the index of where the value would fit on the grid
        index = np.searchsorted(grid, value)
        # ind == len(grid) if value is outside the grid

        # outside grid
        if min(index) == 0:
            # first grid point
            ind_outside_lower = [0, 0, 0, 0]
            ind[index == 0] = ind_outside_lower
            val[index == 0] = grid.iloc[ind_outside_lower]
        if max(index) == edge:
            # last grid point
            ind_outside_upper = [edge - 1, edge - 1, edge - 1, edge - 1]
            ind[index == edge] = ind_outside_upper
            val[index == edge] = grid.iloc[ind_outside_upper]

        good = (index != 0) & (index != edge)
        # max one value below:

        only_one_below = (grid.iloc[index - 1] == grid.iloc[0]).to_numpy()
        # only one grid point above or below
        # [i,i, j,j], will cause a linear interpolation between grid[i] and grid[j]

        if any(only_one_below):
            ind_only_one_below = [0, 0, 1, 1]
            ind[only_one_below] = ind_only_one_below
            val[only_one_below] = grid.iloc[ind_only_one_below]

        only_one_above = np.zeros(value.shape, bool)
        only_one_above[good] = (grid.iloc[index[good]] == grid.iloc[-1]).to_numpy()
        if any(only_one_above):
            # index is not necessary edge-2 due to multiple rows with same mass
            # it can only be true for one specific i, hence index[only_one_above] are all identical
            # use index[only_one_above][0] to select the value
            # np.searched returns the first index with grid>value, so index[only_one_below] is
            ind_only_one_above = [index[only_one_above][0] - 1, index[only_one_above][0] - 1,
                index[only_one_above][0], index[only_one_above][0]]
            ind[only_one_above] = ind_only_one_above
            val[only_one_above] = grid.iloc[ind_only_one_above]

        # select all the cases with two or more grid points on each side
        good = (index != 0) \
               & (only_one_below == False) \
               & (index != edge) \
               & (only_one_above == False)

        index2 = index[good]
        value2 = grid.iloc[index2]
        # next lower value
        value1 = grid.iloc[index2 - 1]
        # index for next lower value
        index1 = np.searchsorted(grid, value1)
        # second lower value
        value0 = grid.iloc[index1 - 1]
        # second lower value
        index0 = np.searchsorted(grid, value0)
        # next higher index, increase value by a small epsilon
        index3 = np.searchsorted(grid, value2 + 1e-10)
        # next higher value
        value3 = grid.iloc[index3]

        # add indices and values to ind and val
        ind[good] = np.array([index0, index1, index2, index3]).T
        val[good] = np.array([value0, value1, value2, value3]).T

        return ind, val

    def combined_grid(self, met_key, mass, iso_ages, props):
        """
        Creates a combined grid for the different stars,
        Uses the isochrone based on each individual metallicity and age
        uses the grid_n4 to find adjacent entries for the masses

        creates a numpy data array which contains the properties
        for the selected grid cells

        Parameters
        ----------
        met_key :ndarray
        mass : ndarray
        iso_ages : ndarray
        props : list

        Returns
        -------
        grid : ndarray
        data : ndarray
        """
        max_mass = np.zeros(len(iso_ages))
        data = np.zeros((len(iso_ages), len(props), 4))
        grid = np.zeros((len(iso_ages), 1, 4))

        # create grid
        for current_met in np.unique(met_key):
            #track_dif_ages = self.isochrones[current_met]
            for current_age in np.unique(iso_ages):
                common_met_and_age = (
                        (iso_ages == current_age)
                        & (met_key == current_met))
                if not common_met_and_age.any():
                    continue
                # get all masses with for the closest age
                #closest_grid_points = track_dif_ages['log10_isochrone_age_yr'] == current_age
                #masses = track_dif_ages[closest_grid_points]["initial_mass"]
                current_iso = self.isochrones_grouped.get_group((current_met, current_age))
                grid_index, grid_masses = self.grid_n4(current_iso.initial_mass, mass[common_met_and_age])
                grid[common_met_and_age, 0, :] = grid_masses
                max_mass[common_met_and_age] = current_iso.initial_mass.max()
                # evaluate the properties at the selected grid points
                for i, item in enumerate(props):
                    #data[common_met_and_age, i] = (
                    #    track_dif_ages[closest_grid_points].iloc[grid_index.ravel()][item]
                    #).values.reshape(grid_index.shape)
                    data[common_met_and_age, i] = (
                        current_iso.iloc[grid_index.ravel()][item]
                        ).values.reshape(grid_index.shape)
        return grid, data

        # Interpolate each property to the proper metallicity

    # for each mass in the mass grid for each eep.

    def interp_props(
            self, met_key, log_iso_ages, mass, props,
            met_key_higher=None, met=None, log_iso_ages_higher=None, log_ages=None,
            age_interpolation='linear', **kwargs
            ):
        """
        Interpolates the isochrones and returns the value at the given
        mass metallicity and age. 
        Determines the closets 2 upper and lower neighbours grid and 
        performs a cubic interpolation using the lagrange polynomials.
        At the edge of isochrone grid it only performs a linear
        interpolation and outside it uses the closest grid point.
        If met_key_higher and/or log_iso_ages_higher it also performs a
        linear interpolation in metallicity and or age 
        
        
        Parameters:
        -----------
        met_sting: ndarray
            numpy array of metallicity dictionary key for each star.
            If met_key_higher is given it is the lower metallicity for the interpolation,
            otherwise the closets metallicity and no interpolation is performed
        log_iso_ages: ndarray
            numpy array of iso_ages as defined in the isochrones,
            If log_iso_ages_higher is given it is the lower log_iso_ages for the interpolation,
            otherwise the closets log_iso_ages and no interpolation is performed

        mass : ndarray
            initial mass for each star

        props : Set 
            list of properties which should be derived from the interpolation
        
        met_key_higher :  ndarray, None, optional
            numpy array of metallicity dictionary keys for each star.
            used us higher value for the interpolation.
            If None, an interpolation is not performed
        met :  ndarray, None, optional
            initial metallicity for each star. used when interpolating the stars
        log_iso_ages_higher : ndarray, None, optional
            numpy array of metallicity dictionary key for each star.
            used us higher value for the interpolation.
            If None, an interpolation is not performed
        log_ages : ndarray, None, optional
            log10(age) for each star. used when interpolating the stars

        age_interpolation : string, optional
            type of age interpolation.
            if 'linear', it performs interpolation in age,
            if 'log', it performs an interpolation in log(age) space
        
        Returns
        -------
        inter_props
            returns a dictionary with interpolated values  for each property that looks like this:
            {'star_mass': ndarray, 'log_Teff': ndarray, 'log_g': ndarray, 'log_L': ndarray,
            'initial_mass': ndarray, '2MASS_Ks': ndarray, '2MASS_J': ndarray, '2MASS_H': ndarray,
            'Bessell_V': ndarray, 'Bessell_I': ndarray}
        """

        if met_key_higher is not None:
            # repeat all stars with different met_key
            met_key = np.hstack([met_key, met_key_higher])
            log_iso_ages = np.hstack([log_iso_ages, log_iso_ages])
            if log_iso_ages_higher is not None:
                log_iso_ages_higher = np.hstack([log_iso_ages_higher, log_iso_ages_higher])
                log_ages = np.hstack([log_ages, log_ages])
            mass = np.hstack([mass, mass])

        if log_iso_ages_higher is not None:
            # repeat all stars with different log_iso_ages
            log_iso_ages = np.hstack([log_iso_ages, log_iso_ages_higher])
            mass = np.hstack([mass, mass])
            met_key = np.hstack([met_key, met_key])
        props.add('phase')
        # get grid data
        grid, data = self.combined_grid(met_key, mass, log_iso_ages, props)
        in_grid = grid[:, 0, 1] != grid[:, 0, 2]
        phase_index = [i for i, j in enumerate(props) if j == 'phase'][0]

        # call interpolation
        result = self.lagrange_poly(mass, grid, data)
        # result[:,phase_index] = np.floor(result[:,phase_index])
        if log_iso_ages_higher is not None:
            # interpolate in age
            # split table in half:
            result_lower, result_higher = np.array_split(result, 2)
            in_grid = np.logical_or(*np.array_split(in_grid, 2))
            log_iso_ages, log_iso_ages_higher = np.array_split(log_iso_ages, 2)

            # Use lagrange polynomial if inside grid
            if age_interpolation == 'log':
                age_lower = log_iso_ages
                age_higher = log_iso_ages_higher
                age = log_ages
            else:
                age_lower = 10 ** log_iso_ages
                age_higher = 10 ** log_iso_ages_higher
                age = 10 ** log_ages
            weight1, weight2 = self.get_weights_from_charon(
                result_lower[:, phase_index], result_higher[:, phase_index], age_lower, age_higher,
                age)
            result = weight1[:, np.newaxis] * result_lower + weight2[:, np.newaxis] * result_higher
            met_key = np.array_split(met_key, 2)[0]

        if met_key_higher is not None:
            # interpolate in met
            result_lower, result_higher = np.array_split(result, 2)
            in_grid = np.logical_or(*np.array_split(in_grid, 2))
            met_key, met_key_higher = np.array_split(met_key, 2)

            weight1, weight2 = self.get_weights_from_charon(
                result_lower[:, phase_index], result_higher[:, phase_index], met_key,
                met_key_higher, met,)
            result = weight1[:, np.newaxis] * result_lower + weight2[:, np.newaxis] * result_higher
        # transform results into a dictionary
        output = {item: result[..., i] for i, item in enumerate(props)}
        return output, in_grid

    def get_weights_from_charon(self, phase1, phase2, x1, x2, x, trans=5.5):
        """ charon transfers the "living" stars at the
        red giant clump to the "death" white Dwarfs"""

        weight1 = np.ones(x.shape) * 0.5
        not_zero = (x1 - x2) != 0
        weight1[not_zero] = (x2 - x)[not_zero] / (x2 - x1)[not_zero]  # is 1 if x=x_higher
        phase_av = weight1 * phase1 + (1 - weight1) * phase2

        # stars in phase transition
        dying = (np.maximum(phase2, phase1) > 5) & (phase1 != phase2)

        # stars already death
        death = phase_av[dying] > trans
        # check which phase contains the larger value
        phase2_is_higher = phase2[dying] >= phase1[dying]
        # combine death & phase2_is_higher
        # phase2_is_higher = True , death = True => weight1 = 0 weight2 = 1
        # phase2_is_higher = False, death = True => weight1 = 1 weight2 = 0
        # phase2_is_higher = True , death = False => weight1 = 1 weight2 = 0
        # phase2_is_higher = False , death = False => weight1 = 0 weight2 = 1
        weight1[dying] = np.logical_xor(phase2_is_higher, death)

        return weight1, 1 - weight1

    def get_evolved_props(
            self, m_init, met, age, props, age_interpolation='linear',
            met_interpolation=True, *kwargs
            ):
        """
        Parameters
        ----------
        m_init : ndarray [Msun]
        met : ndarray [dex]
        age : ndarray [Gyr]
        props : Set
        age_interpolation : string, optional
        met_interpolation : bool, optional
        kwargs : dict, optional

        Returns
        -------
        s_track : dict
        in_grid : ndarray
        """
        # TODO Get a more flexible way to get the closes  the metallicities and ages
        if met_interpolation:
            # get higher and lower grid point for metallicities
            met_index = np.searchsorted(self.file_met, met)
            # uses the closest edge point, if metallicity is outside self.file_met
            met_key = self.file_met[np.maximum(met_index - 1, 0)]
            met_key_higher = self.file_met[np.minimum(met_index, len(self.file_met) - 1)]
        else:
            # estimates the midpoint between two file metallicities
            met_key_higher = None
            file_met_midpoints = (self.file_met[:-1] + self.file_met[1:]) / 2
            # estimates the index for the closets metallicity
            met_index = np.searchsorted(file_met_midpoints, met)
            # estimates ISO string
            met_key = self.file_met[met_index]

        log_ages = np.log10(age) + 9
        if age_interpolation:
            # estimates higher and lower grid point
            log_iso_ages_index = np.searchsorted(self.log_ages, log_ages)
            log_iso_ages = self.log_ages[np.maximum(log_iso_ages_index - 1, 0)]
            log_iso_ages_higher = self.log_ages[
                np.minimum(log_iso_ages_index, len(self.log_ages) - 1)]
        else:
            log_iso_ages_higher = None

            # estimates the closest age in the linear grid
            age_midpoints = (self.iso_ages[:-1] + self.iso_ages[1:]) / 2
            index_age = np.searchsorted(age_midpoints, age)
            closest_age = self.iso_ages[index_age]
            # translate to log10_isochrone_age_yr, ensures the correct rounding
            unique_log_ages_midpoints = (self.log_ages[:-1] + self.log_ages[1:]) / 2
            log_iso_ages = self.log_ages[np.searchsorted(
                unique_log_ages_midpoints, np.log10(closest_age))]

        # Interpolate properties
        s_track, in_grid = self.interp_props(met_key, log_iso_ages, m_init, props,
            met_key_higher=met_key_higher, met=met,
            log_iso_ages_higher=log_iso_ages_higher, log_ages=log_ages,
            age_interpolation=age_interpolation)
        return s_track, in_grid, np.logical_not(in_grid)
