"""
An interpolator for the isochrone grid

Within each isochrone, it uses a cubic interpolation between adjacent grid-points.

It aligns the phases between adjacent ages
and performs a linear interpolation in metallicity and age.

See Klüter & Huston et al (SynthPop Paper I, in prep).
"""

__all__ = ["CharonInterpolator", ]
__author__ = "J. Klüter"
__date__ = "2023-03-01"

from typing import Set
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, RectBivariateSpline
import matplotlib.pyplot as plt
from ._evolution import EvolutionInterpolator

class CharonInterpolator(EvolutionInterpolator):
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
    interpolator_name = 'Charon_Polynomials'

    accept_np_arrays = True

    def __init__(self, isochrones=None, props_no_charon=None, coins=2, **kwargs):
        super().__init__(**kwargs)
        if coins < 2:
            raise ValueError("You needed 2 coins to pay Charon to cross the river Styx.")
        if isochrones is not None:
            self.isochrones = isochrones
            self.isochrones_grouped = isochrones.groupby(["[Fe/H]_init", "log10_isochrone_age_yr"])

        (self.tips, self.tips0_func, self.tips1_func,
        self.tips2_func, self.endp_func) = self.get_tip_func()
        self.log_ages = self.isochrones['log10_isochrone_age_yr'].unique()
        self.lin_ages = 10 ** (self.log_ages - 9)
        self.iso_met = self.isochrones['[Fe/H]_init'].unique()
        self.props_no_charon = {"initial_mass", "[Fe/H]_init", }

        if props_no_charon is not None:
            self.props_no_charon.update(set(props_no_charon))

    def get_tip_func(self):
        # collect phase changes in the isochrone for each met and age
        iso_p0 = self.isochrones[self.isochrones['phase'] == 0]
        iso_p3 = self.isochrones[self.isochrones['phase'] == 3]
        iso_p4 = self.isochrones[self.isochrones['phase'] == 4]
        iso_p6 = self.isochrones[self.isochrones['phase'] == 6]

        grouped0 = iso_p0.groupby(['[Fe/H]_init', 'log10_isochrone_age_yr', ])
        grouped3 = iso_p3.groupby(['[Fe/H]_init', 'log10_isochrone_age_yr', ])
        grouped4 = iso_p4.groupby(['[Fe/H]_init', 'log10_isochrone_age_yr', ])
        grouped6 = iso_p6.groupby(['[Fe/H]_init', 'log10_isochrone_age_yr', ])

        tip0 = grouped0.initial_mass.max()
        tip1 = grouped3.initial_mass.min()
        tip2 = grouped6.initial_mass.min()
        tip2_alt2 = grouped3.initial_mass.max()
        tip2_alt1 = grouped4.initial_mass.max()

        end_point = grouped6.initial_mass.max()
        tip0.name = 'tip0'
        tip1.name = 'tip1'
        tip2.name = 'tip2'
        tip2_alt2.name = 'tip2_alt2'
        tip2_alt1.name = 'tip2_alt1'
        end_point.name = "endp"
        tips = pd.concat([tip1, tip2, tip2_alt1, tip2_alt2, end_point], axis=1)
        tips['width'] = tips.tip2 - tips.tip1

        # fill with alternative estimates for tip2 location
        where = np.isnan(tips.tip2) & (np.logical_not(np.isnan(tips.tip2_alt1)))
        tips.tip2[where] = tips.tip2_alt1[where]

        tips['width'] = tips.tip2 - tips.tip1
        tt = tips.groupby('[Fe/H]_init')
        idx = tt.tip2.idxmin()
        widths = tips.width[idx]
        where = np.isnan(tips.tip2) & (tips.tip1 < 1.1)
        #
        tips.width[where] = widths[idx[tips[where].index.get_level_values(0)]]
        tips.tip2[where] = tips.tip1[where] + tips.width[where]

        where = np.isnan(tips.tip2) & (np.logical_not(np.isnan(tips.tip2_alt2)))
        tips.tip2[where] = tips.tip2_alt2[where]

        # convert table to a 2d regular grid
        mets = tips.index.levels[0]
        log_ages = tips.index.levels[1]
        # ignore stars younger than ~3 million yrs
        log_ages = log_ages[log_ages > 6.5]

        tip0_grid = np.zeros([mets.size, log_ages.size])
        tip1_grid = np.zeros([mets.size, log_ages.size])
        tip2_grid = np.zeros([mets.size, log_ages.size])
        endp_grid = np.zeros([mets.size, log_ages.size])

        for i, met in enumerate(mets):
            sub_table = tips.loc[met]
            if len(sub_table) >= len(log_ages) / 2:
                # fill missing tips
                t1 = interp1d(sub_table.index.values, sub_table.tip2,
                    fill_value='extrapolate')
                t2 = interp1d(sub_table.index.values, sub_table.tip2,
                    fill_value='extrapolate')
                end_not_nan = np.logical_not(np.isnan(sub_table.endp))
                en = interp1d(sub_table.index.values[end_not_nan], sub_table.endp[end_not_nan],
                    fill_value='extrapolate')
                for j, log_age in enumerate(log_ages):
                    tip0_grid[i, j] = tip0.get((met, log_age))
                    tip1_grid[i, j] = tips.tip1.get((met, log_age), t1(log_age))
                    tip2_grid[i, j] = tips.tip2.get((met, log_age), t2(log_age))
                    endp_grid[i, j] = tips.endp.get((met, log_age), en(log_age))
                    if np.isnan(endp_grid[i, j]):
                        endp_grid[i, j] = en(log_age)
        # 2d interpolation of tip locations
        tips0_func = RectBivariateSpline(mets, log_ages, tip0_grid, kx=1, ky=1)
        tips1_func = RectBivariateSpline(mets, log_ages, tip1_grid, kx=1, ky=1)
        tips2_func = RectBivariateSpline(mets, log_ages, tip2_grid, kx=1, ky=1)
        endp_func = RectBivariateSpline(mets, log_ages, endp_grid, kx=1, ky=1)
        tips = pd.concat([tip0, tips], axis=1)
        return tips, tips0_func, tips1_func, tips2_func, endp_func

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
               & np.logical_not(only_one_below) \
               & (index != edge) \
               & np.logical_not(only_one_above)

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

    @staticmethod
    def lagrange_poly(x, grid, fgrid):
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

    def combined_grid(self, mass, iso_met, iso_ages, props):
        """
        Creates a combined grid for the different stars,
        Uses the isochrone based on each individual metallicity and age
        uses the grid_n4 to find adjacent entries for the masses

        creates a numpy data array which contains the properties
        for the selected grid cells

        Parameters
        ----------
        mass : ndarray
        iso_met :ndarray
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
        for current_met in np.unique(iso_met):
            for current_age in np.unique(iso_ages):
                common_met_and_age = (
                        (iso_ages == current_age)
                        & (iso_met == current_met))
                if not common_met_and_age.any():
                    continue
                # get all masses with for the iso_age and iso_met
                current_iso = self.isochrones_grouped.get_group((current_met, current_age))
                grid_index, grid_masses = self.grid_n4(current_iso.initial_mass,
                    mass[common_met_and_age])
                grid[common_met_and_age, 0, :] = grid_masses
                max_mass[common_met_and_age] = current_iso.initial_mass.max()
                # evaluate the properties at the selected grid points
                for i, item in enumerate(props):
                    data[common_met_and_age, i] = (
                        current_iso.iloc[grid_index.ravel()][item]
                        ).values.reshape(grid_index.shape)

        return grid, data, max_mass

        # Interpolate each property to the proper metallicity

    # for each mass in the mass grid for each eep.

    def interp_props(
            self, mass, met_grid, log_age_grid, props
            ):
        """
        Interpolates the isochrones and returns the value at the given
        mass metallicity and age.
        Determines the closets 2 upper and lower neighbours grid and
        performs a cubic interpolation using the lagrange polynomials.
        At the edge of isochrone grid it only performs a linear
        interpolation and outside it uses the closest grid point.

        Parameters:
        -----------
        mass : ndarray
            initial mass for each star

        met_grid: ndarray
            numpy array of metallicity dictionary key for each star.

        log_age_grid: ndarray
            numpy array of iso_ages as defined in the isochrones,

        props : Set
            list of properties which should be derived from the interpolation

        Returns
        -------
        inter_props returns a dictionary with interpolated values  for each property
        """

        grid, data, max_mass = self.combined_grid(mass, met_grid, log_age_grid, props)

        in_grid = grid[:, 0, 1] != grid[:, 0, 2]

        # call interpolation in mass
        result = self.lagrange_poly(mass, grid, data)

        # transform results into a dictionary
        return result, in_grid, max_mass

    def get_modified_mass(self, mass, met, log_age, grid_met, grid_log_age):
        """
        This function estimates the modified mass.
        by shifting the mass such that the phase transition points.

        Parameters
        ----------
        mass
        met
        log_age
        grid_met
        grid_log_age

        Returns
        -------

        """

        tip0_req = self.tips0_func(met, log_age, grid=False)
        tip1_req = self.tips1_func(met, log_age, grid=False)
        tip2_req = self.tips2_func(met, log_age, grid=False)
        end_req = self.endp_func(met, log_age, grid=False)

        tip0_grid = self.tips0_func(grid_met, grid_log_age, grid=False)
        tip1_grid = self.tips1_func(grid_met, grid_log_age, grid=False)
        tip2_grid = self.tips2_func(grid_met, grid_log_age, grid=False)
        end_grid = self.endp_func(grid_met, grid_log_age, grid=False)

        mass_new = (mass - tip1_req) / (tip2_req - tip1_req) \
                    * (tip2_grid - tip1_grid) + tip1_grid
        mass_new_rgb_phase = (mass - tip0_req) / (tip1_req - tip0_req) \
                              * (tip1_grid - tip0_grid) + tip0_grid
        mass_new_wd_phase = (mass - tip2_req) / (end_req - tip2_req) \
                             * (end_grid - tip2_grid) + tip2_grid

        mass_new[mass < tip1_req] = mass_new_rgb_phase[mass < tip1_req]
        mass_new[mass > tip2_req] = mass_new_wd_phase[mass > tip2_req]
        mass_new[mass > end_req] = end_grid[mass > end_req] + 1
        # transition in the last 10%
        weight = (mass - tip0_req * 0.95) / (tip0_req * 0.05)

        weight = np.clip(weight, 0, 1)

        mass_new = (1 - weight) * mass + weight * mass_new

        return mass_new

    def get_adjacent_gridpoints(self, met, log_ages):
        """
        Estimates the adjacent metallicity and age grid-points
        Parameters
        ----------
        met : np.ndarray
            metallicity for each star
        log_ages: np.ndarray
            log10(age) for each star

        Returns
        -------
        iso_met_lower
        iso_met_higher
        log_iso_ages_lower
        log_iso_ages_higher
        iso_ages_lower
        iso_ages_higher
        """
        # get closest grid points
        met_index = np.searchsorted(self.iso_met, met)
        # uses the closest edge point, if metallicity is outside self.file_met
        # get lower and higher metallicity
        iso_met_lower = self.iso_met[np.maximum(met_index - 1, 0)]
        iso_met_higher = self.iso_met[np.minimum(met_index, len(self.iso_met) - 1)]

        log_iso_ages_index = np.searchsorted(self.log_ages, log_ages)
        log_iso_ages_lower = self.log_ages[np.maximum(log_iso_ages_index - 1, 0)]
        log_iso_ages_higher = self.log_ages[
            np.minimum(log_iso_ages_index, len(self.log_ages) - 1)]
        iso_ages_lower = self.lin_ages[np.maximum(log_iso_ages_index - 1, 0)]
        iso_ages_higher = self.lin_ages[np.minimum(log_iso_ages_index, len(self.log_ages) - 1)]

        return (iso_met_lower, iso_met_higher,
        log_iso_ages_lower, log_iso_ages_higher,
        iso_ages_lower, iso_ages_higher)

    @staticmethod
    def get_weights(met, age, iso_met_lower, iso_met_higher, iso_ages_lower, iso_ages_higher):
        """
        Estimates the weights for each grid-point

        Parameters
        ----------
        met
        age
        iso_met_lower
        iso_met_higher
        iso_ages_lower
        iso_ages_higher

        Returns
        -------

        """
        weights_met = np.zeros((len(met), 1))
        diff_grid_points = (iso_met_higher - iso_met_lower) != 0
        weights_met[diff_grid_points, 0] = ((met - iso_met_lower)[diff_grid_points]
                                            / (iso_met_higher - iso_met_lower)[diff_grid_points]
                                            )
        weights_age = np.zeros((len(met), 1))
        diff_grid_points = (iso_ages_higher - iso_ages_lower) != 0
        weights_age[diff_grid_points, 0] = ((age - iso_ages_lower)[diff_grid_points]
                                            / (iso_ages_higher - iso_ages_lower)[diff_grid_points]
                                            )

        weight_low_low = (1 - weights_age) * (1 - weights_met)
        weight_low_high = weights_age * (1 - weights_met)
        weight_high_low = (1 - weights_age) * weights_met
        weight_high_high = weights_age * weights_met

        return weight_low_low, weight_low_high, weight_high_low, weight_high_high

    def get_evolved_props(
            self, m_init: np.ndarray, met: np.ndarray,
            age: np.ndarray, props: Set, inter_age: str = "linear"
            ):
        """
        Parameters
        ----------
        m_init : ndarray [Msun]
            inital masses for the stars
        met : ndarray [dex]
            initial metallicities for the stars
        age : ndarray [Gyr]
            ages for the stars
        props : Set
            list of stellar properties that should be returnd
        inter_age: str ["linear" or "log"]
            performe the age interpolation in age or in log10(age).

        Returns
        -------
        s_track : dict
        in_grid : ndarray
        """

        # calculate log_ages
        log_ages = np.log10(age) + 9

        # get adjacent grid-points
        (iso_met_lower, iso_met_higher,
        log_iso_ages_lower, log_iso_ages_higher,
        iso_ages_lower, iso_ages_higher) = self.get_adjacent_gridpoints(met, log_ages)

        # get weights:
        if inter_age == "linear":
            w1, w2, w3, w4 = self.get_weights(met, age, iso_met_lower, iso_met_higher,
                iso_ages_lower, iso_ages_higher)
        elif inter_age == "log":
            w1, w2, w3, w4 = self.get_weights(met, age, iso_met_lower, iso_met_higher,
                log_iso_ages_lower, log_iso_ages_higher)
        else:
            raise ValueError("inter_age must either be log or linear")

        # split properties
        props_no_charon = props.intersection(self.props_no_charon)
        props_with_charon = props.difference(self.props_no_charon)

        in_grid = None

        # evolve props not_good for charon:
        if props_no_charon:
            p_f1, in_grid1, max_mass1 = self.interp_props(m_init, iso_met_lower, log_iso_ages_lower,
                props_no_charon)
            p_f2, in_grid2, max_mass2 = self.interp_props(m_init, iso_met_lower, log_iso_ages_higher,
                props_no_charon)
            p_f3, in_grid3, max_mass3 = self.interp_props(m_init, iso_met_higher, log_iso_ages_lower,
                props_no_charon)
            p_f4, in_grid4, max_mass4 = self.interp_props(m_init, iso_met_higher, log_iso_ages_higher,
                props_no_charon)

            result = (w1 * p_f1 + w2 * p_f2 + w3 * p_f3 + w4 * p_f4) / (w1 + w2 + w3 + w4)

            in_grid = np.product([in_grid1, in_grid2, in_grid3, in_grid4])
        else:
            result = None

        if props_with_charon:
            props_with_charon.add('phase')
            props_with_charon = list(props_with_charon)
            phase_index = props_with_charon.index('phase')

            # props good for charon:
            mass1 = self.get_modified_mass(m_init, met,
                log_ages, iso_met_lower, log_iso_ages_lower)
            mass2 = self.get_modified_mass(m_init, met,
                log_ages, iso_met_lower, log_iso_ages_higher)
            mass3 = self.get_modified_mass(m_init, met,
                log_ages, iso_met_higher, log_iso_ages_lower)
            mass4 = self.get_modified_mass(m_init, met,
                log_ages, iso_met_higher, log_iso_ages_higher)

            p_f1_charon, in_grid1, max_mass_1 = self.interp_props(
                mass1, iso_met_lower, log_iso_ages_lower, props_with_charon)

            p_f2_charon, in_grid2, max_mass_2 = self.interp_props(
                mass2, iso_met_lower, log_iso_ages_higher, props_with_charon)

            p_f3_charon, in_grid3, max_mass_3 = self.interp_props(
                mass3, iso_met_higher, log_iso_ages_lower, props_with_charon)

            p_f4_charon, in_grid4, max_mass_4 = self.interp_props(
                mass4, iso_met_higher, log_iso_ages_higher, props_with_charon)


            result_charon = ((w1 * p_f1_charon + w2 * p_f2_charon
                              + w3 * p_f3_charon + w4 * p_f4_charon)
                             / (w1 + w2 + w3 + w4))

            in_grid = (in_grid1 | (w1[:, 0] == 0)) \
                      & (in_grid2 | (w2[:, 0] == 0)) \
                      & (in_grid3 | (w3[:, 0] == 0)) \
                      & (in_grid4 | (w4[:, 0] == 0))

        else:
            mass1 = mass2 = mass3 = mass4 = m_init
            result_charon = None

        final_phase = ((m_init > self.tips2_func(met, log_ages, grid=False))
                    | (max_mass_1 < mass1)
                    | (max_mass_2 < mass2)
                    | (max_mass_3 < mass3)
                    | (max_mass_4 < mass4))

        output = {item: result[..., i] for i, item in enumerate(props_no_charon)}
        output.update({item: result_charon[..., i] for i, item in enumerate(props_with_charon)})

        return output, in_grid, final_phase
