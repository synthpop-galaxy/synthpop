"""
This file contains several utility functions.
"""

__all__ = ["solidangle_to_half_cone_angle", "half_cone_angle_to_solidangle",
            "rotation_matrix", "combine_system_mags", "get_primary_mags"]
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]

import numpy as np
import pandas as pd
import warnings
import pdb

def solidangle_to_half_cone_angle(solid_angle):
    return np.arccos(1 - solid_angle / (2. * np.pi))

def half_cone_angle_to_solidangle(cone_angle):
    return (2. * np.pi) * (1 - np.cos(cone_angle))

def rotation_matrix(
        theta_rad: float or np.ndarray or None = None,
        st: float or np.ndarray or None = None, ct: float or np.ndarray or None = None,
        axis: str = ''
        ) -> np.ndarray:
    """
    creates the rotation matrix along a given axis.
    These are:

         | 1    0    0  |
    RX = | 0   cos -sin |
         | 0   sin  cos |

         | cos  0  -sin |
    RY = |  0   1    0  |
         | sin  0   cos |

         | cos -sin   0  |
    RZ = | sin  cos   0  |
         |  0    0    1  |

    Parameters
    ----------
    theta_rad: float or np.ndarray or None [rad]
        rotation angle
    st, ct: float or np.ndarray or None [rad]
        sin(theta) and cos(theta)
    axis: str
        string to specify the axis. either 'x', 'y', or 'z'

    Returns
    -------
    rot_matrix : np.ndarray
        rotation matrix
    """

    if axis not in ['x', 'y', 'z']:
        raise ValueError("axis has to be 'x', 'y', or 'z'!")
    if theta_rad is None and (st is None or ct is None):
        raise ValueError("either theta or ct and st has to be specified")

    if st is None or ct is None:
        st = np.sin(theta_rad)
        ct = np.cos(theta_rad)

    if isinstance(st, np.ndarray):
        one = np.ones(st.shape)
        zero = np.zeros(st.shape)
    else:
        one = 1
        zero = 0

    if axis == 'x':
        return np.array([[one, zero, zero], [zero, ct, -st], [zero, st, ct]])
    if axis == 'y':
        return np.array([[ct, zero, -st], [zero, one, zero], [st, zero, ct]])
    if axis == 'z':
        return np.array([[ct, -st, zero], [st, ct, zero], [zero, zero, one]])


def combine_system_mags(df, comp_df, filters):
    sys_indexer = pd.Index(df['system_idx'])
    comp_rows = sys_indexer.get_indexer(comp_df['system_idx'])

    comp_mags = comp_df[filters].to_numpy()
    main_mags = df[filters].to_numpy()
    main_mags_for_comp = main_mags[comp_rows]

    # Companion-to-primary relative flux ratio: 10^(-0.4 * (m_comp - m_main))
    delta_m = comp_mags - main_mags_for_comp
    comp_to_main_ratio = np.nan_to_num(10.0 ** (-0.4 * delta_m), nan=0.0)

    # Accumulate relative companion flux ratio per primary system
    total_comp_ratio = np.zeros_like(main_mags)
    np.add.at(total_comp_ratio, comp_rows, comp_to_main_ratio)

    # Accumulate absolute companion flux for dark primaries
    comp_abs_flux = np.nan_to_num(10.0 ** (-0.4 * comp_mags), nan=0.0)
    total_comp_abs_flux = np.zeros_like(main_mags)
    np.add.at(total_comp_abs_flux, comp_rows, comp_abs_flux)

    inv_ln10_25 = 2.5 / np.log(10.0)

    with np.errstate(divide='ignore', invalid='ignore'):
        # Precision path using log1p: m_tot = m_main - (2.5/ln10) * ln(1 + ratio)
        combined_with_primary = main_mags - inv_ln10_25 * np.log1p(total_comp_ratio)
        
        # Fallback path for dark primaries (where primary magnitude is NaN)
        combined_dark = np.where(total_comp_abs_flux > 0.0, -2.5 * np.log10(total_comp_abs_flux), np.nan)

        combined_mags = np.where(~np.isnan(main_mags), combined_with_primary, combined_dark)

    df[filters] = combined_mags
    return df

def get_primary_mags(df, comp_df, filters, rtol=1e-12):
    # remove nan mag companions
    comp_clean = comp_df.dropna(subset=[filters[0]])
    if comp_clean.empty:
        return df
    comp_rows = pd.Index(df['system_idx']).get_indexer(comp_clean['system_idx'])

    # magnitude math
    comp_mags = comp_clean[filters].to_numpy()
    sys_mags = df[filters].to_numpy()
    sys_mags_for_comp = sys_mags[comp_rows]
    delta_m = comp_mags - sys_mags_for_comp
    flux_ratio = 10.0 ** (-0.4 * delta_m)
    total_comp_ratio = np.zeros_like(sys_mags)
    np.add.at(total_comp_ratio, comp_rows, flux_ratio)
    inv_ln10_25 = 2.5 / np.log(10.0)

    # mag subtraction
    with np.errstate(invalid='ignore', divide='ignore'):
        remaining_fraction = 1.0 - total_comp_ratio
        valid_primary = (remaining_fraction > rtol) & ~np.isnan(sys_mags)

        mag_correction = -inv_ln10_25 * np.log1p(-total_comp_ratio)
        primary_mags = np.where(valid_primary, sys_mags + mag_correction, np.nan)

    df[filters] = primary_mags
    return df

