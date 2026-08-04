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
    all_systems = df['system_idx'].to_numpy()
    idx_map = pd.Series(np.arange(len(all_systems)), index=all_systems)
    comp_rows = idx_map.loc[comp_df['system_idx']].to_numpy()

    comp_mags = np.ma.masked_invalid(comp_df[filters].to_numpy())
    comp_flux = (10.0 ** (-0.4 * comp_mags)).filled(0.0)

    main_mags = np.ma.masked_invalid(df[filters].to_numpy())
    system_flux = (10.0 ** (-0.4 * main_mags)).filled(0.0)
    np.add.at(system_flux, comp_rows, comp_flux)

    masked_system_flux = np.ma.masked_less_equal(system_flux, 0.0)
    df[filters] = (-2.5 * np.ma.log10(masked_system_flux)).filled(np.nan)

    return df


def get_primary_mags(df, comp_df, filters):
    primary_idxs = (df['n_companions'] > 0).to_numpy()
    if not np.any(primary_idxs):
        return df

    primary_systems = df.loc[primary_idxs, 'system_idx'].to_numpy()
    idx_map = pd.Series(np.arange(len(primary_systems)), index=primary_systems)
    comp_rows = idx_map.loc[comp_df['system_idx']].to_numpy()
    comp_mags = np.ma.masked_invalid(comp_df[filters].to_numpy())
    comp_flux = (10.0 ** (-0.4 * comp_mags)).filled(0.0)

    system_mags = np.ma.masked_invalid(df.loc[primary_idxs, filters].to_numpy())
    primary_flux = (10.0 ** (-0.4 * system_mags)).filled(0.0)
    np.add.at(primary_flux, comp_rows, -comp_flux)

    masked_primary_flux = np.ma.masked_less_equal(primary_flux, 0.0)
    df.loc[primary_idxs, filters] = (-2.5 * np.ma.log10(masked_primary_flux)).filled(np.nan)

    return df

