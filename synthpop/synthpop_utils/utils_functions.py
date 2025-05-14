""" This file contains several utils function """
__all__ = ['solidangle_to_half_cone_angle', 'half_cone_angle_to_solidangle', "rotation_matrix"]
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]


import numpy as np


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
