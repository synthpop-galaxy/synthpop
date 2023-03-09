""" This file contains several utils function """
__all__ = ['solidangle_to_half_cone_angle', 'half_cone_angle_to_solidangle']
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__license__ = "GPLv3"
__version__ = "1.0.0"


import numpy as np


def solidangle_to_half_cone_angle(solid_angle):
    return np.arccos(1 - solid_angle / (2. * np.pi))


def half_cone_angle_to_solidangle(cone_angle):
    return (2. * np.pi) * (1 - np.cos(cone_angle))
