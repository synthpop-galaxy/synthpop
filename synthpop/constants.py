"""
This file contains constants and parameters which should be treated as constant
"""

import os
import numpy as np

# -------- Directories ---------
SYNTHPOP_DIR = os.path.abspath(os.path.dirname(__file__))  # this directory
DEFAULT_MODEL_DIR = f"{SYNTHPOP_DIR}/models"  # location where models are stored
DEFAULT_CONFIG_DIR = f"{SYNTHPOP_DIR}/config_files"
DEFAULT_CONFIG_FILE = f"{DEFAULT_CONFIG_DIR}/_default.synthpop_conf"  # default config files
ISOCHRONES_DIR = f"{SYNTHPOP_DIR}/data/isochrones"  # location where isochrones are stored
EXTINCTIONS_DIR = f"{SYNTHPOP_DIR}/data/extinction"  # location where isochrones are stored
MOMENTS_DIR = f"{SYNTHPOP_DIR}/data/moments"  # location where moment grids are stored

# ------- Data Columns --------
# parameter which needs to be estimated from the isochrones
# the first column MUST be the evolved stellar mass!
# append all further columns at the end of cols
# Note that the Isochrones also need the initial mass, age, and metallicity for the interpolation
REQ_ISO_PROPS = ["star_mass"]

# Columnnames for output table
# You can rename each column, but you must keep the order!
# If you add further items to REQ_ISO_PROPS, you must add the column names at the end.
COL_NAMES = [
        "pop", "iMass", "age",  "Fe/H_initial", "Mass",
        "In_Final_Phase", "Dist", "l", "b",
        "vr_bc", "mul", "mub",  "x", "y", "z",  "U", "V", "W", "VR_LSR", "ExtinctionInMap"]

# ------- Physical Constants --------
c = 3e8  # m/s
G = 6.674e-11  # m^3 kg^-1 s^-2
Msun_kg = 1.989e30 # solar mass in kg
m_per_pc = 3.086e16 # meters per parsec
Mbol_sun = 4.74  # Bolometric magnitude of the sun
DEG2RAD = np.pi / 180
# converts proper motion times distance in mas/yr and kpc to tangential velocity in km/s
MUxD_to_VT = 4.740470446

# -------- Coordinate system --------
# Coordinates of the galactic North Pole in ecliptic coordinates
A_NGP_DEG = 192.8583
D_NGP_DEG = 27.1283
# Galactic longitude of the North Celestial Pole.
L_NCP_DEG = 122.93192
