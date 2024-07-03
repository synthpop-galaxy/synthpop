"""
This file contains constants and parameters which should be treated as constant
"""

import os
import numpy as np

# -------- Directories ---------
SYNTHPOP_DIR = os.path.abspath(os.path.dirname(__file__))  # this directory
DEFAULT_MODEL_DIR = f"{SYNTHPOP_DIR}/models"  # location where models are stored
DEFAULT_CONFIG_DIR = f"{SYNTHPOP_DIR}/config_files/"
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
        "pop", "iMass", "age",  "Fe/H_initial", "Mass", "In_Final_Phase", "Dist", "l", "b",
        "vr_bc", "mul", "mub",  "x", "y", "z",  "U", "V", "W", "VR_LSR", "ExtinctionInMap"]


# ------- Physical Constants --------
c = 3e8  # m/s
G = 6.674e-11  # m^3 kg^-1 s^-2
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

"""
Default values for the position and velocities of the SUN and LSR , 
are set in synthpop_utils.sun_info.SunInfo. 
Can be modified in the Config File

X_SUN = -8.178  # distance from galactic center (kpc)
Y_SUN = 0.0
Z_SUN = 0.017  # distance above galactic plane (kpc)

GAL_DIST = (X_SUN ** 2 + Z_SUN ** 2 + Y_SUN ** 2) ** 0.5  # distance to the galactic center (kpc)

THETA_SUN = np.arcsin(Z_SUN / GAL_DIST)
CT_SUN = np.cos(THETA_SUN)
ST_SUN = np.sin(THETA_SUN)

# rectangular motion of the sun
# from Reid & Brunthaler (2020)
U_SUN = 12.9  # [km/s] = -V_R
V_SUN = 245.6  # [km/s] = V_phi
W_SUN = 7.78  # [km/s] = V_z

# velocity of the local standard of rest
# from Schönrich et al. (2010)
U_LSR = 1.8  # [km/s]
V_LSR = 233.4  # [km/s]
W_LSR = 0.53  # [km/s]

# peculiar motion of the Sun relative to the LSR
V_PEC = np.sqrt((U_SUN - U_LSR) ** 2 + (W_SUN - W_LSR) ** 2 + (W_SUN - W_LSR) ** 2)

B_APX_DEG = 25  # degree direction of the solar apex in galactic coordinates
L_APX_DEG = 53  # degree direction of the solar apex in galactic coordinates
"""
