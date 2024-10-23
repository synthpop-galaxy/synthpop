"""
Takes line of sight, galactic coordinates, distance from the sun, solid angle, and magnitude band as inputs
"""

import numpy as np

#Model view parameters

#Galactic longitude (deg)  (np.arange does not include the end value, so make it one step higher than needed)
l_set = [16.0]

#Galactic latitude (deg)
b_set = [-3.0]
#b_set = [0.0]

#Maximum distance from sun (kpc)
dMax = 10
#Step size, in parsec 
dStep=10

#Field of view solid angle (deg^2)
angleMethod = 1 #1 for angle given, 2 to find angle from stellar mass in region.
#solidAngle =[0.000001]
#solidAngle =[6.09234e-7]
solidAngle =[0.000000609234]
m_stars = 1000


#Limiting mass [min_mass, max_mass]
mass_lims = {'min_mass':0.01,'max_mass':100.}
#Limiting magnitude- [system, band, limit]
maglim = ['UBVRIplus','2MASS_J',1000]

#Magnitude systems and bands to generate for each star
magsys = "UBVRIplus"   #UVBRIplus, JWST, HST_WFPC2, or SDSSugriz for now
bandOpt = 1   #0 to use all bands in the iso file, 1 to specify below
chosenBands = ["2MASS_Ks", "2MASS_J", "2MASS_H", "Bessell_V","Bessell_I"]

#Set of populations to produce (numbers 1-10 inclusive)
populations = range(1,11)
starsPerPop = 100

#Bulge model
alt = 2

#1 for observed magnitude, 0 for absolute
obsmag = 1

#Angle of galactic bar from line of sight toward GC
barang = 29.4

#Output file name
outputLocation = 'outputfiles/calibration'
outputFileType = 'csv'  #for now, pickle or csv

#width of each column in output
columnWidth = 12

#column headers for output
cols = ["iMass", "Mass","pop", "Fe/H", "logL", "logTeff", "logg","logage"]

#constants
y_sun = -8.25 #distance from galactic center (kpc)
x_sun = 0.0 
z_sun = 0.015 #distance above galactic plane (kpc)
R_sun = 8.25 #radius of solar orbit around galactic center (kpc)
Vlsr = 220.0 #radial velocity of sun (km/s)

c = 3 * 10**8 #m/s
G = 6.674 * 10**-11  #m^3 kg^-1 s^-2
Mbol_sun = 4.74

isops = ['star_mass','log_Teff', 'log_g','log_L']
