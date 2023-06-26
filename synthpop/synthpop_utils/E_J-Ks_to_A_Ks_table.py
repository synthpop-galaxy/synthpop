import numpy as np
import os
import pandas as pd
import requests
import urllib.request

"""
This script automates the download and conversion of extinction map values by Surot (2020) (as found at https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/644/A140) from reddening in the J-Ks band to extinction in the A_Ks band. The values are reformatted to follow the structure of marshall_table1.csv and saved to the extinction directory, where they can be used by the surot.py module.

Columns of Surot map data:
    [0] GLON (deg): []Galactic longitude of the node
    [1] GLAT (deg): Galactic latitude of the node
    [2] E(J-Ks) (mag): Reddening
    [3] e_E(J-Ks) (mag): Error on `E(J-Ks)`
    [4] res-avg (arcmin): Average resolution
    [5] res-max (arcmin): Maximum resolution
    [6] Tile: Field ID following VVV designation

Columns of Marshall map data:
    [0] GLON (deg): [0/360] Galactic longitude 
    [1] GLAT (deg): Galactic latitude
    [2] x2all: Chi-squared, all stars 
    [3] x2gts: Chi-squared, giants
    [4] nb: [1,33] Number of bins used
    [5] r1 (kpc): Heliocentric distance of bin1
    [6] ext1 (mag): Ks band extinction of bin1
    [7] e_r1 (kpc): Distance uncertainty of bin1
    [8] e_ext1 (mag): Extinction uncertainty of bin1
    Columns [5,8] repeat for each bin, with a maximum of 33 sets, variable for each value of [0,1]
"""

def a_b_coeff(x):
    '''
    Wavelength-dependent extinction law coefficients as reported by CCM 1989
    Returns a(x), b(x)
    Valid for 0.3 um^-1 <= x <= 1.1 um^-1
    '''
    return 0.574*(x**1.61), -0.527*(x**1.61)

def ext_law(x, R_V):
    '''
    Average R_V-dependent extinction law; returns <A(lambda)/A(V)>
    Valid for 0.3 um^-1 <= x <= 1.1 um^-1
    '''
    I_a, I_b = a_b_coeff(x)
    return I_a + I_b/R_V

def E_to_A(E, x_J, x_K):
    '''
    Converts reddending values in the J - K_s band to extinction values A_K
    Wavenumbers x_J and x_K are in units of inverse microns
    '''
    R_V = 2.5 # reddening parameter approximated for the bulge
    J_law = ext_law(x_J, R_V)
    K_law = ext_law(x_K, R_V)
    return (K_law / (J_law - K_law)) * E

pwd = os.getcwd()
extinction_dir = pwd+'/../modules/extinction'

if os.path.isfile(extinction_dir+'/surot_A_Ks_table1.csv'):
    print('Converted map already exists in the extinction directory.')
    quit()

print('Downloading map file from VizieR...')
map_url = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/644/A140/ejkmap.dat.gz'
map_filename = map_url.split("/")[-1]
with open(map_filename, "wb") as f:
    r = requests.get(map_url)
    f.write(r.content)
    print('Map retrieved.')

print('Reading table...')
E_JKs_map_df = pd.read_fwf(map_filename, header=None, compression='gzip')

print('Reformatting values...')
E_JKs_map = E_JKs_map_df.to_numpy()

# reference wavelength of each filter in micrometers
# source: http://svo2.cab.inta-csic.es/theory/fps/index.php?id=Paranal/VISTA.Ks&&mode=browse&gname=Paranal&gname2=VISTA#filter 
wl_ref_J = 12524.83 * 1e-4
wl_ref_Ks = 21521.52 * 1e-4 

# wavenumber corresponding to the central wavelength of each filter
x_J = 1/wl_ref_J
x_K = 1/wl_ref_Ks
A_Ks_vals = E_to_A(E_JKs_map[:,2], x_J, x_K)

entries = E_JKs_map.shape[0]

surot_reformat = np.zeros((entries, 21))
surot_reformat[:,0] = E_JKs_map[:,0]
surot_reformat[:,1] = E_JKs_map[:,0]
surot_reformat[:,2] = np.zeros(entries)
surot_reformat[:,3] = np.zeros(entries)
surot_reformat[:,4] = 5*np.ones(entries)

# identical extinctions at 4 arbitrary distance bin intervals
for i in range(4):
    surot_reformat[:,(5 + (4*i))] = (0.5 + i*1.5)*np.ones(entries)
    surot_reformat[:,(6 + (4*i))] = A_Ks_vals
    surot_reformat[:,(7 + (4*i))] = np.zeros(entries)
    surot_reformat[:,(8 + (4*i))] = np.zeros(entries)

surot_reformat_df = pd.DataFrame(surot_reformat, columns=None)
surot_reformat_df[4] = surot_reformat_df[4].astype('int')

print('Saving table to extinction directory...')
map_output = 'surot_A_Ks_table1.csv'
surot_reformat_df.to_csv(extinction_dir+'/'+map_output, header=False, index=False)

print('File saved as '+map_output)