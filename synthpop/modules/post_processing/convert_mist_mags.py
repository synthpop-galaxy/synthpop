"""
Postprocessing module to convert magnitude systems for any
filters provided by MIST.
"""

__all__ = ["ConvertMistMags", ]
__author__ = "M.J. Huston"
__date__ = "2024-02-23"

import pandas
import numpy as np
from ._post_processing import PostProcessing

class ConvertMistMags(PostProcessing):
    """
    Postprocessing module to convert magnitude systems for any
    filters provided by MIST. Allowed systems are Vega, AB, and ST.
    
    Attributes
    ----------
    conversions : dict
        dictionary defining which filters to convert to which systems;
        format: {'NEW_SYSTEM':['FILTER1', 'FILTER2']}
    """

    def __init__(self, model, conversions, logger, **kwargs):
        super().__init__(model,logger, **kwargs)
        self.conversions = conversions

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        Perform the magnitude conversions and returns the modified DataFrame.
        """

        # MIST data
        # TODO find a better way to do this from the file
        #  - requires dealing with Roman filter naming discrepancy
        mist_filters = np.array(['CFHT_u', 'CFHT_CaHK', 'CFHT_g', 'CFHT_r', 'CFHT_i_new',
           'CFHT_i_old', 'CFHT_z', 'DECam_u', 'DECam_g', 'DECam_r', 'DECam_i',
           'DECam_z', 'DECam_Y', 'GALEX_FUV', 'GALEX_NUV', 'ACS_HRC_F220W',
           'ACS_HRC_F250W', 'ACS_HRC_F330W', 'ACS_HRC_F344N', 'ACS_HRC_F435W',
           'ACS_HRC_F475W', 'ACS_HRC_F502N', 'ACS_HRC_F550M', 'ACS_HRC_F555W',
           'ACS_HRC_F606W', 'ACS_HRC_F625W', 'ACS_HRC_F658N', 'ACS_HRC_F660N',
           'ACS_HRC_F775W', 'ACS_HRC_F814W', 'ACS_HRC_F850LP',
           'ACS_HRC_F892N', 'ACS_WFC_F435W', 'ACS_WFC_F475W', 'ACS_WFC_F502N',
           'ACS_WFC_F550M', 'ACS_WFC_F555W', 'ACS_WFC_F606W', 'ACS_WFC_F625W',
           'ACS_WFC_F658N', 'ACS_WFC_F660N', 'ACS_WFC_F775W', 'ACS_WFC_F814W',
           'ACS_WFC_F850LP', 'ACS_WFC_F892N', 'WFC3_UVIS_F200LP',
           'WFC3_UVIS_F218W', 'WFC3_UVIS_F225W', 'WFC3_UVIS_F275W',
           'WFC3_UVIS_F280N', 'WFC3_UVIS_F300X', 'WFC3_UVIS_F336W',
           'WFC3_UVIS_F343N', 'WFC3_UVIS_F350LP', 'WFC3_UVIS_F373N',
           'WFC3_UVIS_F390M', 'WFC3_UVIS_F390W', 'WFC3_UVIS_F395N',
           'WFC3_UVIS_F410M', 'WFC3_UVIS_F438W', 'WFC3_UVIS_F467M',
           'WFC3_UVIS_F469N', 'WFC3_UVIS_F475W', 'WFC3_UVIS_F475X',
           'WFC3_UVIS_F487N', 'WFC3_UVIS_F502N', 'WFC3_UVIS_F547M',
           'WFC3_UVIS_F555W', 'WFC3_UVIS_F600LP', 'WFC3_UVIS_F606W',
           'WFC3_UVIS_F621M', 'WFC3_UVIS_F625W', 'WFC3_UVIS_F631N',
           'WFC3_UVIS_F645N', 'WFC3_UVIS_F656N', 'WFC3_UVIS_F657N',
           'WFC3_UVIS_F658N', 'WFC3_UVIS_F665N', 'WFC3_UVIS_F673N',
           'WFC3_UVIS_F680N', 'WFC3_UVIS_F689M', 'WFC3_UVIS_F763M',
           'WFC3_UVIS_F775W', 'WFC3_UVIS_F814W', 'WFC3_UVIS_F845M',
           'WFC3_UVIS_F850LP', 'WFC3_UVIS_F953N', 'WFC3_IR_F098M',
           'WFC3_IR_F105W', 'WFC3_IR_F110W', 'WFC3_IR_F125W', 'WFC3_IR_F126N',
           'WFC3_IR_F127M', 'WFC3_IR_F128N', 'WFC3_IR_F130N', 'WFC3_IR_F132N',
           'WFC3_IR_F139M', 'WFC3_IR_F140W', 'WFC3_IR_F153M', 'WFC3_IR_F160W',
           'WFC3_IR_F164N', 'WFC3_IR_F167N', 'WFPC2_F218W', 'WFPC2_F255W',
           'WFPC2_F300W', 'WFPC2_F336W', 'WFPC2_F439W', 'WFPC2_F450W',
           'WFPC2_F555W', 'WFPC2_F606W', 'WFPC2_F622W', 'WFPC2_F675W',
           'WFPC2_F791W', 'WFPC2_F814W', 'WFPC2_F850LP', 'INT_IPHAS_gR',
           'INT_IPHAS_Ha', 'INT_IPHAS_gI', 'JWST_F070W', 'JWST_F090W',
           'JWST_F115W', 'JWST_F140M', 'JWST_F150W2', 'JWST_F150W',
           'JWST_F162M', 'JWST_F164N', 'JWST_F182M', 'JWST_F187N',
           'JWST_F200W', 'JWST_F210M', 'JWST_F212N', 'JWST_F250M',
           'JWST_F277W', 'JWST_F300M', 'JWST_F322W2', 'JWST_F323N',
           'JWST_F335M', 'JWST_F356W', 'JWST_F360M', 'JWST_F405N',
           'JWST_F410M', 'JWST_F430M', 'JWST_F444W', 'JWST_F460M',
           'JWST_F466N', 'JWST_F470N', 'JWST_F480M', 'LSST_u', 'LSST_g',
           'LSST_r', 'LSST_i', 'LSST_z', 'LSST_y', 'PS_g', 'PS_r', 'PS_i',
           'PS_z', 'PS_y', 'PS_w', 'PS_open', 'SDSS_u', 'SDSS_g', 'SDSS_r',
           'SDSS_i', 'SDSS_z', 'SkyMapper_u', 'SkyMapper_v', 'SkyMapper_g',
           'SkyMapper_r', 'SkyMapper_i', 'SkyMapper_z', 'IRAC_3.6',
           'IRAC_4.5', 'IRAC_5.8', 'IRAC_8.0', 'hsc_g', 'hsc_r', 'hsc_i',
           'hsc_z', 'hsc_y', 'hsc_nb816', 'hsc_nb921', 'Swift_UVW2',
           'Swift_UVM2', 'Swift_UVW1', 'Swift_U', 'Swift_B', 'Swift_V',
           'Bessell_U', 'Bessell_B', 'Bessell_V', 'Bessell_R', 'Bessell_I',
           '2MASS_J', '2MASS_H', '2MASS_Ks', 'Kepler_Kp', 'Kepler_D51',
           'Hipparcos_Hp', 'Tycho_B', 'Tycho_V', 'Gaia_G_DR2Rev',
           'Gaia_BP_DR2Rev', 'Gaia_RP_DR2Rev', 'TESS', 'Gaia_G_EDR3',
           'Gaia_BP_EDR3', 'Gaia_RP_EDR3', 'UKIDSS_Z', 'UKIDSS_Y', 'UKIDSS_J',
           'UKIDSS_H', 'UKIDSS_K', 'VISTA_Z', 'VISTA_Y', 'VISTA_J', 'VISTA_H',
           'VISTA_Ks', 'Washington_C', 'Washington_M', 'Washington_T1',
           'Washington_T2', 'DDO51_vac', 'DDO51_f31', 'Stromgren_u',
           'Stromgren_v', 'Stromgren_b', 'Stromgren_y', 'R062',
           'Z087', 'Y106', 'J129', 'W146',
           'H158', 'F184', 'WISE_W1', 'WISE_W2', 'WISE_W3',
           'WISE_W4', 'SPLUS_uJAVA', 'SPLUS_gSDSS', 'SPLUS_rSDSS',
           'SPLUS_iSDSS', 'SPLUS_zSDSS', 'SPLUS_J0340', 'SPLUS_J0378',
           'SPLUS_J0395', 'SPLUS_J0410', 'SPLUS_J0515', 'SPLUS_J0660',
           'SPLUS_J0861', 'UVIT_F148W', 'UVIT_F154W', 'UVIT_F169M',
           'UVIT_F172M', 'UVIT_F242W', 'UVIT_N219M', 'UVIT_N245M',
           'UVIT_N263M', 'UVIT_N279N'])
        mist_systems = np.array(['AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB',
           'AB', 'AB', 'AB', 'AB', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'AB', 'AB',
           'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB',
           'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB',
           'Vega', 'Vega', 'Vega', 'Vega', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB',
           'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega', 'Vega',
           'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB',
           'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB', 'AB'])
        mist_mag_vega_ab = np.array([ 3.671250e-01, -2.784500e-02, -9.181000e-02,  1.543760e-01,
            3.616700e-01,  3.854910e-01,  5.138370e-01,  3.408350e-01,
           -8.959600e-02,  1.824670e-01,  4.058980e-01,  5.169830e-01,
            5.643160e-01,  2.123122e+00,  1.665903e+00,  1.683496e+00,
            1.495635e+00,  1.182286e+00,  1.151053e+00, -8.659100e-02,
           -9.393600e-02, -8.359600e-02,  2.428300e-02, -7.055000e-03,
            8.079100e-02,  1.601160e-01,  3.765720e-01,  2.788630e-01,
            3.833350e-01,  4.255120e-01,  5.214850e-01,  4.875050e-01,
           -1.016790e-01, -9.787500e-02, -8.368100e-02,  2.451100e-02,
           -6.312000e-03,  8.650000e-02,  1.632580e-01,  3.763820e-01,
            2.787100e-01,  3.879550e-01,  4.236630e-01,  5.195860e-01,
            4.888930e-01,  4.513330e-01,  1.691151e+00,  1.660654e+00,
            1.498659e+00,  1.432770e+00,  1.420732e+00,  1.185760e+00,
            1.155674e+00,  1.567650e-01,  8.952060e-01,  9.930700e-02,
            2.159570e-01, -2.098700e-02, -1.559170e-01, -1.520480e-01,
           -1.547330e-01, -1.572470e-01, -9.709200e-02, -4.930700e-02,
            1.818820e-01, -8.693600e-02, -3.710000e-04, -2.442000e-02,
            3.240820e-01,  8.301300e-02,  1.461920e-01,  1.480780e-01,
            1.577680e-01,  1.895600e-01,  6.204030e-01,  3.255420e-01,
            3.620800e-01,  2.366110e-01,  2.388950e-01,  2.573280e-01,
            2.768730e-01,  3.788000e-01,  3.807460e-01,  4.182570e-01,
            5.002520e-01,  5.208000e-01,  6.141610e-01,  5.615530e-01,
            6.453020e-01,  7.594890e-01,  9.010520e-01,  9.211300e-01,
            9.612040e-01,  1.037116e+00,  9.761090e-01,  9.973750e-01,
            1.078618e+00,  1.076328e+00,  1.253686e+00,  1.251445e+00,
            1.384753e+00,  1.361575e+00,  1.696001e+00,  1.565908e+00,
            1.336668e+00,  1.183504e+00, -1.455440e-01, -8.495400e-02,
           -1.464000e-03,  1.001710e-01,  1.413580e-01,  2.396900e-01,
            4.107450e-01,  4.168170e-01,  5.187360e-01,  1.467310e-01,
            3.586730e-01,  3.938720e-01,  2.778270e-01,  5.114690e-01,
            7.777270e-01,  1.103238e+00,  1.228996e+00,  1.205652e+00,
            1.351846e+00,  1.391996e+00,  1.558015e+00,  1.624210e+00,
            1.675761e+00,  1.780225e+00,  1.801134e+00,  2.116563e+00,
            2.295891e+00,  2.457418e+00,  2.533096e+00,  2.606257e+00,
            2.680412e+00,  2.781285e+00,  2.830344e+00,  3.081301e+00,
            3.072036e+00,  3.170645e+00,  3.207505e+00,  3.333257e+00,
            3.375634e+00,  3.363519e+00,  3.407118e+00,  6.584550e-01,
           -9.177600e-02,  1.450080e-01,  3.636270e-01,  5.075140e-01,
            5.436580e-01, -8.780100e-02,  1.435310e-01,  3.632260e-01,
            5.071250e-01,  5.390140e-01,  1.259080e-01,  1.994100e-01,
            9.310960e-01, -1.005920e-01,  1.424690e-01,  3.558530e-01,
            5.179630e-01,  1.104895e+00,  3.119500e-01, -5.823300e-02,
            1.278690e-01,  4.008760e-01,  5.259820e-01,  2.785147e+00,
            3.258136e+00,  3.750716e+00,  4.391595e+00, -9.538100e-02,
            1.438140e-01,  3.864800e-01,  5.139290e-01,  5.484350e-01,
            4.695170e-01,  5.408240e-01,  1.734154e+00,  1.686980e+00,
            1.510161e+00,  1.013124e+00, -1.158850e-01, -4.726000e-03,
            8.005270e-01, -1.075120e-01,  6.521000e-03,  1.902780e-01,
            4.313720e-01,  8.891760e-01,  1.364157e+00,  1.834505e+00,
            1.226660e-01, -5.223600e-02,  7.642000e-03, -6.085500e-02,
           -1.422800e-02,  1.238910e-01,  6.559000e-02,  3.686690e-01,
            3.637780e-01,  1.289610e-01,  3.023300e-02,  3.733490e-01,
            5.131270e-01,  6.148780e-01,  9.152450e-01,  1.352530e+00,
            1.871458e+00,  5.659320e-01,  6.046970e-01,  9.207850e-01,
            1.360424e+00,  1.817194e+00,  3.324540e-01, -5.143600e-02,
            1.857420e-01,  4.499720e-01, -5.757900e-02, -5.997500e-02,
            1.135368e+00, -1.321360e-01, -1.529660e-01,  2.927000e-03,
            1.370950e-01,  4.873790e-01,  6.537800e-01,  9.583630e-01,
            1.024467e+00,  1.287404e+00,  1.551332e+00,  2.665543e+00,
            3.305247e+00,  5.139422e+00,  6.614602e+00,  1.095256e+00,
           -9.441200e-02,  1.517590e-01,  3.837210e-01,  5.151610e-01,
           -1.326650e-01,  5.846980e-01,  2.552300e-02, -1.544190e-01,
           -5.977600e-02,  3.051520e-01,  5.300250e-01,  2.301383e+00,
            2.318347e+00,  2.014355e+00,  1.913212e+00,  1.628923e+00,
            1.693668e+00,  1.662303e+00,  1.538884e+00,  1.482277e+00])
        mist_mag_vega_st = np.array([-4.2453900e-01, -7.3632100e-01, -3.5765100e-01,  4.4090000e-01,
            1.0558410e+00,  1.1197090e+00,  1.5589800e+00, -4.2018400e-01,
           -3.6639600e-01,  5.2909300e-01,  1.1762060e+00,  1.6340430e+00,
            1.8431400e+00, -6.3168400e-01, -2.1647100e-01, -2.4251600e-01,
           -2.6944000e-02,  1.2369300e-01,  1.3786700e-01, -6.0585800e-01,
           -3.9098500e-01, -2.7089000e-01,  6.5311000e-02, -5.5068000e-02,
            2.3804800e-01,  4.6304000e-01,  7.7688800e-01,  6.8420400e-01,
            1.1136860e+00,  1.2799880e+00,  1.6351120e+00,  1.5463410e+00,
           -6.1269900e-01, -4.0860600e-01, -2.7092600e-01,  6.6147000e-02,
           -5.2478000e-02,  2.5580100e-01,  4.7170100e-01,  7.7677400e-01,
            6.8418100e-01,  1.1264650e+00,  1.2632230e+00,  1.6075170e+00,
            1.5474340e+00,  2.4041900e-01, -2.6088300e-01, -1.5517800e-01,
           -2.8405000e-02,  2.3350000e-03, -1.8331000e-02,  1.2204900e-01,
            1.4342700e-01,  3.0319200e-01,  6.1798000e-02, -6.3893600e-01,
           -5.0721700e-01, -7.2719800e-01, -7.7929600e-01, -6.6343100e-01,
           -4.9436700e-01, -4.9431300e-01, -3.9489200e-01, -2.7209800e-01,
           -7.1910000e-02, -2.7995500e-01, -1.1458000e-02, -9.1758000e-02,
            9.9295800e-01,  2.4057600e-01,  4.2261000e-01,  4.3239600e-01,
            4.6392700e-01,  5.4651500e-01,  1.0130880e+00,  7.2020200e-01,
            7.6242800e-01,  6.6058900e-01,  6.9847600e-01,  7.5241900e-01,
            7.7160100e-01,  1.0946900e+00,  1.1065670e+00,  1.2496560e+00,
            1.4392770e+00,  1.6403220e+00,  1.8176630e+00,  1.8399290e+00,
            2.0697290e+00,  2.3774240e+00,  2.6911290e+00,  2.7283290e+00,
            2.7950470e+00,  2.8865080e+00,  2.8547220e+00,  2.9061790e+00,
            3.0918720e+00,  3.1029280e+00,  3.4882180e+00,  3.4926500e+00,
            3.7673750e+00,  3.7755000e+00, -2.7839700e-01, -5.0464000e-02,
            2.6125000e-02,  1.2263000e-01, -6.6436800e-01, -4.8428900e-01,
           -1.5900000e-02,  2.9772700e-01,  4.0646300e-01,  6.8364000e-01,
            1.2014840e+00,  1.2435290e+00,  1.6287460e+00,  4.2295500e-01,
            7.5381800e-01,  1.1435470e+00,  8.0711300e-01,  1.5965320e+00,
            2.3978710e+00,  3.1502330e+00,  3.6382580e+00,  3.3954630e+00,
            3.7172150e+00,  3.7802540e+00,  4.1961670e+00,  4.2958910e+00,
            4.4761470e+00,  4.6944070e+00,  4.7419690e+00,  5.4172440e+00,
            5.8182820e+00,  6.1477380e+00,  6.3857630e+00,  6.4648860e+00,
            6.6200670e+00,  6.8458480e+00,  6.9330090e+00,  7.4272330e+00,
            7.4340690e+00,  7.6364140e+00,  7.7334860e+00,  7.9690360e+00,
            8.0227370e+00,  8.0352300e+00,  8.1287120e+00, -2.1306600e-01,
           -3.7408000e-01,  4.1836300e-01,  1.0574610e+00,  1.5094480e+00,
            1.7864730e+00, -3.5154200e-01,  4.1384500e-01,  1.0565780e+00,
            1.5062110e+00,  1.7645760e+00,  4.2568800e-01,  7.1531000e-01,
           -5.8240000e-03, -4.3100000e-01,  4.0380800e-01,  1.0362040e+00,
            1.5842270e+00,  1.4326700e-01, -4.6078900e-01, -2.2300200e-01,
            3.7612000e-01,  1.1603680e+00,  1.6401640e+00,  6.8484640e+00,
            7.8331230e+00,  8.8551300e+00,  1.0192009e+01, -3.8523200e-01,
            4.1863100e-01,  1.1237240e+00,  1.6006750e+00,  1.8062330e+00,
            1.3403760e+00,  1.6710340e+00, -3.9078700e-01, -2.4710200e-01,
           -1.2167600e-01,  2.0865000e-02, -6.1572700e-01, -2.4917000e-02,
           -1.1203500e-01, -5.9422100e-01,  1.1340000e-02,  5.6846200e-01,
            1.2496360e+00,  2.6566860e+00,  3.7538420e+00,  4.8150350e+00,
            4.2847300e-01, -1.7811900e-01, -4.3887000e-02, -6.4199600e-01,
           -8.4969000e-02,  4.1011400e-01, -1.1059100e-01,  1.1205610e+00,
            1.1543000e+00,  4.0502200e-01, -1.1983700e-01,  1.1331290e+00,
            1.5499540e+00,  1.9900890e+00,  2.7080050e+00,  3.7294790e+00,
            4.8974220e+00,  1.8362330e+00,  2.0081160e+00,  2.7464110e+00,
            3.7484430e+00,  4.7822630e+00, -4.3285000e-01, -2.0564600e-01,
            5.3108700e-01,  1.2973270e+00, -1.9297400e-01, -2.0044200e-01,
            1.5028000e-01, -7.5607200e-01, -5.0066200e-01, -3.5780000e-03,
            4.2970200e-01,  1.4956060e+00,  2.0907050e+00,  2.8270250e+00,
            3.1324160e+00,  3.5898960e+00,  4.1892780e+00,  6.6104970e+00,
            7.9354230e+00,  1.1856446e+01,  1.4653819e+01,  1.4411200e-01,
           -3.9911900e-01,  4.3972500e-01,  1.1158290e+00,  1.5801570e+00,
           -6.6136400e-01, -2.2384300e-01, -6.8866100e-01, -7.8523800e-01,
           -1.9990700e-01,  7.1537300e-01,  1.5122920e+00, -5.1147800e-01,
           -5.0769500e-01, -6.5505100e-01, -6.0514500e-01, -1.4936100e-01,
           -2.7566700e-01, -8.6757000e-02, -4.0787000e-02,  1.9904000e-02])
            
        for system_new in self.conversions.keys():
            for filt in self.conversions[system_new]:
                if filt in dataframe.keys():
                    idx = np.where(mist_filters==filt)[0][0]
                    system_old = mist_systems[idx]
                    if system_old==system_new:
                        print(filt,'already in system',system_new)
                    elif system_old=='Vega' and system_new=='AB':
                        dataframe[filt] = dataframe[filt] + mist_mag_vega_ab[idx]
                    elif system_old=='Vega' and system_new=='ST':
                        dataframe[filt] = dataframe[filt] + mist_mag_vega_st[idx]
                    elif system_old=='AB' and system_new=='Vega':
                        dataframe[filt] = dataframe[filt] - mist_mag_vega_ab[idx]
                    elif system_old=='ST' and system_new=='Vega':
                        dataframe[filt] = dataframe[filt] - mist_mag_vega_st[idx]
                    elif system_old=='AB' and system_new=='ST':
                        dataframe[filt] = dataframe[filt] - mist_mag_vega_ab[idx] + mist_mag_vega_st[idx]
                    elif system_old=='ST' and system_new=='AB':
                        dataframe[filt] = dataframe[filt] - mist_mag_vega_st[idx] + mist_mag_vega_ab[idx]
                    else:
                        raise ValueError('Invalid magnitude system conversion: '+system_old+
                            ' -> '+system_new)
        return dataframe
