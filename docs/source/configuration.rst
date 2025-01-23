Configuration
=======================

.. note::
    Page under construction

SynthPop is configured via json files with a ``.synthpop_conf`` extension. Each run must have a default configuration, and may optionally have a specific configuration which overrides any default values.

The configuration sets the following: 

* Sightlines
* Coordinate system definition
* Population generation settings
* Extinction
* Isochrones + interpolator
* Photometry settings
* Output contents + formatting

Specific and Default Configurations
-----------------------------------
We provide a ``_default_config.synthpop_conf`` file, which is the default file used if a custom one is not entered.
A specific configuration will override any parameters defined in the default.

The default values are defined in config_files/_default_config.json.
Within the specific config file, it is sufficient to specify only the items 
which differ from the default config file.  
We recommend to include all information you think are 
usefull to understand the meaning of your generated catalog.
For example the config file can look as follows::
    
    {
      "model_base_name":"my_generated_model",
      "l_set": [0, 2],
      "b_set": [-1],
      "solid_angle":5e-7,
      "model_name":"besancon_Robin2003",
      "evolution_class":{
        "name":"MIST", 
        "interpolator":"CharonInterpolator"
        },
      "extinction_map_kwargs": {"name":"Marshall"}
      "extinction_law_kwargs": {"name":"ODonnellCardelli"}
    }

Setting Your Configuration
---------------------------
To set a configuration for a SynthPop run via command line, use::

  python -m synthpop <config_file>

or::

  python -m synthpop --specific_config <config_file> --default_config <default_config_file>

The equivalent process in a script that imports synthpop is::

  import synthpop
  model = synthpop.SynthPop('config_file',**kwargs)
  model.init_population()
  model.process_all()

where ``default_config`` is an optional keyword argument.

Configuration File Contents
-----------------------------

Our ``_default_config.json`` is sorted into different categories, which are not required but may be used for any configuration. Note that any arguments starting with an '#' are ignored, so we use these for comments. The optional category headings are: MANDITORY, SIGHTLINES, SEED, COORDINATE_SYSTEM, POPULATION_GENERATION, EXTINCTION_MAP, ISOCHRONE_INTERPOLATION, PHOTOMETRIC_OUTPUTS, and OUTPUT.

MANDITORY
^^^^^^^^^
**model_name**: (string) model directory containing the population files

**name_for_output**: (string) default directory and base for the output files

example::

    "MANDATORY":{
        "model_name":"Model1",
        "name_for_output":"Model1_v1"
    },

SIGHTLINES
^^^^^^^^^^
**l_set**, **b_set**: lists of Galactic longitude and latitudes for sightlines to generate (degrees)

**l_set_type**, **b_set_type**: three options for how to handle **l_set** and **b_set**:

* "list": create a grid of sightlines, with all l and b values looped over individually, i.e. ((l[0],b[0]),...,(l[0],b[i]),(l[1],b[0]),... 
* "pairs": treats the l and b sets as a list of pairs, i.e. ((l[0],b[0]),..., (l[i],b[i])
* "range": treats the l and b sets as arguments for np.arange([l/b]_set)

**solid_angle**: solid angle of cone FoV in units defined by **solid_angle_unit**

**solid_angle_unit**: angle for solid angle of FoV (options: "deg^2" for square degree, "sr" for steradian)

example::

    "SIGHTLINES":{
            "l_set": [3,4],
            "l_set_type":"list",
            "b_set":[-1,0,1],
            "b_set_type":"list",
            "solid_angle": 1e-2,
            "solid_angle_unit": "deg^2"
    }

SEED
^^^^
**random_seed**

COORDINATE_SYSTEMS
^^^^^^^^^^^^^^^^^^
**sun**: dictionary containing the following values:

* **x**, **y**, **z**: location of the Sun in cartesian Galactic coordinates (kpc)
* **u**, **v**, **w**: motion of the Sun in cartesian Galactic coordinates (km/s)
* **l_apex_deg**, **b_apex_deg**: direction of the Solar apex in Galactic coordinates (degree)

**lsr**: dictionary containing the following values:

* **u_lsr**, **v_lsr**, **w_lsr**: velocity of the local standard of rest (km/s) in cartesian Galactic coordinates

**warp**: dictionary describing the warp of the galaxy, where warp is parametrized as z_w=a(R-R_W)^b sin(phi-phi_W) (see `Chen et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019NatAs...3..320C/abstract>`_) [note: these values can be overwritten in population files]:

* **r_warp**: radius where the warp begins (kpc, R_W in equation below)
* **amp_warp**: warp amplitude (a in equation below) [note: set this to null to use **amp_warp_pos** + **amp_warp_neg**]
* **amp_warp_pos**, **amp_warp_neg**: warp amplitude (a in equation below), for positive and negative phi [note: set these to null to use **amp_warp**]
* **alpha_warp**: exponent of warp power law (b in equation below)
* **phi_warp_deg**, **phi_warp_rad**: angle for line of notes, where one should be specified with indicated unit (degree or radian)

example with solar motion from Reid & Brunthaler (2020), LSR motion from Schönrich et al. (2010), warp from Chen X. et al (2019)::

    "COORDINATE_SYSTEM":{
        "sun": {
            "x": -8.178,
            "y": 0.0,
            "z": 0.017,
            "u": 12.9,
            "v": 245.6,
            "w": 7.78,
            "l_apex_deg": 56.24,
            "b_apex_deg": 22.54
        },
        "lsr":{
            "u_lsr": 1.8,
            "v_lsr": 233.4,
            "w_lsr": 0.53
        },
        "warp": {
            "r_warp": 7.72,
            "amp_warp": 0.060,
            "amp_warp_pos": null,
            "amp_warp_neg": null,
            "alpha_warp": 1.33,
            "phi_warp_deg": 17.5
        }
    },

POPULATION_GENERATION
^^^^^^^^^^^^^^^^^^^^^
**max_distance**: maximum distance for stars in catalog (kpc)

**distance_step_size**: step size for generation of stars as slices in distance (kpc)

**window_type**: dictionary containing the following:

* **window_type**: currently must be set to "cone" [note: other options are planned but not in the immediate future]
* **kwargs**

**mass_lims**: range of initial stellar masses to produce

**N_mc_totmass**: number of random points to sample to estimate average density in a slice

**lost_mass_option**: method to estimate correction for mass loss with four integer options:

* 1: For each population, a test batch of N_av_mass stars is generated and evolved to estimate the total initial stellar mass required to meet the desired present day total stellar mass. These values are saved for all sightlines run with the initialized populations.
* 2: For each population a test batch of N_av_mass stars is generated and evolved to estimate the total initial stellar mass required to meet the desired present day total stellar mass. These value is re-calculated for each sightline run.
* 3: Initially treat the population density as an initial mass density, then add or remove stars as needed.
* 4: Use the precomputed value given by either "av_mass_corr" or "n_star_corr" in each population file to scale the required total initial stellar mass needed to achieve the desired present day total stellar mass.

**N_av_mass**: number of stars to use to estimate average evolved stellar mass

**kinematics_at_the_end**: sets whether to determine stellar masses are evolved at the end of the process, instead of as stars are generated (boolean)

**scale_factor**: scale down number of generated stars as n_generated = (n_total/scale_factor)

**skip_lowmass_stars**: option to skip the generation of low mass stars which cannot be bright enough to reach the magnitude cut [note: improves runtimes and memory usage]

**chunk_size**: for computational feasibility, limit number of stars to evolve at once to this value

example::

    "POPULATION_GENERATION":{
        "max_distance":25,
        "distance_step_size":0.10,
        "window_type":{"window_type":"cone"},
        "mass_lims":{"min_mass":0.08,"max_mass":100},
        "N_mc_totmass":10000,
        "lost_mass_option": 1,
        "N_av_mass":20000,
        "kinematics_at_the_end":true,
        "scale_factor": 1,
        "skip_lowmass_stars": false,
        "chunk_size": 250000
    },

EXTINCTION_MAP
^^^^^^^^^^^^^^

**extinction_map_kwargs**: dictionary containing:

* **name**: name of extinction map module
* **<kwargs>**: any kwargs required or optional for the selected module

**extinction_law_kwargs**: dictionary containing:

* **name**: name of extinction law module
* **R_V**: total to selective extinction ratio [note: only used in select extinction laws, will be ignored if input for others]
* **<kwargs>**: any kwargs required or optional for the selected module

example::

    "EXTINCTION_MAP":
        {
        "extinction_map_kwargs":{
            "name":"MapsFromDustmaps", "dustmap_name":"marshall"
            },
        "extinction_law_kwargs":
            [
            {"name":"SODC", "R_V":2.5}
            ]
        },

ISOCHRONE_INTERPOLATION
^^^^^^^^^^^^^^^^^^^^^^^

**evolution_class**: dictionary containing:
* **name**: name of stellar evolution class
* **interpolator**: name of isochrone interpolator class

example for single evolution class::

    "ISOCHRONE_INTERPOLATION":{
        "evolution_class": {"name":"MIST", "interpolator":"CharonInterpolator"},
    }

example for evolution class determined by initial mass, where we iterate through the list and apply the first appropriate option per star::

    "ISOCHRONE_INTERPOLATION":{
        "evolution_class":[
                {"name":"MIST", "min_mass":0.2, "max_mass":0.3},
                {"name":"MIST", "interpolator":"LagrangeInterpolator","min_mass":0.1, "max_mass":0.7},
                {"name":"MIST", "interpolator":"CharonInterpolator"}
            ]
    }

example for evolution class determined by population::

    "ISOCHRONE_INTERPOLATION":{
        "evolution_class":{
            "default":[
                {"name":"MIST", "interpolator":"CharonInterpolator"}
                ],
            "<population_name>":[
                        {"name":"MIST", "interpolator":"LagrangeInterpolator"}
            ]
    }

PHOTOMETRIC_OUTPUTS
^^^^^^^^^^^^^^^^^^^

**mag_lim**: list containing the band to select on, the magnitude limit in that band, and "keep" or "remove" for whether to drop stars dimmer than the limit

**chosen_bands**: list of filters to include for synthetic photometry

For the MIST evolution module, the following filters are available:

.. list-table::
    :widths: 20 80
    :header-rows: 1

    * - System
      - Filters (use the names as written here)
    * - CFHT
      - CFHT_u, CFHT_CaHK, CFHT_g, CFHT_r, CFHT_i_new, CFHT_i_old, CFHT_z
    * - DECam
      - DECam_u, DECam_g, DECam_r, DECam_i, DECam_z, DECam_Y
    * - GALEX
      - GALEX_FUV, GALEX_NUV
    * - HST_ACSHR
      - ACS_HRC_F220W, ACS_HRC_F250W, ACS_HRC_F330W, ACS_HRC_F344N, ACS_HRC_F435W, ACS_HRC_F475W, ACS_HRC_F502N, ACS_HRC_F550M, ACS_HRC_F555W, ACS_HRC_F606W, ACS_HRC_F625W, ACS_HRC_F658N, ACS_HRC_F660N, ACS_HRC_F775W, ACS_HRC_F814W, ACS_HRC_F850LP, ACS_HRC_F892N
    * - HST_ACSWF
      - ACS_WFC_F435W, ACS_WFC_F475W, ACS_WFC_F502N, ACS_WFC_F550M, ACS_WFC_F555W, ACS_WFC_F606W, ACS_WFC_F625W, ACS_WFC_F658N, ACS_WFC_F660N, ACS_WFC_F775W, ACS_WFC_F814W, ACS_WFC_F850LP, ACS_WFC_F892N
    * - HST_WFC3
      - WFC3_UVIS_F200LP, WFC3_UVIS_F218W, WFC3_UVIS_F225W, WFC3_UVIS_F275W, WFC3_UVIS_F280N, WFC3_UVIS_F300X, WFC3_UVIS_F336W, WFC3_UVIS_F343N, WFC3_UVIS_F350LP, WFC3_UVIS_F373N, WFC3_UVIS_F390M, WFC3_UVIS_F390W, WFC3_UVIS_F395N, WFC3_UVIS_F410M, WFC3_UVIS_F438W, WFC3_UVIS_F467M, WFC3_UVIS_F469N, WFC3_UVIS_F475W, WFC3_UVIS_F475X, WFC3_UVIS_F487N, WFC3_UVIS_F502N, WFC3_UVIS_F547M, WFC3_UVIS_F555W, WFC3_UVIS_F600LP, WFC3_UVIS_F606W, WFC3_UVIS_F621M, WFC3_UVIS_F625W, WFC3_UVIS_F631N, WFC3_UVIS_F645N, WFC3_UVIS_F656N, WFC3_UVIS_F657N, WFC3_UVIS_F658N, WFC3_UVIS_F665N, WFC3_UVIS_F673N, WFC3_UVIS_F680N, WFC3_UVIS_F689M, WFC3_UVIS_F763M, WFC3_UVIS_F775W, WFC3_UVIS_F814W, WFC3_UVIS_F845M, WFC3_UVIS_F850LP, WFC3_UVIS_F953N, WFC3_IR_F098M, WFC3_IR_F105W, WFC3_IR_F110W, WFC3_IR_F125W, WFC3_IR_F126N, WFC3_IR_F127M, WFC3_IR_F128N, WFC3_IR_F130N, WFC3_IR_F132N, WFC3_IR_F139M, WFC3_IR_F140W, WFC3_IR_F153M, WFC3_IR_F160W, WFC3_IR_F164N, WFC3_IR_F167N
    * - HST_WFPC2
      - WFPC2_F218W, WFPC2_F255W, WFPC2_F300W, WFPC2_F336W, WFPC2_F439W, WFPC2_F450W, WFPC2_F555W, WFPC2_F606W, WFPC2_F622W, WFPC2_F675W, WFPC2_F791W, WFPC2_F814W, WFPC2_F850LP
    * - IPHAS
      - INT_IPHAS_gR, INT_IPHAS_Ha, INT_IPHAS_gI
    * - JWST
      - F070W, F090W, F115W, F140M, F150W2, F150W, F162M, F164N, F182M, F187N, F200W, F210M, F212N, F250M, F277W, F300M, F322W2, F323N, F335M, F356W, F360M, F405N, F410M, F430M, F444W, F460M, F466N, F470N, F480M
    * - LSST
      - LSST_u, LSST_g, LSST_r, LSST_i, LSST_z, LSST_y
    * - PanSTARRS
      - PS_g, PS_r, PS_i, PS_z, PS_y, PS_w, PS_open
    * - SDSSugriz
      - SDSS_u, SDSS_g, SDSS_r, SDSS_i, SDSS_z
    * - SkyMapper
      - SkyMapper_u, SkyMapper_v, SkyMapper_g, SkyMapper_r, SkyMapper_i, SkyMapper_z
    * - SPITZER
      - IRAC_3.6, IRAC_4.5, IRAC_5.8, IRAC_8.0
    * - HSC
      - hsc_g, hsc_r, hsc_i, hsc_z, hsc_y, hsc_nb816, hsc_nb921
    * - Swift
      - Swift_UVW2, Swift_UVM2, Swift_UVW1, Swift_U, Swift_B, Swift_V
    * - UBVRIplus
      - Bessell_U, Bessell_B, Bessell_V, Bessell_R, Bessell_I, 2MASS_J, 2MASS_H, 2MASS_Ks, Kepler_Kp, Kepler_D51, Hipparcos_Hp, Tycho_B, Tycho_V, Gaia_G_DR2Rev, Gaia_BP_DR2Rev, Gaia_RP_DR2Rev, Gaia_G_MAW, Gaia_BP_MAWb, Gaia_BP_MAWf, Gaia_RP_MAW, TESS, Gaia_G_EDR3, Gaia_BP_EDR3, Gaia_RP_EDR3
    * - UKIDSS
      - UKIDSS_Z, UKIDSS_Y, UKIDSS_J, UKIDSS_H, UKIDSS_K
    * - VISTA
      - VISTA_Z, VISTA_Y, VISTA_J, VISTA_H, VISTA_Ks
    * - WashDD0uvby
      - Washington_C, Washington_M, Washington_T1, Washington_T2, DDO51_vac, DDO51_f31, Stromgren_u, Stromgren_v, Stromgren_b, Stromgren_y
    * - WFIRST
      - R062, Z087, Y106, J129, W146, H158, F184
    * - WISE
      - WISE_W1, WISE_W2, WISE_W3, WISE_W4
    * - SPLUS
      - SPLUS_uJAVA, SPLUS_gSDSS, SPLUS_rSDSS, SPLUS_iSDSS, SPLUS_zSDSS, SPLUS_J0340, SPLUS_J0378, SPLUS_J0395, SPLUS_J0410, SPLUS_J0515, SPLUS_J0660, SPLUS_J0861
    * - UVIT
      - UVIT_F148W, UVIT_F154W, UVIT_F169M, UVIT_F172M, UVIT_N242W, UVIT_N219M, UVIT_N245M, UVIT_N263M, UVIT_N279N

**eff_wavelengths**: dictionary specifying effective wavelength for each chosen filter [Note: use option {"json_file":"AAA_effective_wavelengths.json"} to load these from a pre-existing file]

**obs_mag**: boolean option to generate observed magnitudes (generates absolute magnitudes if set to false)

**opt_iso_props**: optional stellar properties to save, with original column names from isochrones

For MIST isochrone stellar property options, see `their documentation here <https://waps.cfa.harvard.edu/MIST/README_tables.pdf>`_

**col_names**: columns names for output for the columns determined in **opt_iso_props**

example::

    "PHOTOMETRIC_OUTPUTS":{
        "maglim":["Bessell_I", 18, "remove"],
        "chosen_bands": ["R062","Z087","Y106","J129","W146","H158","F184", "Bessell_U", "Bessell_B", "Bessell_V", "Bessell_R", "Bessell_I", "VISTA_J", "VISTA_H", "VISTA_Ks"],
        "eff_wavelengths": {"json_file":"AAA_effective_wavelengths.json"},
        "obsmag":true,

        "opt_iso_props":["log_L", "log_Teff", "log_g", "[Fe/H]","log_R"],
        "col_names":["logL", "Teff", "logg" ,"Fe/H_evolved","log_radius"]
    },

OUTPUT
^^^^^^
**post_processing_kwargs**: null, or list of dictionaries for postprocessing modules (to be executed in order), where the dictionaries contain:

* **name**: name of postprocessing module
* **<kwargs>**: required or optional kwargs for the given extinction module

**output_location**: path for folder to save output to

**output_filename_pattern**: string describing naming system for output files. Accessible values are model_base_name (str), model_name (str), l_deg (float), b_deg(float), solid_angle (float), date (datetime.date object), time (datetime.time object).

**output_file_type**: list containing output file type and dictionary for additional kwargs, saved via pandas or astropy. valid options: csv, json, html, xml, excel, hdf5, feather, parquet, stata, pickle, sql, fits, vot

**overwrite**: boolean option to overwrite existing output files of the same name

example::

    "OUTPUT":{
        "post_processing_kwargs": [{"name":"ProcessDarkCompactObjects", "remove":true}],
        "output_location":"outputfiles/testing",
        "output_filename_pattern": "{name_for_output}_l{l_deg:.3f}_b{b_deg:.3f}",
        "output_file_type": ["csv",{}],
        "overwrite": true
    }
