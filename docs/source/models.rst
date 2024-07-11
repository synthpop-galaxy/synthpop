Models
======

Model & Population Structure
----------------------------

In the code, each Model option exists as a directory in ``synthpop/modules``.
A model is defined by a collection of populations. 
Each population is further described by the following 5 modules:

1. Population Density
2. Initial Mass Function
3. Age Distribution
4. Metallicity Distribution
5. Kinematics

To avoid ambiguity, a ".popjson" extension should be used for the population files. 
Each files must define the following keywords::

    "name" : "name_of_the_population"
    "imf_func_kwargs" : ...
    "age_func_kwargs" : ...
    "metallicity_func_kwargs" : ...
    "kinematics_func_kwargs" : ...
    "population_density_kwargs" : ...

Each of the kwargs items includes a sub-dictionary 
specifying the module name and any keyword arguments for the module, e.g.::

    "metallicity_func_kwargs" : {
        "name" : "double_gaussian",
        "weight" : 0.323,
        "mean1" : -0.31,
        "std1" : 0.31,
        "mean2" : 0.26,
        "std2" : 0.20
    }

Additional kwargs can be used to define galactic warp (see Configuration info for parametrization), e.g.::

    "warp":{
        "r_warp": 9.8065,
        "amp_warp_pos": 0.61675,
        "amp_warp_neg": 0.18500,
        "phi_warp_rad": -0.1570536732,
        "alpha_warp": 1.0
    }

Available Models
----------------
.. note::
    additional models + descriptions coming soon

synthpop_Huston2024
^^^^^^^^^^^^^^^^^^^

besancon_Robin2003
^^^^^^^^^^^^^^^^^^

GUMS_dr3
^^^^^^^^

GUMS_dr3_mod_dens
^^^^^^^^^^^^^^^^^

Koshimoto2021
^^^^^^^^^^^^^

