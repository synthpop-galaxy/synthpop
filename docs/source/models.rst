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
    Because SynthPop uses different code frameworks than other models (e.g. different numerical techniques, approximations, and configuration options), these models will not exactly reproduce the catalogs and results that come from the cited works if the cited work did not use SynthPop itself.

Huston2024
^^^^^^^^^^^^^^^^^^^
Preliminary version presented by Huston et al. (in prep.), with citations therein to model components drawn from many existing models. This version was used by the Roman Galactic Exoplanet Survey Project Infrastructure Team for some preparations for microlensing science with the Roman Galactic Bulge Time Domain Survey.

Huston2025
^^^^^^^^^^^^^^^^^^^
Final model version presented by Huston et al. (in prep.), with citations therein to model components drawn from many existing models. This version was used by the Roman Galactic Exoplanet Survey Project Infrastructure Team for some preparations for microlensing science with the Roman Galactic Bulge Time Domain Survey. There are two key differences from Huston2024: the inclusion of a nuclear stellar disk and 3-d handling of the dust.

besancon_Robin2003
^^^^^^^^^^^^^^^^^^
SynthPop implementation of the Besancon Galactic model, as presented by `Robin et al. (2003) <https://ui.adsabs.harvard.edu/abs/2003A%26A...409..523R/abstract>`__.

GUMS_dr3
^^^^^^^^
SynthPop implementation of the Gaia Universe Model Snapshot (GUMS), DR3 version, as described in the `Gaia DR3 Documenation, Section 2.2 <https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu2sim/sec_cu2UM/>`__.

GUMS_dr3_mod_dens
^^^^^^^^^^^^^^^^^
SynthPop implementation of the Gaia Universe Model Snapshot (GUMS), DR3 version, as described in the `Gaia DR3 Documenation, Section 2.2 <https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu2sim/sec_cu2UM/>`__. This version is rescaled to better match their catalog with the SynthPop code implementation (see `Klüter & Huston et al. (2025) <https://ui.adsabs.harvard.edu/abs/2024arXiv241118821K/abstract>`__).

Koshimoto2022
^^^^^^^^^^^^^
Galactic model presented by `Koshimoto et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract>`__ and `Koshimoto et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022zndo...6941086K/abstract>`__ (`genstars <https://github.com/nkoshimoto/genstars>`__), which provides multiple model versions. Here, we provide the E+Ex bulge model, the NEED TO UPDATE DISK VERSION, and their nuclear stellar disk model, which is based on `Sormani et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.1857S/abstract>`__.

Additional Populations
----------------
.. note:: We provide some additional population files under **models/spare_populations**, which should not be taken as complete models but may be used for population studies or in combination with other model components.

Cao2013_bulge
^^^^^^^^^^^^^
Galactic bulge density model based on OGLE-III red clump giants from `Cao et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013MNRAS.434..595C/abstract>`__, with additional parameters from other models.

Koshimoto2022_E_bulge_nsd
^^^^^^^^^^^^^^^^^^^^^^^^^
Alternate E bulge model from `Koshimoto et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract>`__, with corresponding NSD implementation from `Sormani et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.1857S/abstract>`__ with the IMF from `Koshimoto et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract>`__ for this bulge parameterization and the age and metallicity used by `genstars <https://github.com/nkoshimoto/genstars>`__.

Koshimoto2022_G_bulge_nsd
^^^^^^^^^^^^^^^^^^^^^^^^^
Alternate G bulge model from `Koshimoto et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract>`__, with corresponding NSD implementation from `Sormani et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.1857S/abstract>`__ with the IMF from `Koshimoto et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract>`__ for this bulge parameterization and the age and metallicity used by `genstars <https://github.com/nkoshimoto/genstars>`__.

Koshimoto2022_GxG_bulge_nsd
^^^^^^^^^^^^^^^^^^^^^^^^^
Alternate G+Gx bulge model from `Koshimoto et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract>`__, with corresponding NSD implementation from `Sormani et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.1857S/abstract>`__ with the IMF from `Koshimoto et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract>`__ for this bulge parameterization and the age and metallicity used by `genstars <https://github.com/nkoshimoto/genstars>`__.

Sormani2022_nsd
^^^^^^^^^^^^^^^
Nuclear Stellar Disk population from `Sormani et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.1857S/abstract>`__.
