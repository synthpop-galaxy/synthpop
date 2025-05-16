Changelog
============

Version 1
---------

v1.0.3
^^^^^^
Date: May 16, 2025

* Cleaned up GULLS post-processing module to prep SynthPop catalogs for GULLS microlensing simulator (@acrisp3)
* New post-processing module: Blending, which computes blended mag within selected blend radius for selected filters (@hustonm)
* Fixed bug in Surot extinction module, where 3-d extinction map could not handle distances <10pc and >30kpc. Now applies 0 extinction at distance<10pc and applies extinction value at 30kpc for distance>30kpc. (@hustonm)
* Updated tutorial for new configuration (@hustonm)
* Added ebfpy to requirements, which is used by some extinction map modules (@hustonm)
* Addressed futurewarnings for Pandas 3.0 (@hustonm)

v1.0.2
^^^^^^
Date: April 24, 2025

* Bug fix: Besancon 2003 kinematic module typo in nan handling (catalog output not affected)
* Configuration added: Finalized huston2025_defaults.synthpop_conf model configuration, to be presented by Huston et al. (in prep.)

v1.0.1
^^^^^^
Date: April 16, 2025

* Bug fix: Surot extinction calculations for project_3d=True (project_3d=False unaffected)
* Useability improvement: automated extinction file download and handling for Galaxia-based extinction options

v1.0.0
^^^^^^
Date: April 3, 2025

First official public release with software description paper `(Klüter & Huston et al., 2025) <https://ui.adsabs.harvard.edu/abs/2024arXiv241118821K/abstract>`_.
