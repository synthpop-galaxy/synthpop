Changelog
============

Version 2
---------

v2.0.0
^^^^^^
Date: TBD

See the Configuration docs page for details on how to implement new functionality.

Major Changes:

* Initial-Final Mass Relation (IFMR): Stellar remnants that have evolved beyond the isochrone grid are now handled in the main body rather than post-processing. They can just be excluded with the `None` option. Either way, their absence or true final masses are properly accounted for in the total stellar mass.
* Multiplicity: The Multiplicity module can be used to generate stellar companions. The output of a SynthPop run is now two catalogs (with distributions no longer being returned), the systems table and the companions table. The companion table will be returned as `None` when multiplicity is turned off. Configuration options allow the systems table to provide either primary star photometry or combined system photometry, while the companions table always has the photometry for each individual companion star.
* SPISEA Integration: The SPISEA stellar populations code can now be used as the basis for stellar evolution in SynthPop. This requires specific SPISEA-based modules to be used for Evolution, IFMR, and Multiplicity. Optionally, SPISEA IMFs may be used, but SynthPop's built-in piecewise power laws can be used as-is. SPISEA exinction laws are also available.
* Improved density handling: We have removed the 'uniform within slice' simplification in SynthPop v1 (see Klueter & Huston et al. 2025). Total stellar mass integration and stellar position drawing is now performed via trapezoidal integration across a 3d grid with user-modifiable resolution.

Other new functionality:

* SynthPop can now generate square or rectangular fields, in addition to circular.
* Photometric systems (i.e. Vega, AB, ST) are now handled in the main body of the code rather than post-processing.
* The code will validate the extinction model before generating a catalog to ensure the extinction law includes the selected photometric filter wavelengths and that the extinction map returns valid values around the outer border of the requested field.

Version 1
---------

v1.1.3
^^^^^^
Date: May 21, 2026

* Fixed a bug in our implementation of the Koshimoto et al. (2021) bulge model, where a unit mismatch cause the streaming motion along the bar to not get applied.

v1.1.2
^^^^^^
Date: May 11, 2026

* New postprocessing module `EstimateRomanExtinction` uses the results of a polynomial fit of spectrum-level simulated per-filter extinctions from A_Ks and colors to more accurately estimate extinction in the Roman filters. Details will be presented in an upcoming paper from Huston et al.
* New postprocessing module `DropColumns` to remove unnecessary columns from the final output tables.
* A few very minor efficiency and logging improvements.

v1.1.1
^^^^^^
Date: October 6, 2025

* Koshimoto bulge kinematic model: bug fix in handling of solid body motion versus triaxial streaming motion. This primarily impacts radial velocities and has a small impact on proper motion in the Galactic l-direction as one goes to higher |l|. One may use the RecalculateKinematics module to resolve this issue in existing catalogs.(@hustonm & @acrisp3)
* PopSyCLE post-processing: removal of hard-coded column number for extinction handling. (@hustonm)

v1.1.0
^^^^^^
Date: August 20, 2025

User interface and requirement changes:

* The configuration keyword `col_names` has been removed in order to avoid issues with inconsistent naming in post-processing. Instead, a user can change their non-magnitude isochrone property column names via the `RenameColumns` post-processing module. If multiple post-processing modules are given, they will be executed in the order they are provided in the configuration. Thus, we recommend using `RenameColumns` last so that no issues arise in other post-processing modules. Default configurations `huston2025_defaults.synthpop_conf` and `_default.synthpop_conf` have been updated accordingly.
* We have replaced the use of `ebfpy` with `h5py`, changing the code requirements.

Post-processing module updates:

* The `ProcessDarkCompactObjects` module in previous versions contained a bug for the `SukhboldN20` and `Raithel18` IFMR options where some neutron stars were incorrectly assigned to be white dwarfs. This has been resolved.
* The `ProcessDarkCompactObjects` module now includes the option to add natal kicks, which is off by default. A user may provide a mean kick velocity for neutron stars and one for black holes. Kick directions are drawn randomly from a uniform sphere, and kick magnitudes are drawn from a Maxwellian distribution around the input mean kick velocities in km/s.
* The `ProcessDarkCompactObjects` module now notes compact object type in the `phase` output column, using 101 for white dwarfs, 102 for neutron stars, and 103 for black holes. The `Dim_Compact_Object_Flag` column has been eliminated.
* The `GullsPostProcessing` module has been modified to be consistent with the changes to `col_names` described above.
* A new `PopsyclePostProcessing` module has been added to output a star catalog in a format compatible with the PopSyCLE microlensing survey simulator. Details in the module explain configuration requirements for compatibility.

Extinction module updates:

* The `Galaxia_3D` module has been renamed to `Galaxia_3d` for simplicity in SynthPop's module loader. 
* The extinction map files required for use in the `Galaxia_3d` and `Surot` extinction modules are now stored via LSU Box in .h5 format for ease of automated download and loading and to eliminate the use of .ebf files, which raised issues for some `numpy` versions.

Other minor changes:

* The calculation of `vr_lsr` has been moved from `Population` to `CoordinatesTransformation` in order to be adaptable to a modified solar frame.

v1.0.5
^^^^^^
Date: July 18, 2025

Significant change -  debugging of Koshimoto2021Bulge model:

* Coordinate system inconsistency fixed
* Postprocessing module RecalculateKinematics added in order to regenerate correct kinematics for existing catalogs

Any catalogs generated using that module prior to SynthPop v1.0.5 needs to be corrected for kinematics to be reliable. In the direction of the Galactic bulge, the issue primarily effects radial velocities, with a minimal impact on proper motions, as it is primarily related to the U velocity component. 

Example code to run the kinematics model correction is shown below.
You'll need to reload your original model configuration with the new post_processing module added, initialize it, then run each catalog through the post_processing function. Note that the `pop_ids` variable can be set to that of the bulge or excluded to rerun all populations' kinematics.

```
model = sp.SynthPop(default_config='huston2025_defaults.synthpop_conf',
                  post_processing_kwargs={"name":"recalculate_kinematics", "pop_ids":[0]})
model.init_populations() 
for f in file_list:
    df = pandas.read_csv(f)
    df = model.post_processing(df)
    df.to_csv(f)
```

Minor changes:

* New postprocessing module EquatorialCoordinates to convert from l,b to RA,Dec, with the option of keeping both
* Small documentation cleanup

v1.0.4
^^^^^^
Date: June 4, 2025

* Fix edge-case bugs in lost_mass_option=3 star generation process (@hustonm)
* Only collect needed effective wavelengths from json file to de-clutter configurations/logs (@hustonm)
* Change fill value to 99 for missing Mbol values in GULLS post-processing (@acrisp3)

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
