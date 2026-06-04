import synthpop 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

synthdata = synthpop.SynthPop('huston2025_defaults.synthpop_conf',
                         chosen_bands = ['Gaia_G_EDR3', 'Gaia_BP_EDR3', 'Gaia_RP_EDR3'],
                         maglim = ['Gaia_G_EDR3', 20, "remove"],
                         post_processing_kwargs=[{"name": "ProcessDarkCompactObjects","remove": True}],
                         name_for_output='mod2test'
                        )

#=================================================
# Proper motion in l
#=================================================

synthdata.init_populations()

l_vals = np.arange(-10, 11, 1)
b_vals = np.arange(-10, 6, 1)
#bin_edges = np.linspace(-30, 30, 100)

#fig, axes = plt.subplots(len(b_vals), len(l_vals), figsize=(6, 6), sharex=True, sharey=True, squeeze=False)
#fig.suptitle("Comparison of Gaia vs. Synthpop Proper Motions (mul / pm_l*)", fontsize=16)

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
Gaia.ROW_LIMIT = -1

stats_list = []

for i, b in enumerate(b_vals):
    for j, l in enumerate(l_vals):
        stats_entry = {'l': l, 'b': b}
        try:
            # GAIA
            center = SkyCoord(l=l*u.deg, b=b*u.deg, frame='galactic').icrs
            radius = 0.1 * u.deg

            job = Gaia.cone_search_async(center, radius)
            result = job.get_results()

            mask = (result['phot_g_mean_mag'] < 20)
            filtered = result[mask]
            valid_pm = (~np.isnan(filtered['pmra'])) & (~np.isnan(filtered['pmdec']))
            filtered = filtered[valid_pm]

            ra_vals = np.array(filtered['ra'], dtype=float)
            dec_vals = np.array(filtered['dec'], dtype=float)
            pmra_vals = np.array(filtered['pmra'], dtype=float)
            pmdec_vals = np.array(filtered['pmdec'], dtype=float)

            coords = SkyCoord(ra=ra_vals*u.deg,
                              dec=dec_vals*u.deg,
                              pm_ra_cosdec=pmra_vals*u.mas/u.yr,
                              pm_dec=pmdec_vals*u.mas/u.yr,
                              frame='icrs')

            gal_coords = coords.transform_to('galactic')
            pm_l_cosb = gal_coords.pm_l_cosb.to(u.mas/u.yr).value

            if len(pm_l_cosb) > 0:
                stats_entry['gaia_median'] = np.median(pm_l_cosb)
                stats_entry['gaia_std'] = np.std(pm_l_cosb)
            else:
                stats_entry['gaia_median'] = np.nan
                stats_entry['gaia_std'] = np.nan

            # SYNTHPOP
            cat2, distr2 = synthdata.process_location(
                l_deg=l, b_deg=b, solid_angle=0.01*np.pi, solid_angle_unit='deg^2'
            )
            mul_vals = cat2['mul'].astype(float).values

            if len(mul_vals) > 0:
                stats_entry['synth_median'] = np.median(mul_vals)
                stats_entry['synth_std'] = np.std(mul_vals)
            else:
                stats_entry['synth_median'] = np.nan
                stats_entry['synth_std'] = np.nan

        except Exception as e:
            stats_entry['gaia_median'] = np.nan
            stats_entry['gaia_std'] = np.nan
            stats_entry['synth_median'] = np.nan
            stats_entry['synth_std'] = np.nan

        stats_list.append(stats_entry)

stats_df = pd.DataFrame(stats_list)
stats_df = stats_df.sort_values(by=['b', 'l']).reset_index(drop=True)

print(stats_df.head())

stats_df.to_csv("l_gaia_vs_synthpop_new_run.csv", index=False)

#=============================================
# Proper motion in b
#=============================================


stats_list = []

for i, b in enumerate(b_vals):
    for j, l in enumerate(l_vals):
        stats_entry = {'l': l, 'b': b}
        try:
            # GAIA
            center = SkyCoord(l=l*u.deg, b=b*u.deg, frame='galactic').icrs
            radius = 0.1 * u.deg

            job = Gaia.cone_search_async(center, radius)
            result = job.get_results()

            mask = (result['phot_g_mean_mag'] < 20)
            filtered = result[mask]
            valid_pm = (~np.isnan(filtered['pmra'])) & (~np.isnan(filtered['pmdec']))
            filtered = filtered[valid_pm]

            ra_vals = np.array(filtered['ra'], dtype=float)
            dec_vals = np.array(filtered['dec'], dtype=float)
            pmra_vals = np.array(filtered['pmra'], dtype=float)
            pmdec_vals = np.array(filtered['pmdec'], dtype=float)

            coords = SkyCoord(ra=ra_vals*u.deg,
                              dec=dec_vals*u.deg,
                              pm_ra_cosdec=pmra_vals*u.mas/u.yr,
                              pm_dec=pmdec_vals*u.mas/u.yr,
                              frame='icrs')

            gal_coords = coords.transform_to('galactic')
            pm_b = gal_coords.pm_b.to(u.mas/u.yr).value

            if len(pm_b) > 0:
                stats_entry['gaia_median'] = np.median(pm_b)
                stats_entry['gaia_std'] = np.std(pm_b)
            else:
                stats_entry['gaia_median'] = np.nan
                stats_entry['gaia_std'] = np.nan

            # SYNTHPOP
            cat2, distr2 = synthdata.process_location(
                l_deg=l, b_deg=b, solid_angle=0.01*np.pi, solid_angle_unit='deg^2'
            )
            mub_vals = cat2['mub'].astype(float).values

            if len(mub_vals) > 0:
                stats_entry['synth_median'] = np.median(mub_vals)
                stats_entry['synth_std'] = np.std(mub_vals)
            else:
                stats_entry['synth_median'] = np.nan
                stats_entry['synth_std'] = np.nan

        except Exception as e:
            stats_entry['gaia_median'] = np.nan
            stats_entry['gaia_std'] = np.nan
            stats_entry['synth_median'] = np.nan
            stats_entry['synth_std'] = np.nan

        stats_list.append(stats_entry)

stats_df = pd.DataFrame(stats_list)
stats_df = stats_df.sort_values(by=['b', 'l']).reset_index(drop=True)

print(stats_df.head())

stats_df.to_csv("b_gaia_vs_synthpop_new_run.csv", index=False)
