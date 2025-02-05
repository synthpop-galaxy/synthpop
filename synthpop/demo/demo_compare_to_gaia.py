'''
Script to compare SynthPop output to GAIA DR3 Universe Model and Source Catalog
'''

# Imports
import sys
import os
from tqdm import tqdm
import pandas
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from scipy import stats

from astropy.coordinates import SkyCoord, Galactic
import astropy.units as u
from astroquery.gaia import Gaia, GaiaClass
from astroquery.vizier import Vizier
DIRNAME = os.path.dirname(__file__)
if DIRNAME == '': DIRNAME = '.'
sys.path.append(os.path.abspath(os.path.join(DIRNAME,'..')))

import synthpop_main as synthpop
from synthpop_utils import half_cone_angle_to_solidangle
from synthpop_utils.coordinates_transformation import CoordTrans
coordtrans = CoordTrans()
DIRNAME = os.path.dirname(__file__)
if DIRNAME == '': DIRNAME = '.'
# Gaia database query setup
GaiaModel = GaiaClass()
GaiaModel.MAIN_GAIA_TABLE = "gaiadr3.gaia_universe_model"
GaiaModel.ROW_LIMIT = 10000000
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
Gaia.ROW_LIMIT = 10000000
GUMS_QUERY = Vizier(row_limit=-1, catalog="VI/137/gum_mw")

gums_ids= {
    1:(1, 0.075),
    2:(1, 0.575),
    3:(1, 1.5),
    4:(1, 2.5),
    5:(1, 4.0),
    6:(1, 6.0),
    7:(1, 8.5),
    9:(2, 12.0),
    10:(2, 10.0),
    8:(3, 14.0),
    0:(4, 10.0)}

GREEN = np.array([1, 143, 53, 255]) / 255
BLUE = np.array([0, 66, 148, 255]) / 255
RED = np.array([225, 0, 35, 255]) / 255
ORANGE = np.array([225, 190, 35, 255])/255

def write_csv_in_chunks(filename, df, chunksize=1000,columns=None):
    """
    Write a pandas DataFrame to a CSV file in chunks of a specific size.

    Args:
        filename (str): Name of the output file.
        df (pandas.DataFrame): DataFrame to write to the output file.
        chunksize (int, optional): Size of the chunks to write. Default is 1000.

    Returns:
        None
    """

    if not os.path.isdir(dirname := os.path.dirname(filename)):
        os.mkdir(dirname)
    with tqdm(total=len(df)) as pbar, open(filename, 'w') as f:
        header = True
        for _, chunk in df.groupby(df.index // chunksize):
            chunk.to_csv(f, header=header, index=False)
            pbar.update(len(chunk))
            header = False

class CompareGaia:
    '''
    Comparison between Gaia data/model and SynthPop model
    '''
    def __init__(
            self, n_locations=10, radius_deg=(1 / np.sqrt(np.pi) / 10 ** 0.5), loc='disk',
            **kwargs
            ):
        np.random.seed(14)
        
        # Set up sightlines
        self.loc = loc
        if loc.lower() == 'disk':
            self.files = [
                f'{DIRNAME}/comp/gaia_synthpop_compare_synth_gums_coor.h5',
                f'{DIRNAME}/comp/gaia_synthpop_compare_synth_gums.h5',
                f'{DIRNAME}/comp/gaia_synthpop_compare_gaia.csv',
                f'{DIRNAME}/comp/gaia_synthpop_compare_gums.csv'
                ]
            self.kwargs = {}
            #self.l_deg = np.random.uniform(20, 90, n_locations).round(6) % 360
            #min_b = 15 # * (self.l_deg < 15)
            #self.b_deg = np.random.uniform(min_b, 45, n_locations).round(6)
            coord_list = np.transpose([[55.976034, 39.194441],[74.121554, 25.267639],[80.929938, 31.166665],[20.563286, 15.176214],[41.681515, 35.194574],[87.032262, 21.300728],[55.91817, 42.976728],[42.27991, 26.227342],[57.743996, 37.572568],[35.487846, 37.89417]])
            self.l_deg =coord_list[0]
            self.b_deg =coord_list[1]
            self.solid_angle = half_cone_angle_to_solidangle(
                radius_deg / 180 * np.pi) * (180 / np.pi) ** 2
            self.radius_deg = radius_deg * u.degree

        elif loc.lower() == 'bulge':
            self.files = [
                f'{DIRNAME}/comp/gaia_synthpop_compare_synth_gums_coor_bulge.h5',
                f'{DIRNAME}/comp/gaia_synthpop_compare_synth_gums_bulge.h5',
                f'{DIRNAME}/comp/gaia_synthpop_compare_gaia_bulge.csv',
                f'{DIRNAME}/comp/gaia_synthpop_compare_gums_bulge.csv'
                ]
            self.kwargs = {}
            '''self.l_deg = np.zeros(n_locations)
            self.b_deg = np.zeros(n_locations)
            while any((self.l_deg **2+self.b_deg **2)<1):
                self.l_deg = np.random.uniform(0, 15, n_locations).round(6)
                self.b_deg = np.random.uniform(-15, 15, n_locations).round(6)'''
            coord_list = np.transpose([[7.70915, 9.194441],[11.597476, -4.732361],[13.056415, 1.166665],[0.120704, -14.823786],[4.646039, 5.194574],[14.364056, -8.699272],[7.696751, 12.976728],[4.774266, -3.772658],[8.087999, 7.572568],[3.318824, 7.89417]])
            self.l_deg = coord_list[0]
            self.b_deg = coord_list[1]
            print(loc)
            self.radius_deg = radius_deg * u.degree
            self.solid_angle = \
                half_cone_angle_to_solidangle(radius_deg / 180 * np.pi) * (180 / np.pi) ** 2

        # SynthPop model objects, initialization
        self.mod1 = synthpop.SynthPop(
            os.path.join(DIRNAME, 'gaia_compare.synthpop_conf'),
            model_name='GUMS_dr3',
            **self.kwargs)
        self.mod2 = synthpop.SynthPop(
            os.path.join(DIRNAME, 'gaia_compare.synthpop_conf'),
            model_name='GUMS_dr3_mod_dens',
            **self.kwargs)
        self.mod1.init_populations()
        self.mod2.init_populations()

        # Generate new catalogs or read existing
        (self.complete_synth1, self.complete_synth2, self.complete_gaia, self.complete_gums
            ) = self.load_data(**kwargs)

        # Filter out stars with no properties
        self.complete_synth1_filtered = self.complete_synth1.loc[self.complete_synth1["logTeff"] >= 0]
        self.complete_synth2_filtered = self.complete_synth2.loc[self.complete_synth2["logTeff"] >= 0]

    # Download data from Gaia
    def download_cat(self):
        print(self.radius_deg)
        complete_gaia_list = []
        complete_gums_list = []
        for i, loc in enumerate(zip(self.l_deg, self.b_deg, )):
            print(i, loc)
            print('Downloading Gaia')
            coord = SkyCoord(l=loc[0] * u.degree, b=loc[1] * u.degree, frame=Galactic)
            job = Gaia.cone_search_async(coord, self.radius_deg)
            gaia_tab = job.get_results()

            print('Downloading Gaia Universe Model', flush=True)
            job = GaiaModel.cone_search_async(coord, self.radius_deg)
            gums_tab = job.get_results()

            complete_gums_list.append(gums_tab.to_pandas())
            complete_gaia_list.append(gaia_tab.to_pandas())

        complete_gaia = pd.concat(complete_gaia_list)
        print(len(complete_gaia))
        complete_gaia_stars = complete_gaia[complete_gaia["classprob_dsc_combmod_star"] >= 0.6]
        print(len(complete_gaia_stars))
        complete_gums = pd.concat(complete_gums_list)

        print('save gaia', flush=True)
        #print(complete_gums)
        write_csv_in_chunks(self.files[2], complete_gaia_stars, columns=['phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','bp_rp'])
        #complete_gaia_stars.to_hdf(self.files[2],key='data',dropna=True, data_columns=['phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','bp_rp'])

        print('save gums', flush=True)
        #print(complete_gums)
        write_csv_in_chunks(self.files[3], complete_gaia_stars, columns=['mag_g','mag_bp','mag_rp','mass','barycentric_distance','population','age'])
        #complete_gums.to_hdf(self.files[3],key='data', data_columns=['mag_g','mag_bp','mag_rp','mass','barycentric_distance','population','age'])

        return complete_gaia_stars, complete_gums

    # Generate SynthPop catalogs
    def run_synthpop(self):
        """ generate the GUMS model"""
        complete_synth_list1 = []
        complete_synth_list2 = []
        for i, loc in enumerate(zip(self.l_deg, self.b_deg, )):
            print(i, loc, flush=True)
            synth_df2, _ = self.mod2.process_location(*loc, self.solid_angle)
            complete_synth_list2.append(synth_df2)
            synth_df1, _ = self.mod1.process_location(*loc, self.solid_angle)
            complete_synth_list1.append(synth_df1)

        complete_synth1 = pd.concat(complete_synth_list1)
        complete_synth2 = pd.concat(complete_synth_list2)

        complete_synth1['ra'], complete_synth1['dec'] = coordtrans.lb_to_ad(
            complete_synth1.l.values, complete_synth1.b.values)
        complete_synth2['ra'], complete_synth2['dec'] = coordtrans.lb_to_ad(
            complete_synth2.l.values, complete_synth2.b.values)

        _, complete_synth1['pmra'], complete_synth1['pmdec'] = coordtrans.uvw_to_vrmuad(
            complete_synth1.l.values, complete_synth1.b.values, complete_synth1.Dist.values,
            complete_synth1.U.values, complete_synth1.V.values, complete_synth1.W.values)

        _, complete_synth2['pmra'], complete_synth2['pmdec'] = coordtrans.uvw_to_vrmuad(
            complete_synth2.l.values, complete_synth2.b.values, complete_synth2.Dist.values,
            complete_synth2.U.values, complete_synth2.V.values, complete_synth2.W.values)

        print('save model', flush=True)
        #write_csv_in_chunks(self.files[0], complete_synth1)
        complete_synth1.to_hdf(self.files[0],key='data')
        print()
        print('save modified model', flush=True)
        #write_csv_in_chunks(self.files[1], complete_synth2)
        complete_synth2.to_hdf(self.files[1],key='data')

        return complete_synth1, complete_synth2

    # Load the data if it exists, or regenerate
    def load_data(self, regenerate=False, redownload=True):
        print("Load Data")
        data = []
        
        # check if generated data exist
        if os.path.isfile(self.files[0]) and os.path.isfile(self.files[1]) and (not regenerate):
            print(f"load {self.files[0]}", end='', flush=True)
            data.append(pandas.read_hdf(self.files[0], key='data'))
            print("\rload GUMS model    ", end='', flush=True)
            data.append(pandas.read_hdf(self.files[1], key='data'))
            # use downloaded data (if exist)
            print("\rSynthpop data loaded", flush=True)
            redownload = False
        else:
            print("generate Data")
            data.extend(self.run_synthpop())

        # check if downloaded data exist
        if os.path.isfile(self.files[2]) and os.path.isfile(self.files[3]) and (not redownload):
            print("load Catalogs")
            data.append(pandas.read_csv(self.files[2]))
            data.append(pandas.read_csv(self.files[3]))
        else:
            data.extend(self.download_cat())

        return data

    # Histogram, Gaia G band
    def plot_hist_G(self):
        colors = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
        axes[0].set_title(self.loc)
        bins = np.arange(10, 21, 0.25)
        hist_synthpop1, _, _ = axes[0].hist(self.complete_synth1['Gaia_G_EDR3'],
            bins=bins, histtype='step', label='SynthPop', color=BLUE)

        c = next(colors)
        hist_synthpop2, _, _ = axes[0].hist(self.complete_synth2['Gaia_G_EDR3'],
            bins=bins, histtype='step', label='SynthPop (re-normalized)', color=GREEN)
        c = next(colors)

        hist_gaia, _, _ = axes[0].hist(self.complete_gaia['phot_g_mean_mag'],
            bins=bins, histtype='step', label='Gaia Stars', color=ORANGE)
        axes[1].stairs(hist_gaia/hist_synthpop1, bins, color=ORANGE)
        axes[1].stairs(hist_gaia/hist_synthpop2, bins,ls=':', color=ORANGE)

        c = next(colors)
        hist_gums, _, _ = axes[0].hist(self.complete_gums['mag_g'],
            bins=bins, histtype='step', label='GUMS Stars', color=RED)
        axes[1].stairs(hist_gums/hist_synthpop1, bins, color=RED)
        axes[1].stairs(hist_gums/hist_synthpop2, bins, ls=':',color=RED)
        axes[1].set_ylim([0.25, 2])
        axes[0].set_yscale('log')
        axes[0].set_ylabel("#")
        axes[0].legend(loc=2)
        axes[1].set_xlabel("G [mag]")
        axes[1].set_ylabel("$N_{x}/N_{SynthPop}$")
        axes[1].axhline(1, ls='--')
        plt.savefig(f"images/{self.loc}_G_hist.pdf")
        
    # Color-magnitude diagram
    def plot_cmd(self):
        which2 = self.complete_synth2['Gaia_G_EDR3'] < 21
        which_gaia = self.complete_gaia['phot_g_mean_mag'] < 21
        which_gums = self.complete_gums['mag_g'] < 21

        plt.figure()
        plt.plot(
            self.complete_synth2.loc[which2, 'Gaia_BP_EDR3']
            - self.complete_synth2.loc[which2, 'Gaia_RP_EDR3'],
            self.complete_synth2.loc[which2, 'Gaia_G_EDR3'],
            '.', ms=0.1, alpha=0.1, label='SynthPop')
        plt.plot(
            self.complete_gaia.loc[which_gaia, 'phot_bp_mean_mag']
            - self.complete_gaia.loc[which_gaia, 'phot_rp_mean_mag'],
            self.complete_gaia.loc[which_gaia, 'phot_g_mean_mag'],
            '.', ms=0.1, alpha=0.1, label='Gaia Stars')
        plt.plot(
            self.complete_gums.loc[which_gums, 'mag_bp']
            - self.complete_gums.loc[which_gums, 'mag_rp'],
            self.complete_gums.loc[which_gums, 'mag_g'],
            '.', ms=0.1, alpha=0.1, label='GUMS Stars')

        plt.legend()
        plt.gca().invert_yaxis()
        plt.xlabel("BP - RP")
        plt.ylabel("G")
        plt.savefig(f"images/{self.loc}_cmds", dpi=300)
        
    # Color-magnitude diagram comparison
    def plot_cmd_diff(self):
        xmin = 0
        xmax = 6
        ymin = 14
        ymax = 21
        xlim = [xmin, xmax]
        ylim = [ymin, ymax]
        which2 = ((self.complete_synth2['Gaia_G_EDR3'] < ymax + 1)
                  & (self.complete_synth2['Gaia_G_EDR3'] > ymin - 1)
                  & ((self.complete_synth2['Gaia_BP_EDR3']
                      - self.complete_synth2['Gaia_RP_EDR3']) > xmin - 1)
                  & ((self.complete_synth2['Gaia_BP_EDR3']
                      - self.complete_synth2['Gaia_RP_EDR3']) < xmax + 1))

        which_gaia = ((self.complete_gaia['phot_g_mean_mag'] < ymax + 1)
                      & (self.complete_gaia['phot_g_mean_mag'] > ymin - 1)
                      & (self.complete_gaia['bp_rp'] < xmax + 1)
                      & (self.complete_gaia['bp_rp'] > xmin - 1))

        which_gums = ((self.complete_gums['mag_g'] < ymax + 1)
                      & (self.complete_gums['mag_g'] > ymin - 1)
                      & ((self.complete_gums['mag_bp']
                          - self.complete_gums['mag_rp']) > xmin - 1)
                      & ((self.complete_gums['mag_bp']
                          - self.complete_gums['mag_rp']) < xmax + 1))

        X, Y = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        values_synthpop = np.vstack(
            [self.complete_synth2.loc[which2, 'Gaia_BP_EDR3']
             - self.complete_synth2.loc[which2, 'Gaia_RP_EDR3'],
                self.complete_synth2.loc[which2, 'Gaia_G_EDR3']]
            )
        values_gaia = np.vstack(
            [self.complete_gaia.loc[which_gaia, 'bp_rp'],
                self.complete_gaia.loc[which_gaia, 'phot_g_mean_mag']]
            )
        values_gums = np.vstack(
            [self.complete_gums.loc[which_gums, 'mag_bp']
             - self.complete_gums.loc[which_gums, 'mag_rp'],
                self.complete_gums.loc[which_gums, 'mag_g']]
            )
        print("crate kernels")
        kernel_synthpop = stats.gaussian_kde(values_synthpop)
        kernel_gums = stats.gaussian_kde(values_gums)
        kernel_gaia = stats.gaussian_kde(values_gaia)
        print("evaluate kernels")
        Z_synthpop = np.reshape(kernel_synthpop(positions).T, X.shape)
        Z_synthpop *= sum(which2)
        Z_gums = np.reshape(kernel_gums(positions).T, X.shape)
        Z_gums *= sum(which_gums)
        Z_gaia = np.reshape(kernel_gaia(positions).T, X.shape)
        Z_gaia *= sum(which_gaia)

        # fig1, ax1, im1 = make_kde_plot((Z_gums)/Z_synthpop, ,
        #    cmap='cool', vmin=-2, vmax=2)
        # fig2, ax2, im2  = make_kde_plot((Z_gaia)/Z_synthpop, [xmin,xmax], [ymin,ymax],
        #    cmap='bwr', vmin=-2, vmax=2)
        fig3, ax3, im3 = plot_lin_diff((Z_gums - Z_synthpop), xlim, ylim, vmin=-18000,
            vmax=6000, cmap='seismic_r')
        ax3.set_title("Gums - SynthPop")

        fig4, ax4, im4 = plot_lin_diff((Z_gaia - Z_synthpop), xlim, ylim, vmin=-18000,
            vmax=6000, cmap='seismic_r')
        ax4.set_title("Gaia - SynthPop")

        fig5, ax5, im5 = plot_lin_diff((Z_gums - Z_gaia), xlim, ylim, vmin=-18000,
            vmax=6000, cmap='seismic_r')
        ax5.set_title("GUMS - Gaia")
        return Z_synthpop, Z_gums, Z_gaia

    # Plot stellar density by distance
    def plot_density(self, maglim=20.5, bin_width=0.15):
        plt.figure('plot_density')
        which2 = self.complete_synth2['Gaia_G_EDR3'] < maglim
        which_gums = self.complete_gums['mag_g'] < maglim
        mass_2 = self.complete_synth2.loc[which2].Mass
        dist_2 = self.complete_synth2.loc[which2].Dist
        solid_angle = half_cone_angle_to_solidangle(self.radius_deg).value * len(self.l_deg)

        mass_gums = self.complete_gums.loc[which_gums].mass
        dist_gums = self.complete_gums.loc[which_gums].barycentric_distance/1000

        _plot_density(mass_2, dist_2, solid_angle, bin_width,
                label=f'Synthpop mag<{maglim}')
        _plot_density(mass_gums, dist_gums, solid_angle, bin_width,
                label=f'Gaia Universe Model  mag<{maglim}')
        plt.ylabel('Density [Msun/kpc^3]')
        plt.xlabel('Dist [kpc]')
        plt.yscale('log')
        plt.legend()

    # Plot stellar density by distance
    def plot_density_hist(self):
        plt.figure('plot_density_hist')
        d = 0
        gums_group = self.complete_gums.groupby(['population', 'age'])

        solid_angle = half_cone_angle_to_solidangle(self.radius_deg).value
        d_kpc = np.linspace(0, 20, 20000)
        bins = np.arange(0, 20, 0.5)
        data = comp_gaia.complete_synth2
        for i, pop in enumerate(self.mod2.populations):
            plt.figure(pop.name)
            color = iter(cm.rainbow(np.linspace(0, 1, len(self.l_deg))))
            data_c = data[data["pop"] == i]
            current_gums = gums_group.get_group(gums_ids[i])
            for loc in zip(self.l_deg, self.b_deg):
                ra_dec = SkyCoord(l=loc[0] * u.degree, b=loc[1] * u.degree, frame=Galactic).icrs
                alpha, delta = ra_dec.ra.degree, ra_dec.dec.degree
                c = next(color)
                rphiz = pop.coord_trans.dlb_to_rphiz(d_kpc, *loc)
                mask = np.sqrt((data_c.l - loc[0]) ** 2 * np.cos(loc[1] * np.pi / 180) ** 2 + (
                            data_c.b - loc[1]) ** 2) <= self.radius_deg.value
                mask_g = np.sqrt(
                    (current_gums.ra - alpha) ** 2 * np.cos(delta * np.pi / 180) ** 2 + (
                                current_gums.dec - delta) ** 2) <= self.radius_deg.value
                if pop.population_density.density_unit == 'number':
                    weights = np.ones(sum(mask))
                    weights_g = np.ones(sum(mask_g))
                if pop.population_density.density_unit == 'mass':
                    weights = data_c[mask].Mass
                    weights_g = current_gums[mask_g].mass
                if pop.population_density.density == 'init_mass':
                    weights = data_c[mask].iMass
                    weights_g = current_gums[mask_g].mass

                plt.plot(d_kpc, pop.population_density.density(*rphiz), c=c, label = loc)
                _dense_hist(weights, data_c[mask].Dist, bins, solid_angle, color=c)
                _dense_hist(weights_g, current_gums[mask_g].barycentric_distance / 1000, bins,
                    solid_angle, color=c, ls=':')
                plt.yscale('log')
        plt.ylabel('Density [Msun/kpc^3]')
        plt.xlabel('Dist [kpc]')

    # Get model density from SynthPop modules, summed over populations
    def get_combined_density(self, pop_ids, d_kpc):
        dens = np.zeros(d_kpc.shape)
        for loc in zip(self.l_deg, self.b_deg):
            for pop_id in pop_ids:
                population = self.mod1.populations[pop_id]
                rphiz = population.coord_trans.dlb_to_rphiz(d_kpc, *loc)
                dens += population.population_density.density(*rphiz)/len(self.b_deg)
        return dens

    # Plot comparison color-magnitude diagrams
    def plot_comp_cmds(self, alpha=0.1, ms=0.1,**kwargs):
        gums = self.complete_gums
        synth = self.complete_synth2
        gaia = self.complete_gaia

        synth_mag = synth.Gaia_G_EDR3 < 21
        gums_mag = gums.mag_g < 21
        gaia_mag = gaia.phot_g_mean_mag < 21
        gums['abp'] = 1.10 * gums.ag / 0.864
        gums['arp'] = 0.629 * gums.ag / 0.864

        if 'E(B-V)' in synth.columns:
            synth['ag'] = 0.864 * synth['E(B-V)'] * 3.1
            synth['abp'] = 1.10 * synth['E(B-V)'] * 3.1
            synth['arp'] = 0.629 * synth['E(B-V)'] * 3.1
        if 'A0' in synth.columns:
            synth['ag'] = 0.864 * synth['A0']
            synth['abp'] = 1.10 * synth['A0']
            synth['arp'] = 0.629 * synth['A0']
        if 'A_Ks' in synth.columns:
            synth['ag'] = 0.864 * synth.A_Ks / 0.125
            synth['abp'] = 1.10 * synth.A_Ks / 0.125
            synth['arp'] = 0.629 * synth.A_Ks / 0.125

        fig, ax = plt.subplots(2, 3, num='comp_cmd', figsize=(12, 8), sharex='row', sharey='row')

        ax[0, 0].plot(synth[synth_mag].Gaia_BP_EDR3 - synth[synth_mag].Gaia_RP_EDR3,
            synth[synth_mag].Gaia_G_EDR3, '.',color=BLUE, ms=ms, alpha=alpha, **kwargs)
        ax[1, 0].plot(
            synth[synth_mag].Gaia_BP_EDR3 - synth[synth_mag].abp - synth[synth_mag].Gaia_RP_EDR3 +
            synth[synth_mag].arp,
            synth[synth_mag].Gaia_G_EDR3 - synth[synth_mag].ag - 5 * np.log10(
                synth[synth_mag].Dist * 100), '.',color=BLUE, ms=ms, alpha=alpha, **kwargs)

        ax[0, 1].plot(gums[gums_mag].mag_bp - gums[gums_mag].mag_rp, gums[gums_mag].mag_g, '.',
            ms=ms, color=RED, alpha=alpha, **kwargs)
        ax[1, 1].plot(
            gums[gums_mag].mag_bp - gums[gums_mag].abp - gums[gums_mag].mag_rp + gums[gums_mag].arp,
            gums[gums_mag].mag_g - 5 * np.log10(gums[gums_mag].barycentric_distance / 10) - gums[
                gums_mag].ag, '.', color=RED, ms=ms, alpha=alpha, **kwargs)

        ax[0, 2].plot(gaia[gaia_mag].phot_bp_mean_mag - gaia[gaia_mag].phot_rp_mean_mag,
            gaia[gaia_mag].phot_g_mean_mag, '.',color=ORANGE, ms=ms, alpha=alpha, **kwargs)
        ax[1, 2].remove()

        ax[0, 0].invert_yaxis()
        ax[1, 0].set_ylim([15, -5])
        ax[1, 0].set_xlim([-1, 7])

        ax[0, 0].set_ylabel('G$_{obs}$ [mag])')
        ax[1, 0].set_ylabel('G$_{abs}$ [mag]')
        ax[0, 0].set_xlabel('BP$_{obs}$ - RP$_{obs}$ [mag]')
        ax[0, 1].set_xlabel('BP$_{obs}$ - RP$_{obs}$ [mag]')
        ax[0, 2].set_xlabel('BP$_{obs}$ - RP$_{obs}$ [mag]')
        ax[1, 0].set_xlabel('BP$_{abs}$ - RP$_{abs}$ [mag]')
        ax[1, 1].set_xlabel('BP$_{abs}$ - RP$_{abs}$ [mag]')
        ax[1, 2].set_xlabel('BP$_{abs}$ - RP$_{abs}$ [mag]')
        ax[0, 0].set_title('SynthPop')
        ax[0, 1].set_title('Gaia Universe Model')
        ax[0, 2].set_title('Gaia DR3')
        plt.savefig(f"images/{self.loc}_comp_cmd.pdf")

    # Plot density histograms by population
    def plot_dens_hist_pop(self, g_lim=20, bin_width=0.2):
        gums = self.complete_gums
        synth1 = self.complete_synth1
        synth2 = self.complete_synth2
        gums_mag_lim = gums.mag_g < g_lim
        synth1_mag_lim = synth1.Gaia_G_EDR3 < g_lim
        synth2_mag_lim = synth2.Gaia_G_EDR3 < g_lim

        solid_angle = half_cone_angle_to_solidangle(self.radius_deg).value * len(self.l_deg)
        gums_corr = gums.mass < 2.2 * gums.age**(-0.4)+0.26
        for gums_pop_id in range(1, 5):
            for corr in ['','_corr', "_cleaned"]:
                pop_gums = (gums['population'] == gums_pop_id)
                if gums_pop_id == 4:
                    title = 'bulge'
                    synth_pop_id = {0}
                elif gums_pop_id == 3:
                    title = 'halo'
                    synth_pop_id = {8}
                elif gums_pop_id == 2:
                    title = 'thick disk'
                    synth_pop_id = {9, 10}
                elif gums_pop_id == 1:
                    title = 'thin disk'
                    synth_pop_id = {1, 2, 3, 4, 5, 6, 7}

                pop_synth1 = synth1['pop'].isin(synth_pop_id)
                pop_synth2 = synth2['pop'].isin(synth_pop_id)

                plt.figure(f'pop_{gums_pop_id}_number{corr}')
                plt.title(title + ' [number_density]')

                dist = np.linspace(synth1[pop_synth1].Dist.min(),
                    synth1[pop_synth1].Dist.max(),
                    10000)
                test_id = next(iter(synth_pop_id))

                dens = self.get_combined_density(synth_pop_id, dist)

                if self.mod1.populations[test_id].population_density.density_unit == 'number':
                    plt.plot(dist, dens, 'k')

                _dense_hist(1, synth1[pop_synth1].Dist, bin_width, solid_angle,
                    color="slategrey", label=f'synthpop')

                heights, _ = _dense_hist(1, synth1[synth1_mag_lim & pop_synth1].Dist, bin_width, solid_angle,
                    color=BLUE, label=f'synthpop (g<{g_lim:.1f}mag)')
                if corr:
                    _dense_hist(1, synth2[synth2_mag_lim & pop_synth2].Dist, bin_width, solid_angle,
                        color=GREEN, label=f'synthpop corr (g<{g_lim:.1f}mag)')
                _dense_hist(1, gums[gums_mag_lim & pop_gums].barycentric_distance / 1000, bin_width, solid_angle,
                    color=RED, label=f'gums (g<{g_lim:.1f}mag)')
                lim = 10**np.floor(np.log10(np.min(heights[heights>1])))
                plt.xlim([None, min(27, plt.xlim()[1])])
                plt.ylim([lim,None])

                plt.xlabel('Distance [kpc]')
                plt.ylabel('Density [1/kpc$^{3}$]')
                plt.yscale('log')
                plt.legend()
                plt.savefig(f'images/{self.loc}_{title}_number{corr}.pdf', dpi=300)

                plt.figure(f'pop_{gums_pop_id}{corr}')
                plt.title(title)

                if self.mod1.populations[test_id].population_density.density_unit == 'mass':
                    plt.plot(dist, dens, 'k')

                _dense_hist(synth1[pop_synth1].Mass, synth1[pop_synth1].Dist, bin_width, solid_angle,
                    color="slategrey", label=f'synthpop')

                _dense_hist(synth1[synth1_mag_lim & pop_synth1].Mass,
                    (synth1[synth1_mag_lim & pop_synth1].Dist), bin_width,  solid_angle, color=BLUE,
                    label=f'synthpop (g<{g_lim:.1f}mag)')
                if corr == '_corr':

                    _dense_hist(synth2[synth2_mag_lim & pop_synth2].Mass,
                        (synth2[synth2_mag_lim & pop_synth2].Dist), bin_width,  solid_angle, color=GREEN,
                        label=f'synthpop_corr (g<{g_lim:.1f}mag)')
                _dense_hist(gums[gums_mag_lim & pop_gums].mass,
                    (gums[gums_mag_lim & pop_gums].barycentric_distance / 1000), bin_width, solid_angle,
                    color=RED, label=f'gums (g<{g_lim:.1f}mag)')
                if corr == '_cleaned':
                    _dense_hist(gums[gums_mag_lim & pop_gums & gums_corr].mass,
                        (gums[gums_mag_lim & pop_gums & gums_corr].barycentric_distance / 1000), bin_width,
                        solid_angle,
                        color=ORANGE, label=f'gums cleand')

                plt.xlabel('Distance [kpc]')
                plt.ylabel('Density [M$_{\odot}$/kpc$^{3}$]')
                plt.xlim([None, min(27, plt.xlim()[1])])
                plt.yscale('log')
                plt.legend()
                plt.savefig(f'images/{self.loc}_{title}{corr}.pdf', dpi=300)
              
def _plot_density(mass, dist, solid_angle, bin_width, label='', **kwargs):
    bins = np.arange(0, max(dist) + bin_width, bin_width)
    return plt.hist(dist, bins=bins,
        weights=mass/(dist**2*bin_width*solid_angle),
        label=label, histtype="step", **kwargs)

def plot_log_ratio(data, xlim, ylim, vmin=1e-1, vmax=1e1, cmap='bwr', **kwargs):
    data[data > 1e3] = np.nan
    fig, ax = plt.subplots()
    im = ax.imshow(np.rot90(data + 1),
        norm=LogNorm(vmin=vmin, vmax=vmax),
        extent=[*xlim, *ylim],
        cmap=cmap,
        **kwargs)
    cb = plt.colorbar(im, ax=ax)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim[::-1])
    ax.set_ylabel("BP - RP [mag]")
    ax.set_xlabel("G [mag]")
    ax.set_aspect("auto")
    return fig, ax, im

def plot_lin_diff(data, xlim, ylim, vmin=None, vmax=None, cmap='bwr', **kwargs):
    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)
    m = vmin / (vmin - vmax)
    func = lambda v: 0.5 * (v / m) * (v <= m) + (1 - 0.5 * (v - 1) / (m - 1)) * (v > m)
    idx = np.linspace(0, 1, 256)
    cmap = cm.get_cmap(cmap)
    cmap_shift = LinearSegmentedColormap.from_list('test', cmap(func(idx)))
    fig, ax = plt.subplots()
    im = ax.imshow(
        np.rot90(data), extent=[*xlim, *ylim],
        vmin=vmin, vmax=vmax, cmap=cmap_shift, **kwargs)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label('$\Delta$ Stars per mag$^2$')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim[::-1])
    ax.set_xlabel("BB - RP [mag]")
    ax.set_ylabel("G [mag]")
    return fig, ax, im

def _dense_hist(weights, dist, bins, solid_angle, **kwargs):
    print(dist)
    if isinstance(bins, (int, float)):
        bins = np.arange(0, max(dist) + bins, bins)
    if isinstance(weights, (int, float)):
        scale = weights
        weights = None
    else:
        scale = 1.
    heights, edges = np.histogram(dist, bins=bins, weights=weights)
    heights = heights.astype(float) * scale
    heights /= solid_angle / 3 * (bins[1:] ** 3 - bins[:-1] ** 3)
    plt.step(edges, np.append(heights, heights[-1]), where='post', **kwargs)
    return heights, edges

if __name__ == "__main__":
    if len(sys.argv) > 1:
        field = sys.argv[1]
    else:
        field = input('select_field (Bulge/Disk): ')

    if len(sys.argv) > 2:
        user_input = sys.argv[2]
    else:
        user_input = input('regenerate (yes/no): ')
    if user_input.lower() == 'yes':
        regen=True
    else:
        regen=False

    if len(sys.argv) > 3:
        redownload = (sys.argv[3])=='yes'
    else:redownload = True
    
    if not os.path.isdir(f'{DIRNAME}/images'):
        os.mkdir(f'{DIRNAME}/images')

    if field.lower() == 'disk':
        comp_gaia = CompareGaia(loc='Disk', regenerate=regen, radius_deg=0.5,
            redownload=redownload, n_locations=10)
    elif field.lower() == 'bulge':
        comp_gaia = CompareGaia(loc='Bulge', regenerate=regen, radius_deg=0.025,
            redownload=redownload, n_locations=10)

    comp_gaia.plot_hist_G()
    comp_gaia.plot_comp_cmds()
    comp_gaia.plot_dens_hist_pop()
    plt.show()
