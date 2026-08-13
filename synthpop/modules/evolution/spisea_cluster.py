"""
Evolution module to store information for the SpiseaGenerator. Not valid for a standard StarGenerator.
"""

__all__ = ["SpiseaCluster", ]
__author__ = "M.J. Huston"
__date__ = "2025-05-28"

try:
    from ._evolution import EvolutionIsochrones, EvolutionInterpolator, EVOLUTION_DIR
except:
    from _evolution import EvolutionIsochrones, EvolutionInterpolator, EVOLUTION_DIR
import numpy as np
import json
from spisea import evolution as spisea_evolution
from spisea import atmospheres as spisea_atmospheres
from spisea import synthetic as spisea_synthetic
import pdb

class SpiseaCluster(EvolutionIsochrones,EvolutionInterpolator):
    """
    Placeholder object to store modules and values for use by the SpiseaGenerator, which
    will generate and evolve stars as binned SPISEA clusters.
    """
    def __init__(self, columns, n_proc=1, photsys=None,
                    spisea_evolution_name="MISTv1", block_spisea_prints=True,
                    spisea_evolution_kwargs={"version":1.2, "synthpop_extension":True},
                    spisea_atm_func_name="get_merged_atmosphere", spisea_wd_atm_func_name="get_wd_atmosphere",
                    min_mass=0, max_mass=1000, effective_wavelengths='pivot',
                    bbh_frac=0.1, **kwargs):
        self.name='SpiseaCluster'
        if not n_proc>=1:
            raise ValueError("n_proc for SpiseaCluster must be at least 1")
        self.n_proc = n_proc
        self.block_spisea_prints=block_spisea_prints

        self.spisea_evolution = getattr(spisea_evolution, spisea_evolution_name)(**spisea_evolution_kwargs)

        self.spisea_atm_func = getattr(spisea_atmospheres, spisea_atm_func_name)
        self.spisea_wd_atm_func = getattr(spisea_atmospheres, spisea_wd_atm_func_name)
        self.min_mass = min_mass
        self.max_mass = max_mass

        self.allowed_non_mag_cols = ["[Fe/H]_init", "log10_isochrone_age_yr", 'phase',
            'star_mass', 'initial_mass', 'log_L', 'log_R', 'log_Teff', 'log_g', 'isWR']

        self.magsys, self.non_mag_cols, self.bands, self.bands_obs_str = self.get_cols(columns)
        with open(f"{EVOLUTION_DIR}/spisea_effective_wavelengths.json") as f:
            all_eff_wavelengths = json.load(f)[effective_wavelengths]
        self.eff_wavelengths = {self.bands[i]:all_eff_wavelengths[self.bands_obs_str[i]] for i in range(len(self.bands))}

        # Deal with magnitude systems
        self.photsys = photsys
        if self.photsys is None:
            self.photsys = 'Vega' # All mags default to Vega in SPISEA
        if self.photsys in ['Vega','AB','ST']:
            self.photsys_dict = {band: self.photsys for band in self.bands}
        elif isinstance(self.photsys, dict):
            self.photsys_dict = {band: 'Vega' for band in self.bands}
            for key, val in self.photsys.items():
                if key not in ['Vega', 'AB','ST']:
                    raise ValueError("Dict input for photsys must be of structure e.g. {'AB':['m_roman_f146','m_roman_f087']} or {'AB':['roman']} "
                        "with photometric systems being Vega, AB, ST, and using any valid filter sets or names from MIST v1.2.")
                for item in val:
                    if item in self.bands:
                        self.photsys_dict.update({item: key})
                    elif item in self.all_spisea_filters:
                        self.photsys_dict.update({band: key for band in self.bands if ('m_'+item+'_' in band)})
                    else:
                        raise ValueError(f"Invalid input in photsys_dict: {item} not found as band or band set.")
        else:
            raise ValueError("photsys must be 'Vega', 'AB', 'ST', None, or dict")

        with open(f"{EVOLUTION_DIR}/spisea_photometric_system_conversions.json") as f:
            all_photsys_convert = json.load(f)
        self.photsys_convert = {}
        for i,band in enumerate(self.bands):
            if self.photsys_dict[band] == 'Vega':
                self.photsys_convert[band] = 0.0
            else:
                self.photsys_convert[band] = all_photsys_convert[self.photsys_dict[band]][self.bands_obs_str[i]]
                                    
        # Binary evolution
        self.bbh_frac=bbh_frac
        
        if spisea_evolution_name=='MISTv1':
            self.feh_list = np.array([-4.0,-3.5,-3.0,-2.5,-2.0,-1.75,-1.5,-1.25,
                                      -1.0,-0.75,-0.5,-0.25,0,0.25,0.5])
            self.log_age_list = np.linspace(5.0,10.3,54)
            self.log_age_list[0] = 5.01
        elif spisea_evolution_name=='MergedBaraffePisaEkstromParsec':
            self.feh_list = np.array([0.0])
            self.log_age_list = np.linspace(6.0,10.1,83)
            self.log_age_list[-1] = 10.9
            Warning(f"Evolution module {spisea_evolution_name} only includes solar metallicy. All stars will be assigned solar metallicity.")
        elif spisea_evolution_name=='COSMIC':
            self.feh_list = np.array([-2.3,-2.0,-1.75,-1.5,-1.25,
                                      -1.0,-0.75,-0.5,-0.25,0,0.176])
            self.log_age_list = np.linspace(5.0,10.3,54)
            if self.bbh_frac<1.0:
                self.bbh_frac=1.0
                Warning("Setting bbh_frac to 1.0 to let COSMIC handle binary evolution.")
            if spisea_atm_func_name != "get_merged_atmosphere_w_bb_supplement":
                self.spisea_atm_func = spisea_atmospheres.get_merged_atmosphere_w_bb_supplement
                Warning("Setting amosphere model to get_merged_atmosphere_w_bb_supplement for COSMIC.")
            self.spisea_wd_atm_func = None
        else:
            raise ValueError("Invalid SPISEA evolution_model. Only MISTv1,  MergedBaraffePisaEkstromParsec,"
                    " and COSMIC are available at this time.")
            
    def get_cols(self, columns):
        with open(f"{EVOLUTION_DIR}/spisea_filters.json") as f:
            self.all_spisea_filters = json.load(f)
        magsys = {}
        all_bands = []
        non_mag_cols = []

        for column in columns:
            column_split = column.split(',')
            # Non-magnitude columns
            if column in self.allowed_non_mag_cols:
                non_mag_cols.append(column)
            # Magnitude systems
            elif column in self.all_spisea_filters.keys():
                magsys = self.all_spisea_filters[column]
                # Magnitude systems with nested categories & need all
                if hasattr(magsys, "keys"):
                    for magsys_subset in magsys.keys():
                        for band in magsys[magsys_subset]:
                            all_bands.append(column+','+magsys_subset+','+band)
                # Magnitude systems with no nested catagories & need all
                else:
                    for band in magsys:
                        all_bands.append(column+','+band)
            # Magnitude systems with nested or band selections
            elif (len(column_split) in [2,3]) and (column_split[0] in self.all_spisea_filters.keys()):
                magsys = self.all_spisea_filters[column_split[0]]
                # Specified nested band set, use all filters
                if (len(column_split)==2) and hasattr(magsys, "keys"):
                    if column_split[1] in magsys.keys():
                        for band in magsys[column_split[1]]:
                            all_bands.append(column+','+band)
                    else:
                        raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
                # Specified nested band set, select filters
                elif len(column_split)==3:
                    if column_split[1] in magsys.keys():
                        if column_split[2] in magsys[column_split[1]]:
                            all_bands.append(column)
                        else:
                            raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
                    else:
                        raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
                elif len(column_split)==2:
                    if column_split[1] in magsys:
                        all_bands.append(column)
                    else:
                        raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
                else:
                    raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
            else:
                raise ValueError('Invalid column '+column+' for SPISEA isochrones.')

        all_bands_obs_str = list(np.unique(all_bands))
        all_bands = ['m_'+spisea_synthetic.get_filter_col_name(band) for band in all_bands_obs_str]

        return None, list(np.unique(non_mag_cols)), all_bands, all_bands_obs_str

    def get_evolved_props(**kwargs):
        raise ValueError("SpiseaCluster Evolution module is only compatible with SpiseaGenerator, not StarGenerator")
        return

def generate_effective_wavelengths_json():
    import pysynphot
    with open(f'{EVOLUTION_DIR}/spisea_filters.json') as f:
        d = json.load(f)
    def do_filt_effs(sys,flts):
        efflamsi = {}
        for flt in flts:
                flt_name = sys+','+flt
                filt = spisea_synthetic.get_filter_info(flt_name)
                obs = pysynphot.Observation(pysynphot.Vega, filt).efflam()
                assert str(filt.waveunits) == 'angstrom'
                efflamsi[flt_name] = obs*1e-4
                #pdb.set_trace()
                print(flt_name, obs*1e-4)
        return efflamsi
    def do_filt_pivs(sys,flts):
        efflamsi = {}
        for flt in flts:
                flt_name = sys+','+flt
                filt = spisea_synthetic.get_filter_info(flt_name)
                obs = filt.pivot()
                assert str(filt.waveunits) == 'angstrom'
                efflamsi[flt_name] = obs*1e-4
                #pdb.set_trace()
                print(flt_name, obs*1e-4)
        return efflamsi
    def do_filt_avgs(sys,flts):
        efflamsi = {}
        for flt in flts:
                flt_name = sys+','+flt
                filt = spisea_synthetic.get_filter_info(flt_name)
                obs = filt.avgwave()
                assert str(filt.waveunits) == 'angstrom'
                efflamsi[flt_name] = obs*1e-4
                #pdb.set_trace()
                print(flt_name, obs*1e-4)
        return efflamsi

    all_effs = {}
    typ = ['vega_eff','pivot','average']
    funcs = [do_filt_effs, do_filt_pivs, do_filt_avgs]
    for tp in range(3):
        efflams = {}
        do_filts = funcs[tp]
        for sys in d:
            flts = d[sys]
            if isinstance(flts,dict):
                for subsys in flts:
                    efflams.update(do_filts(sys+','+subsys,flts[subsys]))
            else:
                efflams.update(do_filts(sys,flts))
        all_effs[typ[tp]] = efflams

    json_object = json.dumps(all_effs, indent=4)
    with open(f"{EVOLUTION_DIR}/spisea_effective_wavelengths.json", "w") as outfile:
        outfile.write(json_object)
    return

def generate_photsys_conversions():
    with open(f'{EVOLUTION_DIR}/spisea_filters.json') as f:
        d = json.load(f)

    def do_filts(sys,flts):
        convs_ab = {}
        convs_st = {}
        for flt in flts:
            flt_name = sys+','+flt
            convs_ab.update({flt_name: spisea_synthetic.calc_ab_vega_filter_conversion(flt_name)})
            convs_st.update({flt_name: spisea_synthetic.calc_st_vega_filter_conversion(flt_name)})
        return convs_ab, convs_st

    all_convs = {"AB":{}, "ST":{}}
    for sys in d:
        flts = d[sys]
        if isinstance(flts,dict):
            for subsys in flts:
                convs_ab, convs_st = do_filts(sys+','+subsys,flts[subsys])
                all_convs["AB"].update(convs_ab)
                all_convs["ST"].update(convs_st)
        else:
            convs_ab, convs_st = do_filts(sys,flts)
            all_convs["AB"].update(convs_ab)
            all_convs["ST"].update(convs_st)

    json_object = json.dumps(all_convs, indent=4)
    with open(f"{EVOLUTION_DIR}/spisea_photometric_system_conversions.json", "w") as outfile:
        outfile.write(json_object)
    return
