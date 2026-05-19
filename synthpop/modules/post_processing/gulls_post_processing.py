"""
Post-processing to put catalog in the Gulls input format.
"""

__all__ = ["GullsPostProcessing", ]
__author__ = "Farzaneh Zohrabi, Ali Crisp"
__date__ = "2025-05-15"

from ._post_processing import PostProcessing
import time
import pandas as pd
import numpy as np
from synthpop.synthpop_utils.coordinates_transformation import lb_to_ad

class GullsPostProcessing(PostProcessing):
    def __init__(self, model, cat_type, **kwargs):
        """
        Parameters:
            model: SynthPop
                to access properties from the model.
                ege model.const, model.params, model.populations[0]
            cat_type:
                settings for the catalog type
                must be one of : lens, source, extra-bright, bright, mid1, mid2, faint

            kwargs: dict
                keyword arguments from the config file
        """
        self.cat_type = cat_type
        self.model = model

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):

        if companion_df is not None:
            raise ValueError("Must run combine_tables postproc before"+ \
                                "gulls_post_processing")
        filtlist = self.model.parms.chosen_bands
        
        # convert l, b to ra, dec
        system_df["l(deg)"] = system_df['l']  # l
        system_df["b(deg)"] = system_df['b']  # b
        ra, dec = lb_to_ad(system_df["l(deg)"], system_df["b(deg)"])
        system_df["RA2000.0"] = ra
        system_df["DEC2000.0"] = dec

        # convert columns
        system_df["logg"] = system_df["log_g"]
        system_df["Mbol"] = -2.5 * system_df["log_L"] + 4.75
        system_df["Teff"] = 10 ** system_df["log_Teff"]
        system_df["[alpha/Fe]"] = 0
        system_df["Radius"] = 10 ** system_df["log_R"]
        system_df["CL"] = system_df["phase"]
        system_df["Vr"] = system_df['vr_bc']

        # For binary systems, get the primary magnitudes
        if 'Is_Binary' in system_df.columns:
            # Copy original magitude columns to combined magnitude columns
            for filt in filtlist:
                system_df[f"combined_{filt}"] = system_df[filt]
        
            # Add in an index column before copying, this let's you map the primaries back later
            # have to do a weird thing to make sure they're unique -- I think this is because of something
            # that carries over from combine_tables maybe? Since SP doesn't save the dataframe between
            # the post-processing steps, I'm guessing those indices never get reset
            if not system_df.index.is_unique:
                system_df = system_df.reset_index(drop=True)
            system_df["original_df_pos"] = np.arange(len(system_df))
            # Create separate tables for primaries and secondaries -- this makes operations more efficient
            # despite creating a need to re-index results, because they're performed across smaller dateframes
            primaries = system_df[system_df['Is_Binary'] == 1].copy()
            secondaries = system_df[system_df['Is_Binary'] == 2].copy()
        
            # Set the index to the primary ID so it can be used in a join
            # Would've been unnecessary if I didn't need to do the reset_index above, but I can't figure
            # out a better way to make sure they get inserted back into the right place w/o looping over rows
            # (and defeating the purpose of all this)
            primaries = primaries.set_index('primary_ID')
            secondaries = secondaries.set_index('primary_ID')
        
            # Combine primary and secondary tables on primary_ID and add columns for system vs secondary
            joined = primaries.join(secondaries, lsuffix='_sys', rsuffix='_sec', how='inner')
            
            # Get indices for mapping values back
            # System
            idx_system = joined['original_df_pos_sys'].to_numpy()
            row_pos_system = system_df.index.get_indexer(idx_system)
            valid_system = (row_pos_system!= -1)
            
            # Secondaries
            idx_secondary = joined['original_df_pos_sec'].to_numpy()
            row_pos_secondary = system_df.index.get_indexer(idx_secondary)
            valid_secondary = (row_pos_secondary != -1)
            
            # Assign the primaries the secondaries' q and map back
            joined['q_sys'] = joined['q_sec']
            q_col_pos = system_df.columns.get_loc('q')
            if valid_system.any():
                system_df.iloc[row_pos_system[valid_system], q_col_pos] = joined['q_sys'].to_numpy()
            
            # Assign the secondaries the primaries' combined_logP and map back
            joined['combined_logP_sec'] = joined['combined_logP_sys']
            logP_col_pos = system_df.columns.get_loc('combined_logP')
            if valid_secondary.any():
                system_df.iloc[row_pos_secondary[valid_secondary], logP_col_pos] = joined['combined_logP_sec'].to_numpy()
        
            # Separate out primary's magnitude and replace original magnitude
            for filt in filtlist:
                # Set columns, more efficient access this way since we added the suffixes on join
                sys_col = f"{filt}_sys"
                sec_col = f"{filt}_sec"
        
                # Assign the secondaries the combined magnitude
                joined[f"combined_{filt}_sec"] = joined[f"{filt}_sys"]
                
                # Get primary's flux -> magnitude
                flux_system = 10.0 ** (-0.4 * joined[sys_col].values)
                flux_secondary = 10.0 ** (-0.4 * joined[sec_col].values)
                
                with np.errstate(divide='ignore', invalid='ignore'):
                    # mag_primary = -2.5 * np.log10(flux_system - flux_secondary)
                    # TODO: write this more similarly to how single stars are done, if possible
                    # If flux_secondary is NaN for a given object, treat the
                    # system flux as just the primary flux. Otherwise use
                    # the difference (system - secondary).
                    flux_primary = flux_system - flux_secondary
                    mag_primary = np.empty_like(flux_primary) # make an array same length as flux array
                    mag_primary.fill(np.nan) # make nan by default
                    sec_nan_mask = np.isnan(flux_secondary) # make mask for where secondary is NaN
                    # If secondary flux is NaN, calculate primary as just the system
                    if sec_nan_mask.any():
                        mag_primary[sec_nan_mask] = -2.5 * np.log10(flux_system[sec_nan_mask])
                    # If secondary flux isn't NaN, calculate primary normally
                    sec_not_nan_mask = ~sec_nan_mask
                    if sec_not_nan_mask.any():
                        mag_primary[sec_not_nan_mask] = -2.5 * np.log10(flux_primary[sec_not_nan_mask])
                    
                # The error state thing makes it a numpy array, need to go back to Series
                mag_primary = pd.Series(mag_primary, index=joined[sys_col].index)
                
                # Forcing the calculation results in some infs, so need to get
                # rid of those. Also need to make sure any cases where the
                # secondary is the non-dark object are accounted for. Doing this
                # weird pandas masking gives the same behavior as the subtract_magnitudes
                # function in utils_functions, but in a pandas-friendly way that
                # let's me keep everything "vectorized" or whatever
                mag_primary = (mag_primary.replace([np.inf, -np.inf], np.nan)
                                          .mask(np.isclose(joined[sys_col],
                                                           joined[sec_col],
                                                           equal_nan=True))
                              )
                

                # This + conversion back to mag could be one operation, but I think doing it this way
                # keeps it vectorized longer?
                # flux_primary = flux_system - flux_secondary
                # # If flux > 0, convert to magnitude; nan if not
                # mag_primary = np.where(flux_primary > 0,
                #                         -2.5*np.log10(flux_primary),
                #                         np.nan)
                
        
                # I don't fully understand why, but according to StackOverflow, I need to do this to
                # stop it from complaining about equal len keys and values when trying to map back the primaries
                # idx_system = joined['original_df_pos_sys'].to_numpy()
                # row_pos_system = system_df.index.get_indexer(idx_system)
                # valid_system = (row_pos_system != -1)
                filt_col_pos_system = system_df.columns.get_loc(filt)
                if valid_system.any():
                    system_df.iloc[row_pos_system[valid_system], filt_col_pos_system] = mag_primary 
                
                # idx_secondary = joined['original_df_pos_sec'].to_numpy()
                # row_pos_secondary = system_df.index.get_indexer(idx_secondary)
                # valid_secondary = (row_pos_secondary != -1)
                comb_col_pos_secondary = system_df.columns.get_loc(f"combined_{filt}")
                if valid_secondary.any():
                    system_df.iloc[row_pos_secondary[valid_secondary], comb_col_pos_secondary] = joined[f"combined_{filt}_sec"]
   
        # reduce data frame to the needed columns
        cols = filtlist + ["mul", "mub", "Vr", "U", "V", "W", "iMass", 
                "CL", "age", "Teff", "logg", "pop", "Mass", "Mbol", "Radius", 
                "[Fe/H]","l", "b", "RA2000.0", "DEC2000.0", 
                "Dist", "x", "y", "z", "A_Ks", "[alpha/Fe]"]

        # gulls needs additional columns when handling binaries
        # if running catalogs to include in a gulls run that includes binaries,
        # need to include these columns.
        # HOWEVER, it does not need them and will not be able to process them
        # if running gulls simulations that are only single-sources or single
        # (non-planetary) lenses. So, don't want to just universally include these.
        if 'Is_Binary' in system_df.columns:
            binary_cols = ['Is_Binary', 'primary_ID', 'ID',
                           'total_mass', 'q', 'combined_logP', 'eccentricity']
            for filt in filtlist:
                binary_cols.append(f"combined_{filt}")
            
            cols = cols + binary_cols
            
        system_df = system_df[cols]
        
        # compute an approximate magnitude for Roman F213 from 2MASS Ks
        k213 = system_df["2MASS_Ks"] + 1.834505
        # Get index of F184 band and insert K213 after that
        idx = system_df.columns.get_loc("F184")+1
        system_df.insert(idx, "K213", k213)
        
        # Also add to the primary magnitudes if binary
        if 'Is_Binary' in system_df.columns:
            combined_k213 = system_df["combined_2MASS_Ks"] + 1.834505
            idx = system_df.columns.get_loc("combined_F184") + 1
            system_df.insert(idx, "combined_K213", combined_k213)
        
        
        # Replace NaNs with 99 for magnitude columns, and a non-physical
        # very small value for other parameters.
        for filt in filtlist:
            system_df.loc[:, filt] = system_df.loc[:,filt].fillna(99)
        system_df.loc[:, "K213"] = system_df.loc[:,"K213"].fillna(99)
        system_df.loc[:, "Mbol"] = system_df.loc[:,"Mbol"].fillna(99)
        if 'Is_Binary' in system_df.columns:
            for filt in filtlist:
                system_df.loc[:, f"combined_{filt}"] = system_df.loc[:, f"combined_{filt}"].fillna(99)
            system_df.loc[:, "combined_K213"] = system_df.loc[:, "combined_K213"].fillna(99)
        system_df = system_df.fillna(value=2e-50)
        system_df = system_df.replace(-np.inf, 2e-50) # takes care of -inf from combining logL when it's a NaN single or both-NaN binary


        # Impose magnitude lower limits based on Roman expectations
        # Lens has no limits, source and extra bright only have faint limits
        # All faint limits should be handled by the configuration file, so are
        # not set here.
        if self.cat_type == "lens":
            pass

        elif self.cat_type == "source":
            pass
            
        elif self.cat_type == "extra-bright":
            pass
            
        elif self.cat_type == "bright":
            system_df = system_df.loc[system_df["W146"] > 12]

        elif self.cat_type == "mid1":
            system_df = system_df.loc[system_df["W146"] > 17]

        elif self.cat_type == "mid2":
            system_df = system_df.loc[system_df["W146"] > 20]

        elif self.cat_type == "faint":
            system_df = system_df.loc[system_df["W146"] > 24]

        return system_df, companion_df
