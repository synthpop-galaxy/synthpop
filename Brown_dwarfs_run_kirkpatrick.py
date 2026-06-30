#!/usr/bin/env python3
# ============================================================
# generate_ogle_bd_kirkpatrick.py
#
# OGLE rerun for SP-H25 with the Kirkpatrick IMF and min_mass = 0.01
# - patches the SynthPop config
# - downloads/loads OGLE tables
# - creates subfs_inmap.csv
# - generates source catalogs
# - generates lens catalogs
# - computes mulens_rates_ogle_0tE300_<tag>.txt
#
# Assumes these imports already work in your environment:
#   import fetch_data
#   from mulens_rates import microlensing_calculations
#   import synthpop as sp
#
# NOTE:
#   Replace the default --config path below with your actual
#   Kirkpatrick IMF config filename.
# ============================================================

from pathlib import Path
import argparse
import json
import warnings

import numpy as np
import pandas as pd
import synthpop as sp
import fetch_data
from mulens_rates import microlensing_calculations

warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate OGLE source/lens catalogs and microlensing rates for a brown-dwarf min-mass rerun with the Kirkpatrick IMF."
    )
    parser.add_argument(
        "--config",
        default="./synthpop/config_files/huston2025_defaults_min001.synthpop_conf",
        help="Path to the Kirkpatrick IMF SynthPop config file."
    )
    parser.add_argument(
        "--tag",
        default="kirkpatrick_min001",
        help="Run tag appended to output names."
    )
    parser.add_argument(
        "--min-mass",
        type=float,
        default=0.01,
        help="Minimum mass to write into the patched config."
    )
    parser.add_argument(
        "--output-root",
        default="outputfiles",
        help="Root directory for output catalogs."
    )
    return parser.parse_args()


def run_catalogs(mod, fields_df, label, ulim=10000, llim=3000, alim=4000):
    mod.init_populations()
    solang = 1e-5

    for i in fields_df.index:
        l = fields_df.loc[i, "GLON"]
        b = fields_df.loc[i, "GLAT"]

        df, _ = mod.process_location(l_deg=l, b_deg=b, solid_angle=solang)
        leng = len(df)

        if leng > ulim or leng < llim:
            print(f"{label} rerun: idx={i}, l={l:.3f}, b={b:.3f}, N={leng}")
            solang = solang * alim / leng
            df, _ = mod.process_location(l_deg=l, b_deg=b, solid_angle=solang)
            leng = len(df)

        solang = solang * alim / leng
        print(f"{label} done: idx={i}, l={l:.3f}, b={b:.3f}, N={leng}")


def main():
    args = parse_args()

    BASE_CONF = Path(args.config).expanduser().resolve()
    if not BASE_CONF.exists():
        raise FileNotFoundError(f"Could not find config file: {BASE_CONF}")

    RUN_TAG = args.tag
    NEW_MIN_MASS = args.min_mass

    MODEL_NAME = f"Huston2025_{RUN_TAG}"
    ROOT_OUT = Path(args.output_root).expanduser().resolve() / f"ogle_chips_{RUN_TAG}"
    SRC_OUT = ROOT_OUT / "src"
    LENS_OUT = ROOT_OUT / "lens"
    RATES_OUT = Path.cwd() / f"mulens_rates_ogle_0tE300_{RUN_TAG}.txt"
    PATCHED_CONF = BASE_CONF.parent / f"{BASE_CONF.stem}_{RUN_TAG}{BASE_CONF.suffix}"
    SUBFS_CSV = Path.cwd() / f"subfs_inmap_{RUN_TAG}.csv"

    chosen_bands = [
        "Bessell_U", "Bessell_B", "Bessell_V", "Bessell_R", "Bessell_I",
        "2MASS_J", "2MASS_H", "2MASS_Ks"
    ]

    print("Base config:", BASE_CONF)
    print("Run tag:", RUN_TAG)
    print("New min mass:", NEW_MIN_MASS)
    print("Output root:", ROOT_OUT)
    print("Rates file:", RATES_OUT)

    # --------------------------------------------------------
    # Patch config
    # --------------------------------------------------------
    conf = json.loads(BASE_CONF.read_text())
    conf["POPULATION_GENERATION"]["mass_lims"]["min_mass"] = NEW_MIN_MASS
    PATCHED_CONF.write_text(json.dumps(conf, indent=2))
    print("Patched config written to:", PATCHED_CONF)

    # --------------------------------------------------------
    # Ensure OGLE data directory exists
    # --------------------------------------------------------
    DATA_DIR = Path("data/apjsab426b")
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Load OGLE tables and create subfs_inmap
    # --------------------------------------------------------
    ogle_surfdens, ogle_fields, ogle_rates = fetch_data.ogle_mroz2019(
        location=str(DATA_DIR) + "/"
    )

    ogle_surfdens = ogle_surfdens.copy()
    subfs_inmap = ogle_surfdens[
        (np.abs(ogle_surfdens.GLON) < 10)
        & (ogle_surfdens.GLAT < 5)
        & (ogle_surfdens.GLAT > -10)
    ].copy()

    subfs_inmap.to_csv(SUBFS_CSV, index=True)
    print(f"{len(subfs_inmap)} chips inside the extinction-map region")
    print("Saved field list to:", SUBFS_CSV)

    # --------------------------------------------------------
    # Make output directories
    # --------------------------------------------------------
    SRC_OUT.mkdir(parents=True, exist_ok=True)
    LENS_OUT.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Build source and lens models
    # --------------------------------------------------------
    src_mod = sp.SynthPop(
        str(PATCHED_CONF),
        maglim=["Bessell_I", 21, "remove"],
        chosen_bands=chosen_bands,
        post_processing_kwargs=[{"name": "ProcessDarkCompactObjects", "remove": True}],
        name_for_output=MODEL_NAME,
        output_location=str(SRC_OUT),
        skip_lowmass_stars=True,
        output_file_type="h5",
    )

    lens_mod = sp.SynthPop(
        str(PATCHED_CONF),
        maglim=["Bessell_I", 99, "keep"],
        chosen_bands=chosen_bands,
        post_processing_kwargs=[{"name": "ProcessDarkCompactObjects", "remove": False}],
        name_for_output=MODEL_NAME,
        output_location=str(LENS_OUT),
        output_file_type="h5",
    )

    # --------------------------------------------------------
    # Generate catalogs
    # --------------------------------------------------------
    print("\n=== Generating source catalogs ===")
    run_catalogs(src_mod, subfs_inmap, "SRC")

    print("\n=== Generating lens catalogs ===")
    run_catalogs(lens_mod, subfs_inmap, "LENS")

    # --------------------------------------------------------
    # Compute microlensing rates
    # --------------------------------------------------------
    print("\n=== Computing microlensing rates ===")
    chip_rows = []

    output_cols = None
    for i in subfs_inmap.index:
        l = subfs_inmap.loc[i, "GLON"]
        b = subfs_inmap.loc[i, "GLAT"]

        f_lens = LENS_OUT / f"{MODEL_NAME}_l{l:2.3f}_b{b:2.3f}.h5"
        f_src = SRC_OUT / f"{MODEL_NAME}_l{l:2.3f}_b{b:2.3f}.h5"

        try:
            dat, output_cols = microlensing_calculations.mulens_stats(
                l, b,
                str(f_lens),
                str(f_src),
                nsd=True,
                field_id=i,
                tE_range=[0, 300],
            )
            chip_rows.append(dat)
            print(f"RATE ok: idx={i}, l={l:.3f}, b={b:.3f}")
        except Exception as e:
            print(f"RATE fail: idx={i}, l={l:.3f}, b={b:.3f} -> {e}")
            if output_cols is None:
                raise
            chip_rows.append([l, b] + list(np.repeat(np.nan, len(output_cols) - 2)))

    chip_rates = pd.DataFrame(chip_rows, columns=output_cols)
    chip_rates.to_csv(RATES_OUT, index=False)

    print("\nDone.")
    print("Patched config:", PATCHED_CONF)
    print("Source catalogs:", SRC_OUT)
    print("Lens catalogs:", LENS_OUT)
    print("Rates table:", RATES_OUT)


if __name__ == "__main__":
    main()