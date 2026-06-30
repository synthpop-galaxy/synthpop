#!/usr/bin/env python3
# ============================================================
# generate_ogle_kirkpatrick_min001.py
#
# Kirkpatrick IMF, min_mass=0.01 OGLE workflow:
#   1) patch/check SynthPop config
#   2) generate 3451 OGLE source/lens chip catalogs
#   3) compute mulens_rates_ogle_0tE300_kirkpatrick_min001.txt
#   4) make Fig. 5-style four-statistic plot, residual plots,
#      baseline comparison plots, and statistics CSVs
#
# Typical run:
#   python generate_ogle_kirkpatrick_min001.py
#
# Plot only after rates already exist:
#   python generate_ogle_kirkpatrick_min001.py --only-plot
# ============================================================

from pathlib import Path
import argparse
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import synthpop as sp
import fetch_data
from mulens_rates import microlensing_calculations

warnings.filterwarnings("ignore")


# ============================================================
# Settings / args
# ============================================================
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--config",
        default="./synthpop/config_files/huston2025_defaults_min001_kirkpatrick_min001.synthpop_conf",
        help="Kirkpatrick SynthPop config file.",
    )
    p.add_argument("--tag", default="kirkpatrick_min001")
    p.add_argument("--min-mass", type=float, default=0.01)
    p.add_argument("--output-root", default="outputfiles")
    p.add_argument("--data-dir", default="data/apjsab426b")
    p.add_argument("--baseline-rates", default="mulens_rates_ogle_0tE300.txt")
    p.add_argument("--skip-catalogs", action="store_true")
    p.add_argument("--skip-rates", action="store_true")
    p.add_argument("--only-plot", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--source-maglim", type=float, default=21.0)
    p.add_argument("--lens-maglim", type=float, default=99.0)
    p.add_argument("--solid-angle-start", type=float, default=1e-5)
    p.add_argument("--lower-n", type=float, default=3000.0)
    p.add_argument("--upper-n", type=float, default=10000.0)
    p.add_argument("--target-n", type=float, default=4000.0)
    return p.parse_args()


PROPS = ["eventrate_area", "eventrate_source", "avg_tau", "avg_t", "avg_logt"]
METRICS = [
    {
        "key": "avg_tau",
        "label": r"$\tau$",
        "ylabel": r"$\tau$ (10$^{-6}$)",
        "ratio_label": r"$\tau_{\rm model}/\tau_{\rm OGLE}$",
        "resid_label": r"$\tau$ fractional residual",
        "model_scale": 1e6,
        "ogle_col": "tau",
        "ogle_err": "e_tau",
        "unit": "1e-6",
        "ylim": (0, 3.5),
        "ratio_ylim": (0.5, 2.5),
    },
    {
        "key": "eventrate_source",
        "label": r"$\Gamma_*$",
        "ylabel": r"$\Gamma_*$ (10$^{-6}$ yr$^{-1}$)",
        "ratio_label": r"$\Gamma_{*,\rm model}/\Gamma_{*,\rm OGLE}$",
        "resid_label": r"$\Gamma_*$ fractional residual",
        "model_scale": 1e6,
        "ogle_col": "gam",
        "ogle_err": "e_gam",
        "unit": "1e-6 yr^-1 star^-1",
        "ylim": (0, 35),
        "ratio_ylim": (0.5, 2.5),
    },
    {
        "key": "eventrate_area",
        "label": r"$\Gamma$",
        "ylabel": r"$\Gamma$ (deg$^{-2}$ yr$^{-1}$)",
        "ratio_label": r"$\Gamma_{\rm model}/\Gamma_{\rm OGLE}$",
        "resid_label": r"$\Gamma$ fractional residual",
        "model_scale": 1.0,
        "ogle_col": "gam-deg2",
        "ogle_err": "e_gam-deg2",
        "unit": "deg^-2 yr^-1",
        "ylim": None,
        "ratio_ylim": (0.5, 2.8),
    },
    {
        "key": "avg_t",
        "label": r"$\overline{t_E}$",
        "ylabel": r"$\overline{t_E}$ (days)",
        "ratio_label": r"$\overline{t_E}_{\rm model}/\overline{t_E}_{\rm OGLE}$",
        "resid_label": r"$\overline{t_E}$ fractional residual",
        "model_scale": 1.0,
        "ogle_col": "tE-mean",
        "ogle_err": "e_tE-mean",
        "unit": "days",
        "ylim": None,
        "ratio_ylim": (0.5, 1.5),
    },
]


def tau_fit(b):
    return 1.36 * np.exp(0.39 * (3 - np.abs(b)))


def gamma_star_fit(b):
    return 13.4 * np.exp(0.49 * (3 - np.abs(b)))


# ============================================================
# Config and OGLE tables
# ============================================================
def patch_config(base_conf, tag, min_mass):
    base_conf = Path(base_conf).expanduser().resolve()
    if not base_conf.exists():
        raise FileNotFoundError(f"Could not find config file: {base_conf}")

    conf = json.loads(base_conf.read_text())
    old_min = conf["POPULATION_GENERATION"]["mass_lims"].get("min_mass")
    conf["POPULATION_GENERATION"]["mass_lims"]["min_mass"] = min_mass

    patched = base_conf.parent / f"{base_conf.stem}_patched_{tag}{base_conf.suffix}"
    patched.write_text(json.dumps(conf, indent=2))

    print("Base config:   ", base_conf)
    print("Patched config:", patched)
    print(f"min_mass:      {old_min} -> {min_mass}")
    return patched


def load_ogle_tables(data_dir):
    data_dir = Path(data_dir).expanduser().resolve()
    data_dir.mkdir(parents=True, exist_ok=True)

    ogle_surfdens, ogle_fields, ogle_rates = fetch_data.ogle_mroz2019(
        location=str(data_dir) + "/"
    )
    ogle_surfdens = ogle_surfdens.copy()
    ogle_fields = ogle_fields.copy()
    ogle_rates = ogle_rates.copy()

    ogle_rates["RAdeg"] = ogle_fields["RAdeg"]
    ogle_rates["DEdeg"] = ogle_fields["DEdeg"]
    ogle_rates.sort_values("GLAT", inplace=True)

    ogle_rates["in_lims"] = (
        (np.abs(ogle_rates["GLON"]) < 10)
        & (ogle_rates["GLAT"] < 5)
        & (ogle_rates["GLAT"] > -10)
    )
    ogle_rates["use_for_latdep"] = (
        (np.abs(ogle_rates["GLON"]) < 3)
        & (ogle_rates["Nevents"] > 1)
        & (np.abs(ogle_rates["GLAT"]) > 0.7)
    )

    subfs_inmap = ogle_surfdens[
        (np.abs(ogle_surfdens.GLON) < 10)
        & (ogle_surfdens.GLAT < 5)
        & (ogle_surfdens.GLAT > -10)
    ].copy()

    return ogle_surfdens, ogle_fields, ogle_rates, subfs_inmap


# ============================================================
# Generation and rates
# ============================================================
def expected_catalog_path(out_dir, model_name, l, b):
    return Path(out_dir) / f"{model_name}_l{l:2.3f}_b{b:2.3f}.h5"


def build_models(patched_conf, model_name, src_out, lens_out, args):
    chosen_bands = [
        "Bessell_U", "Bessell_B", "Bessell_V", "Bessell_R", "Bessell_I",
        "2MASS_J", "2MASS_H", "2MASS_Ks",
    ]

    src_mod = sp.SynthPop(
        str(patched_conf),
        maglim=["Bessell_I", args.source_maglim, "remove"],
        chosen_bands=chosen_bands,
        post_processing_kwargs=[{"name": "ProcessDarkCompactObjects", "remove": True}],
        name_for_output=model_name,
        output_location=str(src_out),
        skip_lowmass_stars=True,
        output_file_type="h5",
    )

    lens_mod = sp.SynthPop(
        str(patched_conf),
        maglim=["Bessell_I", args.lens_maglim, "keep"],
        chosen_bands=chosen_bands,
        post_processing_kwargs=[{"name": "ProcessDarkCompactObjects", "remove": False}],
        name_for_output=model_name,
        output_location=str(lens_out),
        output_file_type="h5",
    )
    return src_mod, lens_mod


def run_catalogs(mod, fields_df, label, args):
    mod.init_populations()
    solang = args.solid_angle_start

    for n_done, idx in enumerate(fields_df.index, start=1):
        l = fields_df.loc[idx, "GLON"]
        b = fields_df.loc[idx, "GLAT"]

        df, _ = mod.process_location(l_deg=l, b_deg=b, solid_angle=solang)
        n = len(df)
        if n <= 0:
            raise RuntimeError(f"{label}: zero rows for idx={idx}, l={l:.3f}, b={b:.3f}")

        if n > args.upper_n or n < args.lower_n:
            print(f"{label} rerun: {n_done}/{len(fields_df)} idx={idx}, l={l:.3f}, b={b:.3f}, N={n}")
            solang = solang * args.target_n / n
            df, _ = mod.process_location(l_deg=l, b_deg=b, solid_angle=solang)
            n = len(df)
            if n <= 0:
                raise RuntimeError(f"{label}: zero rows after rerun for idx={idx}, l={l:.3f}, b={b:.3f}")

        solang = solang * args.target_n / n
        print(f"{label} done:  {n_done}/{len(fields_df)} idx={idx}, l={l:.3f}, b={b:.3f}, N={n}")


def compute_rates(subfs_inmap, lens_out, src_out, model_name, rates_out):
    rows = []
    output_cols = None

    for n_done, idx in enumerate(subfs_inmap.index, start=1):
        l = subfs_inmap.loc[idx, "GLON"]
        b = subfs_inmap.loc[idx, "GLAT"]
        f_lens = expected_catalog_path(lens_out, model_name, l, b)
        f_src = expected_catalog_path(src_out, model_name, l, b)

        if not f_lens.exists():
            raise FileNotFoundError(f"Missing lens catalog for idx={idx}: {f_lens}")
        if not f_src.exists():
            raise FileNotFoundError(f"Missing source catalog for idx={idx}: {f_src}")

        try:
            dat, output_cols = microlensing_calculations.mulens_stats(
                l,
                b,
                str(f_lens),
                str(f_src),
                nsd=True,
                field_id=idx,
                tE_range=[0, 300],
            )
            rows.append(dat)
            print(f"RATE ok:   {n_done}/{len(subfs_inmap)} idx={idx}, l={l:.3f}, b={b:.3f}")
        except Exception as e:
            print(f"RATE fail: {n_done}/{len(subfs_inmap)} idx={idx}, l={l:.3f}, b={b:.3f} -> {e}")
            if output_cols is None:
                raise
            failed = {col: np.nan for col in output_cols}
            for col, val in {"l": l, "b": b, "field_id": idx, "f_lens": str(f_lens), "f_src": str(f_src)}.items():
                if col in failed:
                    failed[col] = val
            rows.append([failed[col] for col in output_cols])

    chip_rates = pd.DataFrame(rows, columns=output_cols)
    chip_rates.to_csv(rates_out, index=False)
    print("Saved rates table:", rates_out)
    return chip_rates


# ============================================================
# Aggregation / statistics
# ============================================================
def aggregate_chip_rates_to_fields(ogle_fields, ogle_rates, subfs_inmap, chip_rates):
    if len(chip_rates) != len(subfs_inmap):
        raise ValueError(f"Rates/subfield length mismatch: rates={len(chip_rates)}, subfields={len(subfs_inmap)}")

    chip_rates = chip_rates.copy()
    chip_rates.index = subfs_inmap.index
    sim_fields = ogle_fields[["GLON", "GLAT"]].copy()

    for idx in ogle_fields.index:
        chips = np.array([str(idx) in str(idx_chip) for idx_chip in subfs_inmap.index])
        if np.sum(chips) > 16:
            for prop in PROPS:
                sim_fields.loc[idx, prop] = np.nanmean(pd.to_numeric(chip_rates.loc[chips, prop], errors="coerce"))

    ogle_rates_sel = ogle_rates[ogle_rates["use_for_latdep"]].copy()
    sp_dat = sim_fields.loc[ogle_rates_sel.index].copy()
    return sp_dat, ogle_rates_sel


def metric_arrays(sp_dat, ogle_rates_sel, m):
    model = pd.to_numeric(sp_dat[m["key"]], errors="coerce").to_numpy() * m["model_scale"]
    ogle = pd.to_numeric(ogle_rates_sel[m["ogle_col"]], errors="coerce").to_numpy()
    err = pd.to_numeric(ogle_rates_sel[m["ogle_err"]], errors="coerce").to_numpy()
    return model, ogle, err


def weighted_mean_model_over_ogle(model, ogle, err):
    ratio = model / ogle
    ratio_err = err * model / ogle**2

    # Matches the notebook convention: omit final selected OGLE point.
    ratio = ratio[:-1]
    ratio_err = ratio_err[:-1]

    good = np.isfinite(ratio) & np.isfinite(ratio_err) & (ratio_err > 0)
    if not np.any(good):
        return np.nan
    return np.average(ratio[good], weights=1 / ratio_err[good]**2)


def write_stats(sp_dat, ogle_rates_sel, tag, stats_csv, residual_csv):
    stats_rows = []
    residuals = pd.DataFrame({
        "field_id": ogle_rates_sel.index.astype(str),
        "GLON": ogle_rates_sel["GLON"].to_numpy(),
        "GLAT": ogle_rates_sel["GLAT"].to_numpy(),
    })

    for m in METRICS:
        model, ogle, err = metric_arrays(sp_dat, ogle_rates_sel, m)
        ratio = model / ogle
        frac_resid = ratio - 1.0
        abs_resid = model - ogle

        stats_rows.append({
            "tag": tag,
            "metric": m["key"],
            "unit": m["unit"],
            "n_fields": int(np.sum(np.isfinite(ratio))),
            "weighted_mean_model_over_ogle": weighted_mean_model_over_ogle(model, ogle, err),
            "mean_model_over_ogle": np.nanmean(ratio),
            "median_model_over_ogle": np.nanmedian(ratio),
            "mean_fractional_residual_model_minus_ogle": np.nanmean(frac_resid),
            "median_fractional_residual_model_minus_ogle": np.nanmedian(frac_resid),
            "rms_fractional_residual": np.sqrt(np.nanmean(frac_resid**2)),
            "mean_absolute_residual": np.nanmean(abs_resid),
            "median_absolute_residual": np.nanmedian(abs_resid),
        })

        residuals[f"{m['key']}_model"] = model
        residuals[f"{m['key']}_ogle"] = ogle
        residuals[f"{m['key']}_model_over_ogle"] = ratio
        residuals[f"{m['key']}_frac_resid"] = frac_resid
        residuals[f"{m['key']}_abs_resid"] = abs_resid

    stats = pd.DataFrame(stats_rows)
    stats.to_csv(stats_csv, index=False)
    residuals.to_csv(residual_csv, index=False)

    print("\nWeighted mean ratios, same convention as notebook:")
    for _, row in stats.iterrows():
        print(f"{row['metric']:18s} = {row['weighted_mean_model_over_ogle']:.6g}")
    print("Saved comparison statistics:", stats_csv)
    print("Saved field residuals:       ", residual_csv)
    return stats, residuals


def write_baseline_stats(new_rates, baseline_rates_path, tag, out_csv):
    baseline_rates_path = Path(baseline_rates_path).expanduser().resolve()
    if not baseline_rates_path.exists():
        print("No baseline rates table found; skipping baseline comparison:", baseline_rates_path)
        return None

    base = pd.read_csv(baseline_rates_path).copy()
    new = new_rates.copy()
    for df in [base, new]:
        df["_l_key"] = pd.to_numeric(df["l"], errors="coerce").round(3)
        df["_b_key"] = pd.to_numeric(df["b"], errors="coerce").round(3)

    keep = ["avg_tau", "eventrate_source", "eventrate_area", "avg_t", "avg_logt"]
    if "frac_lowmass" in base.columns and "frac_lowmass" in new.columns:
        keep.append("frac_lowmass")

    cmp = base[["_l_key", "_b_key", "l", "b"] + keep].merge(
        new[["_l_key", "_b_key", "l", "b"] + keep],
        on=["_l_key", "_b_key"],
        suffixes=("_baseline", f"_{tag}"),
        how="inner",
    )

    rows = []
    for metric in keep:
        old = pd.to_numeric(cmp[f"{metric}_baseline"], errors="coerce").to_numpy()
        newv = pd.to_numeric(cmp[f"{metric}_{tag}"], errors="coerce").to_numpy()
        ratio = newv / old
        rows.append({
            "metric": metric,
            "n_matched_chips": int(np.sum(np.isfinite(ratio))),
            f"mean_{tag}_over_baseline": np.nanmean(ratio),
            f"median_{tag}_over_baseline": np.nanmedian(ratio),
            f"mean_percent_change_{tag}_minus_baseline": 100 * np.nanmean(ratio - 1),
            f"median_percent_change_{tag}_minus_baseline": 100 * np.nanmedian(ratio - 1),
            f"mean_absolute_difference_{tag}_minus_baseline": np.nanmean(newv - old),
            f"median_absolute_difference_{tag}_minus_baseline": np.nanmedian(newv - old),
        })

    stats = pd.DataFrame(rows)
    stats.to_csv(out_csv, index=False)
    cmp_out = out_csv.with_name(out_csv.stem.replace("statistics", "matched_chips") + out_csv.suffix)
    cmp.to_csv(cmp_out, index=False)
    print("Saved baseline comparison statistics:", out_csv)
    print("Saved baseline matched-chip table:  ", cmp_out)
    return stats, cmp


# ============================================================
# Plots
# ============================================================
def savefig(fig, pdf_path):
    pdf_path = Path(pdf_path)
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(pdf_path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    print("Saved plot:", pdf_path)
    print("Saved plot:", pdf_path.with_suffix(".png"))


def make_fig5(sp_dat, ogle_rates_sel, out_pdf):
    b_pts = np.arange(-7, 7.0001, 0.1)
    plt.rcParams.update({"font.size": 18})

    fig = plt.figure(figsize=(12, 12))
    outer = gridspec.GridSpec(2, 1, figure=fig, hspace=0.12)
    gs_top = outer[0].subgridspec(2, 2, height_ratios=[3, 1], hspace=0, wspace=0.3)
    gs_bottom = outer[1].subgridspec(2, 2, height_ratios=[3, 1], hspace=0, wspace=0.3)

    ax = [[None, None], [None, None], [None, None], [None, None]]
    ax[0][0] = fig.add_subplot(gs_top[0, 0])
    ax[1][0] = fig.add_subplot(gs_top[1, 0], sharex=ax[0][0])
    ax[0][1] = fig.add_subplot(gs_top[0, 1])
    ax[1][1] = fig.add_subplot(gs_top[1, 1], sharex=ax[0][1])
    ax[2][0] = fig.add_subplot(gs_bottom[0, 0])
    ax[3][0] = fig.add_subplot(gs_bottom[1, 0], sharex=ax[2][0])
    ax[2][1] = fig.add_subplot(gs_bottom[0, 1])
    ax[3][1] = fig.add_subplot(gs_bottom[1, 1], sharex=ax[2][1])

    main_axes = [ax[0][0], ax[0][1], ax[2][0], ax[2][1]]
    ratio_axes = [ax[1][0], ax[1][1], ax[3][0], ax[3][1]]

    for i, m in enumerate(METRICS):
        model, ogle, err = metric_arrays(sp_dat, ogle_rates_sel, m)
        x = ogle_rates_sel["GLAT"].to_numpy()
        ratio = model / ogle
        ratio_err = err * model / ogle**2

        main_axes[i].plot(x, model, "o", label="SP-H25 Kirkpatrick min=0.01")
        main_axes[i].errorbar(x, ogle, yerr=err, linestyle="none", marker="o", c="k", label="OGLE-IV")
        main_axes[i].set_ylabel(m["ylabel"])
        if m["ylim"] is not None:
            main_axes[i].set_ylim(*m["ylim"])

        if m["key"] == "avg_tau":
            main_axes[i].plot(b_pts, tau_fit(b_pts), label="OGLE Fit")
            main_axes[i].legend(fontsize=12)
        elif m["key"] == "eventrate_source":
            main_axes[i].plot(b_pts, gamma_star_fit(b_pts), label="OGLE Fit")

        ratio_axes[i].errorbar(x, ratio, yerr=ratio_err, linestyle="none", marker="o", c="k")
        ratio_axes[i].axhline(1, color="gray", linestyle="--")
        ratio_axes[i].set_ylabel(m["ratio_label"])
        ratio_axes[i].set_ylim(*m["ratio_ylim"])

    ax[3][0].set_xlabel("b (deg.)")
    ax[3][1].set_xlabel("b (deg.)")
    fig.tight_layout()
    savefig(fig, out_pdf)
    plt.close(fig)


def make_ogle_residual_plot(sp_dat, ogle_rates_sel, out_pdf):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    axes = axes.ravel()
    x = ogle_rates_sel["GLAT"].to_numpy()

    for ax, m in zip(axes, METRICS):
        model, ogle, err = metric_arrays(sp_dat, ogle_rates_sel, m)
        resid = model / ogle - 1
        resid_err = err * model / ogle**2
        ax.errorbar(x, resid, yerr=resid_err, linestyle="none", marker="o", c="k")
        ax.axhline(0, color="gray", linestyle="--")
        ax.set_ylabel(m["resid_label"])

    axes[2].set_xlabel("b (deg.)")
    axes[3].set_xlabel("b (deg.)")
    fig.suptitle("Kirkpatrick min=0.01 model − OGLE fractional residuals", y=0.98)
    fig.tight_layout()
    savefig(fig, out_pdf)
    plt.close(fig)


def make_baseline_difference_plot(sp_dat, base_sp_dat, out_pdf):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    axes = axes.ravel()
    x = sp_dat["GLAT"].to_numpy()

    for ax, m in zip(axes, METRICS):
        y = pd.to_numeric(sp_dat[m["key"]], errors="coerce").to_numpy()
        y0 = pd.to_numeric(base_sp_dat[m["key"]], errors="coerce").to_numpy()
        ax.plot(x, y / y0 - 1, "o")
        ax.axhline(0, color="gray", linestyle="--")
        ax.set_ylabel(m["label"] + r" Kirkpatrick/baseline $-1$")

    axes[2].set_xlabel("b (deg.)")
    axes[3].set_xlabel("b (deg.)")
    fig.suptitle("Kirkpatrick min=0.01 − baseline fractional differences", y=0.98)
    fig.tight_layout()
    savefig(fig, out_pdf)
    plt.close(fig)


# ============================================================
# Main
# ============================================================
def main():
    args = parse_args()
    if args.only_plot:
        args.skip_catalogs = True
        args.skip_rates = True

    tag = args.tag
    model_name = f"Huston2025_{tag}"
    root_out = Path(args.output_root).expanduser().resolve() / f"ogle_chips_{tag}"
    src_out = root_out / "src"
    lens_out = root_out / "lens"
    rates_out = Path.cwd() / f"mulens_rates_ogle_0tE300_{tag}.txt"
    subfs_csv = Path.cwd() / f"subfs_inmap_{tag}.csv"

    fig5_pdf = Path.cwd() / f"fig5_ogle_latdep_{tag}.pdf"
    ogle_resid_pdf = Path.cwd() / f"fig5_ogle_fractional_residuals_{tag}.pdf"
    baseline_diff_pdf = Path.cwd() / f"fig5_baseline_fractional_differences_{tag}.pdf"
    stats_csv = Path.cwd() / f"fig5_comparison_statistics_{tag}.csv"
    residual_csv = Path.cwd() / f"fig5_ogle_field_residuals_{tag}.csv"
    baseline_stats_csv = Path.cwd() / f"fig5_baseline_comparison_statistics_{tag}.csv"

    print("Run tag:       ", tag)
    print("Model name:    ", model_name)
    print("Output root:   ", root_out)
    print("Rates table:   ", rates_out)

    patched_conf = patch_config(args.config, tag, args.min_mass)
    _, ogle_fields, ogle_rates, subfs_inmap = load_ogle_tables(args.data_dir)
    subfs_inmap.to_csv(subfs_csv, index=True)
    print(f"{len(subfs_inmap)} chips inside the extinction-map region")
    print("Saved subfield table:", subfs_csv)

    src_out.mkdir(parents=True, exist_ok=True)
    lens_out.mkdir(parents=True, exist_ok=True)

    if not args.skip_catalogs:
        src_mod, lens_mod = build_models(patched_conf, model_name, src_out, lens_out, args)
        print("\n=== Generating source catalogs ===")
        run_catalogs(src_mod, subfs_inmap, "SRC", args)
        print("\n=== Generating lens catalogs ===")
        run_catalogs(lens_mod, subfs_inmap, "LENS", args)
    else:
        print("\nSkipping catalog generation.")

    if not args.skip_rates:
        print("\n=== Computing microlensing rates ===")
        chip_rates = compute_rates(subfs_inmap, lens_out, src_out, model_name, rates_out)
    else:
        print("\nSkipping rate calculation.")
        if not rates_out.exists():
            raise FileNotFoundError(f"Rates file does not exist: {rates_out}")
        chip_rates = pd.read_csv(rates_out)

    sp_dat, ogle_rates_sel = aggregate_chip_rates_to_fields(ogle_fields, ogle_rates, subfs_inmap, chip_rates)
    print(f"\nUsing {len(sp_dat)} plotted OGLE fields")
    print(f"Using {len(chip_rates)} chip-level rates from {rates_out}")

    write_stats(sp_dat, ogle_rates_sel, tag, stats_csv, residual_csv)
    write_baseline_stats(chip_rates, args.baseline_rates, tag, baseline_stats_csv)

    if "frac_lowmass" in chip_rates.columns and "eventrate_source" in chip_rates.columns and "n_source" in chip_rates.columns:
        w = pd.to_numeric(chip_rates["eventrate_source"], errors="coerce") * pd.to_numeric(chip_rates["n_source"], errors="coerce")
        frac = pd.to_numeric(chip_rates["frac_lowmass"], errors="coerce")
        good = np.isfinite(w) & np.isfinite(frac) & (w > 0)
        if np.any(good):
            event_weighted_frac = np.sum(w[good] * frac[good]) / np.sum(w[good])
            print("\nEvent-weighted frac_lowmass:")
            print(f"  {event_weighted_frac:.6f} ({100 * event_weighted_frac:.2f}%)")

    if not args.no_plots:
        make_fig5(sp_dat, ogle_rates_sel, fig5_pdf)
        make_ogle_residual_plot(sp_dat, ogle_rates_sel, ogle_resid_pdf)

        baseline_path = Path(args.baseline_rates).expanduser().resolve()
        if baseline_path.exists():
            base_chip_rates = pd.read_csv(baseline_path)
            if len(base_chip_rates) == len(subfs_inmap):
                base_sp_dat, _ = aggregate_chip_rates_to_fields(ogle_fields, ogle_rates, subfs_inmap, base_chip_rates)
                make_baseline_difference_plot(sp_dat, base_sp_dat, baseline_diff_pdf)
            else:
                print(
                    "Skipping baseline field-difference plot because baseline length "
                    f"({len(base_chip_rates)}) != subfield length ({len(subfs_inmap)})."
                )

    print("\nDone.")
    print("Key outputs:")
    print("  rates:                 ", rates_out)
    print("  Fig. 5 style plot:     ", fig5_pdf)
    print("  OGLE residual plot:    ", ogle_resid_pdf)
    print("  statistics CSV:        ", stats_csv)
    print("  field residuals CSV:   ", residual_csv)
    print("  baseline stats CSV:    ", baseline_stats_csv)


if __name__ == "__main__":
    main()
