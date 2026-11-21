"""
pkpd-composition/data_explore.py

Raw-data overview of the public Eigenmann et al. (2016) LXF PK/PD dataset,
as redistributed in the BSD-3-Clause repository DavAug/ErlotinibGefitinib.

Data files (downloaded into ./data):
    PKPD_ErloAndGefi_LXF.csv   combined tumour-volume + dosing + PK parameters
    PK_LXF_gefi.csv            gefitinib plasma concentrations
    PK_LXF_erlo.csv            erlotinib plasma concentrations

Run:
    python3 data_explore.py
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PALETTE = ["#4D4D4D", "#0072B2", "#D55E00", "#009E73"]
GREY = "#4D4D4D"


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.0,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.6,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
    })


def panel_label(ax, letter: str) -> None:
    ax.text(-0.16, 1.07, f"({letter})", transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="left")


def main() -> None:
    setup_style()
    os.makedirs("figs", exist_ok=True)

    df = pd.read_csv("data/PKPD_ErloAndGefi_LXF.csv", sep=";", na_values=["."])
    # gefitinib = DRUG 2
    g = df[(df["DRUG"] == 2) & (df["YTYPE"] == 2)].copy()
    g = g.sort_values(["#ID", "TIME"])

    pk = pd.read_csv("data/PK_LXF_gefi.csv", sep=";", na_values=["."])
    pko = pk[pk["YTYPE"] == 1].copy()
    pko["phase"] = np.where(pko["TIME"] < 13, "day 10", "day 16")

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0])

    # (a) tumour growth curves
    ax_a = fig.add_subplot(gs[0, 0])
    dose_colors = {0.0: PALETTE[0], 6.25: PALETTE[1], 25.0: PALETTE[2], 100.0: PALETTE[3]}
    for dose in [0.0, 6.25, 25.0, 100.0]:
        sub = g[g["DOSE GROUP"] == dose]
        c = dose_colors[dose]
        for _, s in sub.groupby("#ID"):
            ax_a.plot(s["TIME"], s["Y"], color=c, alpha=0.18, lw=0.8)
        stats = sub.groupby("TIME")["Y"].agg(["mean", "sem"]).reset_index()
        ax_a.errorbar(stats["TIME"], stats["mean"], yerr=stats["sem"],
                      color=c, lw=2.0, marker="o", ms=3.5, capsize=2,
                      label=f"{dose:g} mg/kg")
    ax_a.set_xlabel("time (days)")
    ax_a.set_ylabel("tumour volume (mm$^3$)")
    ax_a.set_xlim(0, 30)
    ax_a.legend(frameon=False, fontsize=6.8, ncol=2)
    ax_a.set_title("Gefitinib tumour growth inhibition (LXF A677)")
    panel_label(ax_a, "a")

    # (b) PK concentrations
    ax_b = fig.add_subplot(gs[0, 1])
    phases = ["day 10", "day 16"]
    xpos = {6.25: 0, 25.0: 1, 100.0: 2}
    for dose in [6.25, 25.0, 100.0]:
        sub = pko[pko["DOSE GROUP"] == dose]
        c = dose_colors[dose]
        for i, ph in enumerate(phases):
            vals = sub[sub["phase"] == ph]["Y"].values
            x = xpos[dose] + (i - 0.5) * 0.28
            ax_b.scatter(np.full_like(vals, x, dtype=float), vals,
                         color=c, s=18, alpha=0.8, edgecolor="white", linewidth=0.3)
            if len(vals):
                ax_b.plot([x - 0.06, x + 0.06], [np.mean(vals)] * 2,
                          color="black", lw=1.4)
    ax_b.set_xticks([0, 1, 2])
    ax_b.set_xticklabels(["6.25", "25", "100"])
    ax_b.set_xlabel("dose (mg/kg)")
    ax_b.set_ylabel("plasma concentration (ng/mL)")
    ax_b.set_yscale("log")
    ax_b.set_title("Gefitinib exposure at steady state")
    panel_label(ax_b, "b")

    # (c) body-weight change
    ax_c = fig.add_subplot(gs[1, 0])
    bw = df[df["DRUG"] == 2].sort_values(["#ID", "TIME"]).copy()
    bw["baseline"] = bw.groupby("#ID")["BW"].transform("first")
    bw["delta"] = bw["BW"] - bw["baseline"]
    for dose in [0.0, 6.25, 25.0, 100.0]:
        sub = bw[bw["DOSE GROUP"] == dose]
        stats = sub.groupby("TIME")["delta"].agg(["mean", "sem"]).reset_index()
        ax_c.plot(stats["TIME"], stats["mean"], color=dose_colors[dose],
                  lw=1.8, marker="o", ms=3, label=f"{dose:g} mg/kg")
        ax_c.fill_between(stats["TIME"],
                          stats["mean"] - stats["sem"],
                          stats["mean"] + stats["sem"],
                          color=dose_colors[dose], alpha=0.12, lw=0)
    ax_c.axhline(0, color=GREY, lw=0.8, ls=":")
    ax_c.set_xlabel("time (days)")
    ax_c.set_ylabel("body-weight change (g)")
    ax_c.set_xlim(0, 30)
    ax_c.legend(frameon=False, fontsize=6.8, ncol=2)
    ax_c.set_title("Tolerability")
    panel_label(ax_c, "c")

    # (d) day-30 tumour volume and TGI
    ax_d = fig.add_subplot(gs[1, 1])
    last_day = g["TIME"].max()
    last = g[g["TIME"] == last_day]
    means = last.groupby("DOSE GROUP")["Y"].mean()
    sems = last.groupby("DOSE GROUP")["Y"].sem()
    doses = [0.0, 6.25, 25.0, 100.0]
    x = np.arange(len(doses))
    vals = [means.get(d, np.nan) for d in doses]
    errs = [sems.get(d, np.nan) for d in doses]
    colors = [dose_colors[d] for d in doses]
    ax_d.bar(x, vals, yerr=errs, color=colors, width=0.62,
             capsize=3, edgecolor="white", linewidth=0.8)
    vehicle = means.get(0.0, np.nan)
    for i, d in enumerate(doses):
        if d != 0 and np.isfinite(vehicle):
            tgi = (1 - means[d] / vehicle) * 100
            ax_d.text(i, vals[i] + errs[i] + 60, f"TGI\n{tgi:+.0f}%",
                      ha="center", fontsize=6.2)
    ax_d.set_xticks(x)
    ax_d.set_xticklabels([f"{d:g}" for d in doses])
    ax_d.set_xlabel("dose (mg/kg)")
    ax_d.set_ylabel(f"tumour volume at day {int(last_day)} (mm$^3$)")
    ax_d.set_title("Day-30 tumour volume and TGI")
    panel_label(ax_d, "d")

    fig.tight_layout()
    fig.savefig("figs/fig6_external_dataset.pdf")
    fig.savefig("figs/fig6_external_dataset.png")
    print("wrote figs/fig6_external_dataset.pdf/.png")

    # summary table
    print("\nDay-30 tumour volume by dose:")
    print(last.groupby("DOSE GROUP")["Y"].agg(["mean", "sem", "count"]).to_string())
    print("\nPK concentration by dose/phase:")
    print(pko.groupby(["DOSE GROUP", "phase"])["Y"].agg(["mean", "sem", "count"]).to_string())


if __name__ == "__main__":
    main()
