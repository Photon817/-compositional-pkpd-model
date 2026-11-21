"""
pkpd-composition/figures.py

Paper-style figures for the compositional PK/PD case study.

Run:
    python3 figures.py

Outputs (in ./figs):
    fig1_framework.pdf/png
    fig2_pk_occupancy.pdf/png
    fig3_signal_efficacy.pdf/png
    fig4_dose_schedule_surface.pdf/png
    fig5_lumping.pdf/png
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

from params import GEFITINIB, PD_GEFITINIB, PKParams, PDParams
from models import (
    baseline_pd,
    dose_response_curve,
    free_drug_uM,
    hill,
    kd_eff,
    lumping_analysis,
    pk_matrix,
    reduced_pk_profiles,
    simulate_regimen,
)


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
GREY = "#4D4D4D"


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "lines.linewidth": 1.6,
        "lines.markersize": 4,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
    })


def panel_label(ax, letter: str, x: float = -0.18, y: float = 1.08) -> None:
    s = f"({letter})"
    if hasattr(ax, "text2D"):
        ax.text2D(x, y, s, transform=ax.transAxes,
                  fontsize=11, fontweight="bold", va="top", ha="left")
    else:
        ax.text(x, y, s, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top", ha="left")


def save(fig, name: str, outdir: str = "figs") -> None:
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"))
    print(f"wrote {outdir}/{name}.pdf/.png")


# ---------------------------------------------------------------------------
# Helper: summary metrics from a regimen
# ---------------------------------------------------------------------------

def regimen_metrics(sim: Dict[str, np.ndarray], pk: PKParams,
                    window_h: float = 24.0) -> Dict[str, float]:
    t = sim["t"]
    mask = t >= t[-1] - window_h
    C = sim["D_uM"][mask]
    occ = sim["occupancy"][mask]
    perk = sim["pERK_inhibition"][mask]
    apop = sim["apoptosis_fold"][mask]
    return {
        "Cmax_uM": float(np.max(C)),
        "Cmin_uM": float(np.min(C)),
        "AUC_uM_h": float(np.trapezoid(C, t[mask])),
        "occ_max": float(np.max(occ)),
        "occ_min": float(np.min(occ)),
        "perk_max": float(np.max(perk)),
        "perk_min": float(np.min(perk)),
        "apop_max": float(np.max(apop)),
        "apop_auc": float(np.trapezoid(apop, t[mask])),
    }


# ---------------------------------------------------------------------------
# Figure 2: PK and target occupancy
# ---------------------------------------------------------------------------

def fig2_pk_occupancy(pk: PKParams, pd: PDParams) -> None:
    doses = [100.0, 250.0, 500.0]
    colors = [PALETTE[0], PALETTE[1], PALETTE[2]]
    labels = [f"{int(d)} mg q24h" for d in doses]

    sims = [simulate_regimen(pk, pd, dose_mg=d, tau_h=24.0,
                             n_days=14.0, t_step=0.2) for d in doses]
    metrics = [regimen_metrics(s, pk) for s in sims]

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    # (a) total plasma concentration
    for sim, c, lab in zip(sims, colors, labels):
        ax_a.plot(sim["t"] / 24.0, sim["C_total"], color=c, lw=1.3, label=lab)
    ax_a.set_xlabel("time (days)")
    ax_a.set_ylabel("total plasma conc. (mg/L)")
    ax_a.set_xlim(0, 14)
    ax_a.legend(frameon=True, facecolor="white", edgecolor="none",
                framealpha=0.9, loc="upper left", ncol=1, fontsize=6.5)
    ax_a.set_title("Total drug exposure")
    panel_label(ax_a, "a")

    # (b) free concentration, with Kd_eff line
    for sim, c, lab in zip(sims, colors, labels):
        ax_b.plot(sim["t"] / 24.0, sim["D_uM"], color=c, lw=1.3, label=lab)
    ax_b.axhline(kd_eff(pd), color=GREY, ls="--", lw=1.0)
    ax_b.text(0.8, kd_eff(pd) * 1.08, r"$K_d^{eff}$", color=GREY, fontsize=8)
    ax_b.set_xlabel("time (days)")
    ax_b.set_ylabel("free plasma conc. ($\\mu$M)")
    ax_b.set_xlim(0, 14)
    ax_b.set_title("Free drug and target binding")
    panel_label(ax_b, "b")

    # (c) occupancy
    for sim, c, lab in zip(sims, colors, labels):
        ax_c.plot(sim["t"] / 24.0, sim["occupancy"], color=c, lw=1.3, label=lab)
    ax_c.axhline(0.5, color=GREY, ls=":", lw=1.0)
    ax_c.set_xlabel("time (days)")
    ax_c.set_ylabel("target occupancy")
    ax_c.set_ylim(0, 0.85)
    ax_c.set_xlim(0, 14)
    ax_c.set_title("EGFR occupancy")
    panel_label(ax_c, "c")

    # (d) summary bars: Cmax, occupancy max, pERK inhibition max
    x = np.arange(len(doses))
    w = 0.26
    vals_cmax = [m["Cmax_uM"] for m in metrics]
    vals_occ = [m["occ_max"] for m in metrics]
    vals_perk = [m["perk_max"] for m in metrics]
    ax_d.bar(x - w, vals_cmax, w, color=PALETTE[0], label="Cmax free ($\\mu$M)")
    ax_d.bar(x, vals_occ, w, color=PALETTE[1], label="max occupancy")
    ax_d.bar(x + w, vals_perk, w, color=PALETTE[2], label="max pERK inh.")
    ax_d.set_xticks(x)
    ax_d.set_xticklabels([f"{int(d)}" for d in doses])
    ax_d.set_xlabel("dose (mg q24h)")
    ax_d.set_ylabel("value")
    ax_d.set_ylim(0, 1.0)
    ax_d.legend(frameon=False, fontsize=6.5, loc="upper right")
    ax_d.set_title("Exposure-response summary")
    panel_label(ax_d, "d")

    fig.tight_layout()
    save(fig, "fig2_pk_occupancy")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: signal and efficacy
# ---------------------------------------------------------------------------

def fig3_signal_efficacy(pk: PKParams, pd: PDParams) -> None:
    doses = [100.0, 250.0, 500.0]
    colors = [PALETTE[0], PALETTE[1], PALETTE[2]]
    labels = [f"{int(d)} mg q24h" for d in doses]

    sims = [simulate_regimen(pk, pd, dose_mg=d, tau_h=24.0,
                             n_days=14.0, t_step=0.2) for d in doses]

    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0])

    # (a) steady-state dose-response
    ax_a = fig.add_subplot(gs[0, 0])
    D = np.logspace(-3, 0.3, 300)
    perk = np.array([1.0 - pd_steady_state_y(pd, d) / baseline_pd(pd)[3] for d in D])
    ax_a.semilogx(D, perk, color=PALETTE[0], lw=2.0)
    ax_a.axhline(0.5, color=GREY, ls=":", lw=1.0)
    ax_a.axvspan(0.05, 0.15, color=PALETTE[1], alpha=0.12, lw=0)
    ax_a.text(0.055, 0.78, "clinical\nCmax range", fontsize=6.8, color=PALETTE[1])
    ax_a.set_xlabel("free drug ($\\mu$M)")
    ax_a.set_ylabel("pERK inhibition")
    ax_a.set_ylim(0, 1.0)
    ax_a.set_title("Steady-state signal response")
    panel_label(ax_a, "a")

    # (b) pERK time courses at steady state (last 72 h)
    ax_b = fig.add_subplot(gs[0, 1])
    for sim, c, lab in zip(sims, colors, labels):
        t = sim["t"]
        mask = t >= t[-1] - 72.0
        ax_b.plot((t[mask] - t[mask][0]) / 24.0, sim["pERK_inhibition"][mask],
                  color=c, lw=1.4, label=lab)
    ax_b.set_xlabel("time in steady-state window (days)")
    ax_b.set_ylabel("pERK inhibition")
    ax_b.set_ylim(0, 1.0)
    ax_b.legend(frameon=False, fontsize=6.8, loc="upper right")
    ax_b.set_title("pERK dynamics")
    panel_label(ax_b, "b")

    # (c) peak/trough band vs dose
    ax_c = fig.add_subplot(gs[1, 0])
    dose_grid = np.array([25, 50, 75, 100, 150, 200, 250, 350, 500, 750, 1000])
    peak = []
    trough = []
    for d in dose_grid:
        sim = simulate_regimen(pk, pd, dose_mg=d, tau_h=24.0,
                               n_days=14.0, t_step=0.2)
        m = regimen_metrics(sim, pk)
        peak.append(m["perk_max"])
        trough.append(m["perk_min"])
    peak = np.array(peak)
    trough = np.array(trough)
    ax_c.fill_between(dose_grid, trough, peak, color=PALETTE[0], alpha=0.16, lw=0)
    ax_c.plot(dose_grid, peak, "o-", color=PALETTE[0], lw=1.8, label="peak")
    ax_c.plot(dose_grid, trough, "s-", color=PALETTE[0], alpha=0.55, lw=1.2, label="trough")
    ax_c.axvline(250, color=GREY, ls="--", lw=1.0)
    ax_c.set_xscale("log")
    ax_c.set_xticks([25, 50, 100, 250, 500, 1000])
    ax_c.set_xticklabels(["25", "50", "100", "250", "500", "1000"])
    ax_c.set_xlabel("dose (mg q24h)")
    ax_c.set_ylabel("pERK inhibition")
    ax_c.set_ylim(0, 1.0)
    ax_c.legend(frameon=False, fontsize=6.8, loc="upper left")
    ax_c.set_title("Daily peak-trough window")
    panel_label(ax_c, "c")

    # (d) apoptosis fold vs dose
    ax_d = fig.add_subplot(gs[1, 1])
    dr = dose_response_curve(pk, pd, doses_mg=dose_grid, n_days=14.0)
    ax_d.plot(dr["dose_mg"], dr["apoptosis_max"], "o-", color=PALETTE[3], lw=1.8)
    ax_d.axvline(250, color=GREY, ls="--", lw=1.0)
    ax_d.text(260, 1.1, "clinical\ndose", fontsize=6.8, color=GREY)
    ax_d.set_xscale("log")
    ax_d.set_xticks([25, 50, 100, 250, 500, 1000])
    ax_d.set_xticklabels(["25", "50", "100", "250", "500", "1000"])
    ax_d.set_xlabel("dose (mg q24h)")
    ax_d.set_ylabel("max apoptosis fold")
    ax_d.set_title("Apoptotic output")
    panel_label(ax_d, "d")

    fig.tight_layout()
    save(fig, "fig3_signal_efficacy")
    plt.close(fig)


def pd_steady_state_y(pd: PDParams, D_uM: float) -> float:
    from models import pd_steady_state
    return pd_steady_state(pd, D_uM)["Y"]


# ---------------------------------------------------------------------------
# Figure 4: dose-schedule-efficacy surface
# ---------------------------------------------------------------------------

def fig4_surface(pk: PKParams, pd: PDParams) -> None:
    doses = np.array([50, 65, 80, 100, 125, 150, 180, 220, 250, 300, 350,
                      425, 500, 650, 800, 1000])
    intervals = np.array([6.0, 8.0, 12.0, 24.0, 48.0])
    Z = np.zeros((len(intervals), len(doses)))
    for i, tau in enumerate(intervals):
        for j, d in enumerate(doses):
            sim = simulate_regimen(pk, pd, dose_mg=d, tau_h=tau,
                                   n_days=14.0, t_step=0.5)
            Z[i, j] = regimen_metrics(sim, pk)["apop_auc"]

    fig = plt.figure(figsize=(7.2, 5.6))

    # (a) 3D surface
    ax_a = fig.add_axes([0.01, 0.50, 0.41, 0.44], projection="3d")
    X, Y = np.meshgrid(np.log10(doses), intervals)
    ax_a.plot_surface(X, Y, Z, cmap="viridis", linewidth=0,
                      antialiased=True, alpha=0.95)
    ax_a.set_xlabel("dose (mg)", labelpad=3, fontsize=7)
    ax_a.set_ylabel("interval (h)", labelpad=3, fontsize=7)
    ax_a.set_zlabel("apoptosis AUC", labelpad=5, fontsize=7)
    ax_a.set_xticks(np.log10([50, 100, 250, 500, 1000]))
    ax_a.set_xticklabels(["50", "100", "250", "500", "1000"], fontsize=6)
    ax_a.set_yticks([6, 12, 24, 48])
    ax_a.set_yticklabels(["6", "12", "24", "48"], fontsize=6)
    ax_a.tick_params(pad=0)
    ax_a.view_init(elev=22, azim=-125)

    # (b) contour map
    ax_b = fig.add_axes([0.55, 0.56, 0.34, 0.33])
    cs = ax_b.contourf(np.log10(doses), intervals, Z, levels=12, cmap="viridis")
    ax_b.contour(np.log10(doses), intervals, Z, levels=6, colors="white",
                 linewidths=0.6, alpha=0.7)
    ax_b.plot(np.log10(250), 24, "o", color=PALETTE[1], ms=7,
              markeredgecolor="white", markeredgewidth=0.8)
    ax_b.text(np.log10(250) + 0.04, 28.5, "clinical\nregimen", color=PALETTE[1],
              fontsize=6.5)
    ax_b.set_xticks(np.log10([50, 100, 250, 500, 1000]))
    ax_b.set_xticklabels(["50", "100", "250", "500", "1000"])
    ax_b.set_yticks(intervals)
    ax_b.set_yticklabels([f"{int(v)}" for v in intervals])
    ax_b.set_xlabel("dose (mg)")
    ax_b.set_ylabel("dosing interval (h)")
    cb = fig.colorbar(cs, ax=ax_b, fraction=0.046, pad=0.03)
    cb.set_label("apoptosis AUC", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    # (c) interface choice: free vs total vs a 2x-scaled free interface
    ax_c = fig.add_axes([0.08, 0.08, 0.36, 0.33])
    dose_scan = np.array([50, 100, 150, 250, 350, 500, 750, 1000])
    apop_free, apop_total, apop_half = [], [], []
    for d in dose_scan:
        sim = simulate_regimen(pk, pd, dose_mg=d, tau_h=24.0, n_days=14.0)
        apop_free.append(regimen_metrics(sim, pk)["apop_max"])
        sim_t = _simulate_with_interface_scale(pk, pd, d, scale=1.0 / pk.fu)
        apop_total.append(regimen_metrics(sim_t, pk)["apop_max"])
        sim_h = _simulate_with_interface_scale(pk, pd, d, scale=2.0)
        apop_half.append(regimen_metrics(sim_h, pk)["apop_max"])
    ax_c.plot(dose_scan, apop_free, "o-", color=PALETTE[2], label="free drug (correct)")
    ax_c.plot(dose_scan, apop_total, "s--", color=PALETTE[1], label="total drug (wrong)")
    ax_c.plot(dose_scan, apop_half, "^:", color=PALETTE[3], label="2x free (wrong)")
    ax_c.set_xscale("log")
    ax_c.set_xticks([50, 100, 250, 500, 1000])
    ax_c.set_xticklabels(["50", "100", "250", "500", "1000"])
    ax_c.set_xlabel("dose (mg q24h)")
    ax_c.set_ylabel("max apoptosis fold")
    ax_c.set_ylim(0.8, 3.2)
    ax_c.legend(frameon=False, fontsize=6.2, loc="lower right")
    ax_c.set_title("Interface choice changes the prediction", fontsize=8.5, pad=5)

    # (d) same total daily dose, different schedules
    ax_d = fig.add_axes([0.58, 0.08, 0.36, 0.33])
    schedules = [(250.0, 24.0), (125.0, 12.0), (250.0 / 3.0, 8.0), (62.5, 6.0)]
    labels = ["250 q24h", "125 q12h", "83 q8h", "62.5 q6h"]
    cols = [PALETTE[0], PALETTE[1], PALETTE[2], PALETTE[3]]
    rows = []
    for (dose, tau), lab, c in zip(schedules, labels, cols):
        sim = simulate_regimen(pk, pd, dose_mg=dose, tau_h=tau,
                               n_days=14.0, t_step=0.2)
        t = sim["t"]
        mask48 = t >= t[-1] - 48.0
        ax_d.plot(t[mask48] - t[mask48][0], sim["pERK_inhibition"][mask48],
                  color=c, lw=1.4, label=lab)
        mask24 = t >= t[-1] - 24.0
        y = sim["pERK_inhibition"][mask24]
        rows.append((lab, float(np.mean(y >= 0.15) * 24.0), float(np.mean(y))))
    ax_d.axhline(0.15, color=GREY, ls="--", lw=1.0)
    ax_d.text(0.0, 0.151, "threshold 0.15", fontsize=5.5, color=GREY)
    ax_d.set_xlabel("time in steady-state window (h)")
    ax_d.set_ylabel("pERK inhibition")
    ax_d.set_ylim(0.12, 0.22)
    ax_d.legend(frameon=False, fontsize=5.8, loc="upper right", ncol=2)
    ax_d.set_title("Same daily dose, different schedules", fontsize=8.5, pad=5)
    txt = "\n".join(f"{lab}: {hours:.1f} h/day" for lab, hours, _ in rows)
    ax_d.text(0.02, 0.04, txt, transform=ax_d.transAxes, fontsize=5.4,
              va="bottom", ha="left")

    # panel letters and titles
    fig.text(0.005, 0.975, "(a)", fontsize=11, fontweight="bold", va="top")
    fig.text(0.22, 0.975, "Dose-schedule-efficacy surface", fontsize=9.5,
             ha="center", va="top")
    fig.text(0.535, 0.975, "(b)", fontsize=11, fontweight="bold", va="top")
    fig.text(0.72, 0.975, "Iso-effect landscape", fontsize=9.5,
             ha="center", va="top")
    fig.text(0.005, 0.445, "(c)", fontsize=11, fontweight="bold", va="top")
    fig.text(0.535, 0.445, "(d)", fontsize=11, fontweight="bold", va="top")

    save(fig, "fig4_dose_schedule_surface")
    plt.close(fig)


def figS1_sensitivity(pk: PKParams, pd: PDParams) -> None:
    """One-at-a-time parameter sensitivity of the apoptosis AUC."""
    base = _apop_auc_for(pk, pd, 250, 24, 1.0, 1.0, 1.0, 1.0, 1.0)
    factors = ["CL", "Vc", "fu", "Kd", "ka"]
    lo, hi = [], []
    for fac in factors:
        vals = []
        for mult in (0.5, 2.0):
            kw = dict(cl=1.0, vc=1.0, fu=1.0, kd=1.0, ka=1.0)
            kw[{"CL": "cl", "Vc": "vc", "fu": "fu", "Kd": "kd", "ka": "ka"}[fac]] = mult
            vals.append(_apop_auc_for(pk, pd, 250, 24, **kw))
        lo.append((min(vals) / base - 1.0) * 100.0)
        hi.append((max(vals) / base - 1.0) * 100.0)
    order = np.argsort([abs(h - l) for l, h in zip(lo, hi)])
    ypos = np.arange(len(factors))
    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    ax.barh(ypos, [lo[i] for i in order], color=PALETTE[0], height=0.6, label="-50%")
    ax.barh(ypos, [hi[i] for i in order], color=PALETTE[1], height=0.6, label="+50%")
    ax.axvline(0, color=GREY, lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels([factors[i] for i in order])
    ax.set_xlabel("change in apoptosis AUC (%)")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    fig.tight_layout()
    save(fig, "figS1_sensitivity")
    plt.close(fig)


def _simulate_with_interface_scale(pk: PKParams, pd: PDParams,
                                   dose: float, tau: float = 24.0,
                                   scale: float = 1.0) -> Dict[str, np.ndarray]:
    """Simulate the composite but multiply the PD input by `scale`.

    This is the interface-ablation experiment: scale=1 is free drug,
    scale=1/fu is total drug, scale=2 is a two-fold overestimate.
    """
    # easiest: simulate the PK/PD with a modified copy of the PK params
    # whose fu is multiplied by scale; the PK itself is unchanged.
    from dataclasses import replace
    pk_scaled = replace(pk, fu=pk.fu * scale)
    return simulate_regimen(pk_scaled, pd, dose_mg=dose, tau_h=tau, n_days=14.0)


def _apop_auc_for(pk: PKParams, pd: PDParams, dose: float, tau: float,
                  cl: float = 1.0, vc: float = 1.0, fu: float = 1.0,
                  kd: float = 1.0, ka: float = 1.0) -> float:
    from dataclasses import replace
    pk2 = replace(pk, CL=pk.CL * cl, Vc=pk.Vc * vc, fu=pk.fu * fu, ka=pk.ka * ka)
    # scale Kd by adjusting koff (keeping kon fixed): Kd_eff = (koff+kint)/kon
    # simpler: scale the PD input by 1/kd via the interface, which is
    # equivalent to scaling Kd for the occupancy curve.
    sim = simulate_regimen(pk2, pd, dose_mg=dose, tau_h=tau, n_days=14.0)
    # approximate Kd scaling by post-processing the occupancy is not exact;
    # instead, rescale the drug input to the PD module.
    if abs(kd - 1.0) > 1e-9:
        sim = _simulate_with_interface_scale(pk2, pd, dose, tau, scale=1.0 / kd)
    return regimen_metrics(sim, pk2)["apop_auc"]


# ---------------------------------------------------------------------------
# Figure 5: exact lumpability of the PK model
# ---------------------------------------------------------------------------

def fig5_lumping(pk: PKParams, pd: PDParams) -> None:
    info = lumping_analysis(pk)
    t = np.linspace(0.0, 72.0, 721)
    red = reduced_pk_profiles(pk, 250.0, t, y0=(250.0 / pk.Vc, 0.0))

    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05])

    # (a) full vs exact vs naive, linear 0-24 h + log inset
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.semilogy(t, red["C_full"], color=PALETTE[0], lw=2.2,
                  label="full 2-compartment")
    ax_a.semilogy(t, red["C_exact"], color=PALETTE[1], lw=1.4, ls="--",
                  label="exact slow-mode")
    ax_a.semilogy(t, red["C_naive"], color=PALETTE[2], lw=1.4, ls=":",
                  label="naive 1-compartment")
    ax_a.set_xlim(0, 72)
    ax_a.set_ylim(1e-5, 2e-3)
    ax_a.set_xlabel("time (h)")
    ax_a.set_ylabel("central conc. (mg/L)")
    ax_a.legend(frameon=False, fontsize=6.5, loc="upper right")
    ax_a.set_title("PK reduction: exact vs naive")
    ax_a.text(0.02, 0.03,
              "commutator residual\nexact: 3.1e-17\nnaive: 2.7e-2",
              transform=ax_a.transAxes, fontsize=6.0, va="bottom")
    panel_label(ax_a, "a")

    # (b) eigenvalue separation over Q and Vp
    q_grid = np.linspace(5, 200, 80)
    vp_grid = np.linspace(50, 1500, 80)
    QQ, VV = np.meshgrid(q_grid, vp_grid)
    RATIO = np.zeros_like(QQ)
    for i in range(QQ.shape[0]):
        for j in range(QQ.shape[1]):
            from dataclasses import replace
            p = replace(pk, Q=QQ[i, j], Vp=VV[i, j])
            A = pk_matrix(p)
            ev = np.sort(np.abs(np.linalg.eigvals(A)))
            RATIO[i, j] = ev[0] / ev[1] if ev[1] > 0 else np.nan
    ax_b = fig.add_subplot(gs[0, 1])
    cs = ax_b.contourf(QQ, VV, np.log10(RATIO), levels=20, cmap="viridis")
    ax_b.contour(QQ, VV, RATIO, levels=[0.1], colors=PALETTE[1], linewidths=1.6)
    ax_b.plot(pk.Q, pk.Vp, "o", color=PALETTE[1], ms=7,
              markeredgecolor="white", markeredgewidth=0.8)
    ax_b.text(pk.Q + 4, pk.Vp + 60, "clinical PK", color=PALETTE[1], fontsize=6.5)
    ax_b.set_xlabel("inter-compartmental clearance Q (L/h)")
    ax_b.set_ylabel("peripheral volume Vp (L)")
    cb = fig.colorbar(cs, ax=ax_b, fraction=0.046, pad=0.03)
    cb.set_label("log10($|\\lambda_{slow}|/|\\lambda_{fast}|$)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    from matplotlib.lines import Line2D
    ax_b.legend(handles=[Line2D([0], [0], color=PALETTE[1], lw=1.6,
                                label="10x separation")],
                frameon=False, fontsize=6.5, loc="lower right")
    ax_b.set_title("Eigenvalue separation")
    panel_label(ax_b, "b")

    # (c) trajectory error vs Q
    ax_c = fig.add_subplot(gs[1, 0])
    q_scan = np.linspace(5, 150, 40)
    err_exact = []
    err_naive = []
    for q in q_scan:
        from dataclasses import replace
        p = replace(pk, Q=q)
        tt = np.linspace(0, 72, 361)
        rr = reduced_pk_profiles(p, 250.0, tt, y0=(250.0 / p.Vc, 0.0))
        denom = np.trapezoid(rr["C_full"], tt)
        err_exact.append(np.trapezoid(np.abs(rr["C_full"] - rr["C_exact"]), tt) / denom)
        err_naive.append(np.trapezoid(np.abs(rr["C_full"] - rr["C_naive"]), tt) / denom)
    ax_c.semilogy(q_scan, np.maximum(err_exact, 1e-18), color=PALETTE[0],
                  label="exact slow-mode")
    ax_c.semilogy(q_scan, err_naive, color=PALETTE[2], label="naive total-amount")
    ax_c.axvline(pk.Q, color=GREY, ls="--", lw=1.0)
    ax_c.axvspan(5, 20, color=GREY, alpha=0.10, lw=0)
    ax_c.text(7, 0.55, "poor\nseparation", fontsize=6.2, color=GREY)
    ax_c.set_xlabel("Q (L/h)")
    ax_c.set_ylabel("relative concentration error")
    ax_c.legend(frameon=False, fontsize=6.5)
    ax_c.set_title("Trajectory error of the reduction")
    panel_label(ax_c, "c")

    # (d) downstream error: delta pERK inhibition
    ax_d = fig.add_subplot(gs[1, 1])
    from scipy.integrate import solve_ivp
    from models import composite_rhs, baseline_pd
    t2 = np.linspace(0, 72, 361)
    y0 = np.array([0.0, 500.0, 0.0, *baseline_pd(pd)])
    sol_full = solve_ivp(composite_rhs, (0, 72), y0, t_eval=t2,
                         args=(pk, pd), method="LSODA", rtol=1e-8, atol=1e-10)
    def naive_rhs(tt, y, pk, pd):
        G, A1, R, DR, X, Y, A = y
        C = A1 / pk.Vd
        D = pk.fu * C * 1000.0 / pk.MW
        dG = -pk.ka * G
        dA1 = pk.ka * G - (pk.CL / pk.Vd) * A1
        dR = pd.ksyn - pd.kdeg * R - pd.kon * D * R + pd.koff * DR
        dDR = pd.kon * D * R - pd.koff * DR - pd.kint * DR
        dX = pd.kx * hill(R, pd.KR, pd.nX) - pd.kdx * X
        dY = pd.ky * hill(X, pd.KX, pd.nY) - pd.kdy * Y
        dA = pd.ka_apop * (1.0 - hill(Y, pd.KY, pd.nA)) - pd.kda * A
        return [dG, dA1, dR, dDR, dX, dY, dA]
    y0n = np.array([0.0, 500.0, *baseline_pd(pd)])
    sol_naive = solve_ivp(naive_rhs, (0, 72), y0n, t_eval=t2,
                          args=(pk, pd), method="LSODA", rtol=1e-8, atol=1e-10)
    perk_full = 1.0 - sol_full.y[6] / baseline_pd(pd)[3]
    perk_naive = 1.0 - sol_naive.y[5] / baseline_pd(pd)[3]
    delta = perk_naive - perk_full
    ax_d.fill_between(t2, 0, delta, color=PALETTE[2], alpha=0.25, lw=0)
    ax_d.plot(t2, delta, color=PALETTE[2], lw=1.6)
    ax_d.axhline(0, color=GREY, lw=0.8)
    ax_d.set_xlabel("time (h)")
    ax_d.set_ylabel("$\\Delta$ pERK inhibition (naive - full)")
    ax_d.set_title("Downstream error of the naive reduction")
    panel_label(ax_d, "d")

    fig.tight_layout()
    save(fig, "fig5_lumping")
    plt.close(fig)
    print("eigenvalues:", info["eig"])
    print("slow left eigenvector:", info["l_slow"])
    print("exact residual:", info["residual_exact"])
    print("naive residual:", info["residual_naive"])


# ---------------------------------------------------------------------------
# Figure 1: framework schematic
# ---------------------------------------------------------------------------

def fig1_framework(pk: PKParams, pd: PDParams) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # --- PK module (left) -------------------------------------------------
    pk_box = FancyBboxPatch((0.03, 0.42), 0.30, 0.46, boxstyle="round,pad=0.012",
                            fc="#EAF2F8", ec=PALETTE[0], lw=1.2)
    ax.add_patch(pk_box)
    ax.text(0.18, 0.84, "PK module", ha="center", va="center",
            fontsize=9, fontweight="bold", color=PALETTE[0])
    # compartments
    ax.add_patch(FancyBboxPatch((0.06, 0.62), 0.08, 0.12, boxstyle="round,pad=0.006",
                                fc="white", ec=PALETTE[0], lw=0.9))
    ax.text(0.10, 0.68, "gut", ha="center", va="center", fontsize=7)
    ax.add_patch(FancyBboxPatch((0.18, 0.62), 0.10, 0.12, boxstyle="round,pad=0.006",
                                fc="white", ec=PALETTE[0], lw=0.9))
    ax.text(0.23, 0.68, "central", ha="center", va="center", fontsize=7)
    ax.add_patch(FancyBboxPatch((0.18, 0.46), 0.10, 0.10, boxstyle="round,pad=0.006",
                                fc="white", ec=PALETTE[0], lw=0.9))
    ax.text(0.23, 0.51, "periph.", ha="center", va="center", fontsize=7)
    ax.annotate("", xy=(0.18, 0.68), xytext=(0.14, 0.68),
                arrowprops=dict(arrowstyle="->", lw=0.9, color=PALETTE[0]))
    ax.annotate("", xy=(0.23, 0.62), xytext=(0.23, 0.56),
                arrowprops=dict(arrowstyle="<->", lw=0.9, color=PALETTE[0]))
    ax.annotate("", xy=(0.10, 0.74), xytext=(0.05, 0.80),
                arrowprops=dict(arrowstyle="->", lw=1.2, color=PALETTE[0]))
    ax.text(0.035, 0.82, "dose", fontsize=6.5, color=PALETTE[0])

    # --- interface --------------------------------------------------------
    ax.add_patch(FancyBboxPatch((0.37, 0.55), 0.10, 0.12, boxstyle="round,pad=0.008",
                                fc="#FFF4E6", ec=PALETTE[1], lw=1.4))
    ax.text(0.42, 0.61, "$C_{free}$", ha="center", va="center",
            fontsize=9, color=PALETTE[1], fontweight="bold")
    ax.annotate("", xy=(0.37, 0.61), xytext=(0.33, 0.61),
                arrowprops=dict(arrowstyle="->", lw=1.4, color=PALETTE[1]))
    ax.text(0.42, 0.71, "interface", ha="center", fontsize=6.5, color=PALETTE[1])

    # --- PD module (right) ------------------------------------------------
    pd_box = FancyBboxPatch((0.50, 0.30), 0.47, 0.58, boxstyle="round,pad=0.012",
                            fc="#E8F8F5", ec=PALETTE[2], lw=1.2)
    ax.add_patch(pd_box)
    ax.text(0.735, 0.84, "PD module", ha="center", va="center",
            fontsize=9, fontweight="bold", color=PALETTE[2])
    # cascade boxes
    labels = ["EGFR", "MEK", "ERK", "apoptosis"]
    xs = [0.55, 0.66, 0.77, 0.88]
    for x, lab in zip(xs, labels):
        ax.add_patch(FancyBboxPatch((x, 0.60), 0.08, 0.10,
                                    boxstyle="round,pad=0.006",
                                    fc="white", ec=PALETTE[2], lw=0.9))
        ax.text(x + 0.04, 0.65, lab, ha="center", va="center", fontsize=6.8)
    for x0, x1 in zip(xs[:-1], xs[1:]):
        ax.annotate("", xy=(x1, 0.65), xytext=(x0 + 0.08, 0.65),
                    arrowprops=dict(arrowstyle="->", lw=0.9, color=PALETTE[2]))
    ax.annotate("", xy=(0.55, 0.65), xytext=(0.47, 0.61),
                arrowprops=dict(arrowstyle="->", lw=1.4, color=PALETTE[1]))
    ax.text(0.735, 0.42, "drug binding → occupancy → signal → effect",
            ha="center", fontsize=6.5, color=GREY)

    # --- cospan equation --------------------------------------------------
    ax.text(0.5, 0.16, "composite  =  PK  $\\otimes_{C_{free}}$  PD",
            ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.5, 0.06, "pushout on the scalar interface  +  shared-coordinate dynamics",
            ha="center", va="center", fontsize=7.5, color=GREY)

    fig.tight_layout()
    save(fig, "fig1_framework")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    setup_style()
    pk, pd = GEFITINIB, PD_GEFITINIB
    fig1_framework(pk, pd)
    fig2_pk_occupancy(pk, pd)
    fig3_signal_efficacy(pk, pd)
    fig4_surface(pk, pd)
    fig5_lumping(pk, pd)
    figS1_sensitivity(pk, pd)


if __name__ == "__main__":
    main()
