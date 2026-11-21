"""
pkpd-composition/calibrate_tgi.py

Calibrate and externally validate the tumour-growth output layer of the
PK/PD model against the public LXF A677 gefitinib dataset (Eigenmann et al.
2016; redistributed under BSD-3-Clause).

Primary model
-------------
Exposure is the *measured* mean steady-state plasma concentration for each
dose group, applied as a step during the dosing window (days 3-16):

    C(t) = C_ss(dose)   for 3 <= t <= 16
         = 0            otherwise

Tumour growth (Eigenmann):

    dV/dt = 2 lam0 lam1 V / (2 lam0 V + lam1)  -  kappa C(t) V

This avoids propagating the systematic bias of the per-subject PK
reconstruction into the TGI fit (see the PK diagnostic printed at the end).

Calibration
-----------
1. Fit lam0, lam1 to the vehicle group (kappa = 0).
2. Fit kappa to the 100 mg/kg group.
3. Validate on 6.25 and 25 mg/kg (no fitting).

Run:
    python3 calibrate_tgi.py
"""

from __future__ import annotations

import os
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, minimize_scalar

PALETTE = ["#4D4D4D", "#0072B2", "#D55E00", "#009E73"]
GREY = "#4D4D4D"

# Mean steady-state plasma concentrations measured in the dataset (ng/mL).
C_SS = {6.25: 45.0, 25.0: 267.0, 100.0: 1703.0}
DOSE_START = 3.0
DOSE_END = 16.0


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_data():
    df = pd.read_csv("data/PKPD_ErloAndGefi_LXF.csv", sep=";", na_values=["."])
    g = df[(df["DRUG"] == 2) & (df["YTYPE"] == 2)].copy()
    g = g.sort_values(["#ID", "TIME"])
    params = df[df["DRUG"] == 2].groupby("#ID")[["KA", "V", "KE"]].first()
    return df, g, params


# ---------------------------------------------------------------------------
# Tumour model
# ---------------------------------------------------------------------------

def exposure(dose: float, t: float) -> float:
    if dose <= 0.0:
        return 0.0
    return C_SS[dose] if DOSE_START <= t <= DOSE_END else 0.0


def simulate_tumour(V0: float, dose: float, t_eval: np.ndarray,
                    lam0: float, lam1: float, kappa: float) -> np.ndarray:
    """Eigenmann growth + linear drug kill; V stays non-negative."""
    def rhs(t, V):
        v = max(V[0], 1e-6)
        growth = 2.0 * lam0 * lam1 * v / (2.0 * lam0 * v + lam1)
        kill = kappa * exposure(dose, t) * v
        return [growth - kill]

    sol = solve_ivp(rhs, (0.0, float(t_eval[-1])), [V0], t_eval=t_eval,
                    method="LSODA", rtol=1e-7, atol=1e-9)
    return np.maximum(sol.y[0], 0.0)


def fit_vehicle(g: pd.DataFrame):
    obs = (g[g["DOSE GROUP"] == 0.0]
           .groupby("TIME")["Y"].agg(["mean", "sem"]).reset_index())
    t = obs["TIME"].values
    y = obs["mean"].values
    sem = np.maximum(obs["sem"].values, 1.0)
    V0 = y[0]

    def predict(params):
        lam0, lam1 = params
        return simulate_tumour(V0, 0.0, t, lam0, lam1, 0.0)

    def resid(params):
        return (predict(params) - y) / sem

    res = least_squares(resid, [0.1, 50.0],
                        bounds=([1e-4, 1.0], [2.0, 5000.0]))
    lam0, lam1 = res.x
    return lam0, lam1, {"t": t, "obs": y, "sem": sem, "pred": predict(res.x)}


def group_prediction(g: pd.DataFrame, dose: float,
                     lam0: float, lam1: float, kappa: float):
    sub = g[g["DOSE GROUP"] == dose]
    times = np.array(sorted(sub["TIME"].unique()))
    preds = []
    for sid in sub["#ID"].unique():
        s = sub[sub["#ID"] == sid].sort_values("TIME")
        V0 = float(s["w0"].iloc[0])
        preds.append(simulate_tumour(V0, dose, times, lam0, lam1, kappa))
    pred_mean = np.mean(preds, axis=0)
    obs = sub.groupby("TIME")["Y"].agg(["mean", "sem"]).reindex(times)
    return times, obs["mean"].values, obs["sem"].values, pred_mean


def fit_kappa(g: pd.DataFrame, lam0: float, lam1: float) -> float:
    times, obs, sem, _ = group_prediction(g, 100.0, lam0, lam1, 0.0)
    sem = np.maximum(sem, 5.0)

    def loss(logk):
        kappa = 10.0 ** logk
        _, _, _, pred = group_prediction(g, 100.0, lam0, lam1, kappa)
        return float(np.sum(((pred - obs) / sem) ** 2))

    grid = np.linspace(-5.0, -2.0, 31)
    vals = [loss(x) for x in grid]
    best = grid[int(np.argmin(vals))]
    res = minimize_scalar(loss, bounds=(best - 0.3, best + 0.3),
                          method="bounded")
    return 10.0 ** res.x


def metrics(obs: np.ndarray, pred: np.ndarray) -> Dict[str, float]:
    return {
        "RMSE": float(np.sqrt(np.mean((pred - obs) ** 2))),
        "MAPE": float(np.mean(np.abs(pred - obs) / np.maximum(obs, 1.0)) * 100.0),
        "day30_obs": float(obs[-1]),
        "day30_pred": float(pred[-1]),
    }


# ---------------------------------------------------------------------------
# PK diagnostic: how well does the per-subject PK reconstruction match?
# ---------------------------------------------------------------------------

def dose_events(sub: pd.DataFrame) -> List[Tuple[float, float]]:
    ev = []
    for _, row in sub[sub["DOSE"].notna()].iterrows():
        t0, d = float(row["TIME"]), float(row["DOSE"])
        ev.append((t0, d))
        if pd.notna(row["ADDL"]) and pd.notna(row["II"]):
            for k in range(1, int(row["ADDL"]) + 1):
                ev.append((t0 + float(row["II"]) * k, d))
    return sorted(ev)


def simulate_pk(ka: float, V: float, ke: float,
                events: Sequence[Tuple[float, float]],
                t_end: float, t_grid: np.ndarray) -> np.ndarray:
    def rhs(t, y):
        return [-ka * y[0], ka * y[0] - ke * y[1]]

    y = np.array([0.0, 0.0])
    t = 0.0
    ts = [0.0]
    cs = [0.0]
    for td, dose in events:
        if td > t_end:
            break
        if td > t:
            sol = solve_ivp(rhs, (t, td), y,
                            t_eval=np.linspace(t, td, max(2, int((td - t) * 20) + 1)),
                            method="LSODA", rtol=1e-8, atol=1e-10)
            ts.extend(sol.t[1:].tolist())
            cs.extend((sol.y[1] / V)[1:].tolist())
            y = sol.y[:, -1]
            t = td
        y[0] += dose
    if t < t_end:
        sol = solve_ivp(rhs, (t, t_end), y,
                        t_eval=np.linspace(t, t_end, max(2, int((t_end - t) * 20) + 1)),
                        method="LSODA", rtol=1e-8, atol=1e-10)
        ts.extend(sol.t[1:].tolist())
        cs.extend((sol.y[1] / V)[1:].tolist())
    return np.interp(t_grid, ts, cs)


def pk_diagnostic(df: pd.DataFrame, params: pd.DataFrame) -> pd.DataFrame:
    pk = pd.read_csv("data/PK_LXF_gefi.csv", sep=";", na_values=["."])
    pk = pk.merge(params, on="#ID", how="left")
    rows = []
    for sid in pk["#ID"].unique():
        sub = df[(df["#ID"] == sid) & (df["DRUG"] == 2)].sort_values("TIME")
        if sub.empty:
            continue
        p = params.loc[sid]
        tgrid = np.linspace(0.0, 17.0, 1701)
        C = simulate_pk(float(p["KA"]), float(p["V"]), float(p["KE"]),
                        dose_events(sub), 17.0, tgrid)
        obs = pk[(pk["#ID"] == sid) & (pk["YTYPE"] == 1)]
        for _, r in obs.iterrows():
            pred = np.interp(r["TIME"], tgrid, C)
            rows.append({"dose": r["DOSE GROUP"], "obs": r["Y"],
                         "pred": pred, "ratio": r["Y"] / pred})
    res = pd.DataFrame(rows)
    return (res.groupby("dose")["ratio"]
            .agg(["median", "mean", "count"]).reset_index())


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def panel_label(ax, letter: str) -> None:
    ax.text(-0.17, 1.07, f"({letter})", transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="left")


def make_figure(results, lam0, lam1, kappa) -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Hiragino Sans GB", "Arial Unicode MS",
                           "PingFang HK", "DejaVu Sans"],
        "axes.unicode_minus": False,
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

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05])
    colors = {0.0: PALETTE[0], 6.25: PALETTE[1], 25.0: PALETTE[2], 100.0: PALETTE[3]}

    # (a) vehicle
    ax_a = fig.add_subplot(gs[0, 0])
    r = results[0.0]
    ax_a.errorbar(r["t"], r["obs"], yerr=r["sem"], color=colors[0.0],
                  marker="o", ms=3.5, lw=0, capsize=2, elinewidth=1.0,
                  label="实测")
    ax_a.plot(r["t"], r["pred"], color=colors[0.0], lw=2.0, label="模型")
    ax_a.set_xlabel("时间（天）")
    ax_a.set_ylabel("肿瘤体积（mm$^3$）")
    ax_a.set_xlim(0, 30)
    ax_a.legend(frameon=False, fontsize=7)
    ax_a.set_title("Vehicle 生长（校准）")
    panel_label(ax_a, "a")

    # (b) 100 mg/kg
    ax_b = fig.add_subplot(gs[0, 1])
    r = results[100.0]
    ax_b.axvspan(DOSE_START, DOSE_END, color=GREY, alpha=0.08, lw=0)
    ax_b.errorbar(r["t"], r["obs"], yerr=r["sem"], color=colors[100.0],
                  marker="o", ms=3.5, lw=0, capsize=2, elinewidth=1.0,
                  label="实测")
    ax_b.plot(r["t"], r["pred"], color=colors[100.0], lw=2.0, label="模型")
    ax_b.text(DOSE_START + 0.3, ax_b.get_ylim()[1] * 0.9, "给药",
              fontsize=6.5, color=GREY)
    ax_b.set_xlabel("时间（天）")
    ax_b.set_ylabel("肿瘤体积（mm$^3$）")
    ax_b.set_xlim(0, 30)
    ax_b.legend(frameon=False, fontsize=7)
    ax_b.set_title("100 mg/kg（校准）")
    panel_label(ax_b, "b")

    # (c) validation
    ax_c = fig.add_subplot(gs[1, 0])
    for dose in [6.25, 25.0]:
        r = results[dose]
        ax_c.errorbar(r["t"], r["obs"], yerr=r["sem"], color=colors[dose],
                      marker="o", ms=3, lw=0, capsize=2, elinewidth=0.8,
                      alpha=0.9, label=f"{dose:g} mg/kg 实测")
        ax_c.plot(r["t"], r["pred"], color=colors[dose], lw=1.6, ls="--",
                  label=f"{dose:g} mg/kg 模型")
    ax_c.set_xlabel("时间（天）")
    ax_c.set_ylabel("肿瘤体积（mm$^3$）")
    ax_c.set_xlim(0, 30)
    ax_c.legend(frameon=False, fontsize=6.2, ncol=2)
    ax_c.set_title("亚治疗剂量（验证）")
    panel_label(ax_c, "c")

    # (d) day-30 response
    ax_d = fig.add_subplot(gs[1, 1])
    doses = [0.0, 6.25, 25.0, 100.0]
    x = np.arange(len(doses))
    obs30 = [results[d]["day30_obs"] for d in doses]
    pred30 = [results[d]["day30_pred"] for d in doses]
    w = 0.36
    ax_d.bar(x - w / 2, obs30, w, color=[colors[d] for d in doses],
             edgecolor="white", linewidth=0.8, label="实测")
    ax_d.bar(x + w / 2, pred30, w, color=[colors[d] for d in doses],
             alpha=0.45, edgecolor="white", linewidth=0.8, label="模型")
    vehicle = obs30[0]
    for i, d in enumerate(doses):
        if d != 0:
            tgi_obs = (1 - obs30[i] / vehicle) * 100
            tgi_pred = (1 - pred30[i] / vehicle) * 100
            ax_d.text(i - w / 2, obs30[i] + 45, f"{tgi_obs:+.0f}%",
                      ha="center", fontsize=5.4)
            ax_d.text(i + w / 2, pred30[i] + 45, f"{tgi_pred:+.0f}%",
                      ha="center", fontsize=5.4)
    ax_d.text(0.02, 0.96, "TGI（相对 vehicle）", transform=ax_d.transAxes,
              fontsize=6.0, va="top", color=GREY)
    ax_d.set_xticks(x)
    ax_d.set_xticklabels([f"{d:g}" for d in doses])
    ax_d.set_xlabel("剂量（mg/kg）")
    ax_d.set_ylabel("第 30 天肿瘤体积（mm$^3$）")
    ax_d.set_title("第 30 天响应")
    panel_label(ax_d, "d")

    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    fig.savefig("figs/fig7_tgi_calibration.pdf")
    fig.savefig("figs/fig7_tgi_calibration.png")
    print("wrote figs/fig7_tgi_calibration.pdf/.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    df, g, params = load_data()
    lam0, lam1, vehicle_fit = fit_vehicle(g)
    print(f"vehicle fit: lam0={lam0:.4f} 1/day, lam1={lam1:.2f} mm3/day")

    kappa = fit_kappa(g, lam0, lam1)
    print(f"fitted kappa = {kappa:.3e} mL/(ng day)")

    results = {}
    for dose in [0.0, 6.25, 25.0, 100.0]:
        t, obs, sem, pred = group_prediction(g, dose, lam0, lam1, kappa)
        m = metrics(obs, pred)
        results[dose] = {"t": t, "obs": obs, "sem": sem, "pred": pred, **m}
        print(f"dose {dose:6.2f}: RMSE={m['RMSE']:7.1f}  "
              f"MAPE={m['MAPE']:6.1f}%  "
              f"day30 obs={m['day30_obs']:7.1f} pred={m['day30_pred']:7.1f}")

    make_figure(results, lam0, lam1, kappa)

    pd.DataFrame([
        {"parameter": "lam0", "value": lam0, "unit": "1/day",
         "source": "vehicle fit"},
        {"parameter": "lam1", "value": lam1, "unit": "mm3/day",
         "source": "vehicle fit"},
        {"parameter": "kappa", "value": kappa, "unit": "mL/(ng day)",
         "source": "100 mg/kg fit, measured C_ss"},
    ]).to_csv("figs/fig7_parameters.csv", index=False)
    print("wrote figs/fig7_parameters.csv")

    print("\nPK reconstruction diagnostic (observed / predicted median):")
    print(pk_diagnostic(df, params).to_string(index=False))


if __name__ == "__main__":
    main()
