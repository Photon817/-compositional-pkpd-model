"""
pkpd-composition/models.py

Composite PK/PD model, dosing simulation, and exact-lumpability analysis.

State vector
------------
    y = [G, A1, A2, R, DR, X, Y, A]

    G   drug amount in the gut depot (mg)
    A1  drug amount in the central compartment (mg)
    A2  drug amount in the peripheral compartment (mg)
    R   free EGFR (uM)
    DR  drug-bound EGFR (uM)
    X   phospho-MEK (normalized)
    Y   phospho-ERK (normalized)
    A   apoptosis marker (normalized)

The PK module is a two-compartment model with first-order oral absorption.
The PD module is a reduced EGFR -> MEK -> ERK -> apoptosis cascade.  The
interface between the two modules is the free drug concentration

    C_free = fu * A1 / Vc        (mg/L)
    D      = C_free * 1000 / MW  (uM)

which is the shared variable of the decorated cospan.  In the composite
ODE, D is not a state: it is the image of the PK state under the interface
map.  This is the concrete meaning of "gluing on a scalar interface".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from params import PDParams, PKParams


# ---------------------------------------------------------------------------
# Basic helpers
# ---------------------------------------------------------------------------

def hill(x: float, K: float, n: int) -> float:
    """Hill function x^n / (K^n + x^n), safe at x=0."""
    if x <= 0.0:
        return 0.0
    xn = x ** n
    return xn / (K ** n + xn)


def free_drug_mgL(A1: float, pk: PKParams) -> float:
    """Free drug concentration in mg/L from central amount in mg."""
    return pk.fu * A1 / pk.Vc


def free_drug_uM(A1: float, pk: PKParams) -> float:
    """Free drug concentration in uM from central amount in mg."""
    return free_drug_mgL(A1, pk) * 1000.0 / pk.MW


def kd_eff(pd: PDParams) -> float:
    """Effective steady-state occupancy half-saturation.

    With kint = kdeg the total receptor pool is conserved, and the
    occupancy is D / (D + (koff + kint)/kon).
    """
    return (pd.koff + pd.kint) / pd.kon


# ---------------------------------------------------------------------------
# PD steady state at a fixed free drug concentration
# ---------------------------------------------------------------------------

def pd_steady_state(pd: PDParams, D_uM: float) -> Dict[str, float]:
    """Analytic steady state of the PD module for constant free drug D."""
    K = kd_eff(pd)
    R = pd.R_total / (1.0 + D_uM / K)
    DR = pd.R_total - R
    X = pd.kx * hill(R, pd.KR, pd.nX) / pd.kdx
    Y = pd.ky * hill(X, pd.KX, pd.nY) / pd.kdy
    A = pd.ka_apop * (1.0 - hill(Y, pd.KY, pd.nA)) / pd.kda
    return {"R": R, "DR": DR, "X": X, "Y": Y, "A": A}


def baseline_pd(pd: PDParams) -> Tuple[float, float, float, float, float]:
    s = pd_steady_state(pd, 0.0)
    return s["R"], s["DR"], s["X"], s["Y"], s["A"]


# ---------------------------------------------------------------------------
# Composite right-hand side
# ---------------------------------------------------------------------------

def composite_rhs(t: float, y: Sequence[float], pk: PKParams, pd: PDParams):
    """Right-hand side of the composite PK/PD model.

    The PK part is written in amounts (mg); the PD part uses the free
    drug concentration (uM) obtained from the central amount.
    """
    G, A1, A2, R, DR, X, Y, A = y
    D = free_drug_uM(A1, pk)

    # PK: oral absorption + two-compartment distribution + elimination
    dG = -pk.ka * G
    dA1 = pk.ka * G - pk.k10 * A1 - pk.k12 * A1 + pk.k21 * A2
    dA2 = pk.k12 * A1 - pk.k21 * A2

    # PD: drug-receptor binding
    dR = pd.ksyn - pd.kdeg * R - pd.kon * D * R + pd.koff * DR
    dDR = pd.kon * D * R - pd.koff * DR - pd.kint * DR

    # PD: MEK -> ERK -> apoptosis
    dX = pd.kx * hill(R, pd.KR, pd.nX) - pd.kdx * X
    dY = pd.ky * hill(X, pd.KX, pd.nY) - pd.kdy * Y
    dA = pd.ka_apop * (1.0 - hill(Y, pd.KY, pd.nA)) - pd.kda * A

    return [dG, dA1, dA2, dR, dDR, dX, dY, dA]


def rhs_for_lumping(t: float, x: Sequence[float], pk: PKParams):
    """Two-compartment PK only (amounts A1, A2), for the T6 analysis."""
    A1, A2 = x
    return [
        -pk.k10 * A1 - pk.k12 * A1 + pk.k21 * A2,
        pk.k12 * A1 - pk.k21 * A2,
    ]


# ---------------------------------------------------------------------------
# Dosing simulation
# ---------------------------------------------------------------------------

def simulate_regimen(
    pk: PKParams,
    pd: PDParams,
    dose_mg: float | None = None,
    tau_h: float | None = None,
    n_days: float = 14.0,
    t_step: float = 0.25,
) -> Dict[str, np.ndarray]:
    """Simulate a repeated oral dosing regimen.

    Returns a dict of arrays: t, G, A1, A2, C_total, C_free, D_uM,
    occupancy, R, DR, X, Y, A, pERK_inhibition, apoptosis_fold.
    """
    if dose_mg is None:
        dose_mg = pk.dose_mg
    if tau_h is None:
        tau_h = pk.tau_h

    t_end = n_days * 24.0
    y0 = np.array([0.0, 0.0, 0.0, *baseline_pd(pd)], dtype=float)

    t_all: List[float] = []
    y_all: List[np.ndarray] = []

    t = 0.0
    y = y0.copy()
    dose_times = np.arange(0.0, t_end + 1e-9, tau_h)

    for k, td in enumerate(dose_times):
        # bolus into the gut at the dose time
        if abs(t - td) > 1e-9:
            # integrate from the current time to the dose time
            seg_end = td
            n = max(2, int(np.ceil((seg_end - t) / t_step)) + 1)
            t_eval = np.linspace(t, seg_end, n)
            sol = solve_ivp(
                composite_rhs, (t, seg_end), y, t_eval=t_eval,
                args=(pk, pd), method="LSODA", rtol=1e-8, atol=1e-10,
            )
            t_all.extend(sol.t.tolist())
            y_all.extend(sol.y.T.tolist())
            y = sol.y[:, -1]
            t = seg_end
        y = y.copy()
        y[0] += dose_mg * pk.F

        seg_end = min(td + tau_h, t_end)
        if seg_end <= t + 1e-9:
            continue
        n = max(2, int(np.ceil((seg_end - t) / t_step)) + 1)
        t_eval = np.linspace(t, seg_end, n)
        sol = solve_ivp(
            composite_rhs, (t, seg_end), y, t_eval=t_eval,
            args=(pk, pd), method="LSODA", rtol=1e-8, atol=1e-10,
        )
        t_all.extend(sol.t.tolist())
        y_all.extend(sol.y.T.tolist())
        y = sol.y[:, -1]
        t = seg_end

    t_arr = np.asarray(t_all)
    y_arr = np.asarray(y_all)
    # de-duplicate time points (segment boundaries)
    keep = np.concatenate([[True], np.diff(t_arr) > 1e-9])
    t_arr = t_arr[keep]
    y_arr = y_arr[keep]

    G, A1, A2, R, DR, X, Y, A = (y_arr[:, i] for i in range(8))
    C_total = A1 / pk.Vc
    C_free = pk.fu * C_total
    D_uM = C_free * 1000.0 / pk.MW
    occupancy = DR / (R + DR + 1e-30)

    base = baseline_pd(pd)
    pERK_inh = 1.0 - Y / base[3]
    apoptosis_fold = A / base[4]

    return {
        "t": t_arr, "G": G, "A1": A1, "A2": A2,
        "C_total": C_total, "C_free": C_free, "D_uM": D_uM,
        "occupancy": occupancy, "R": R, "DR": DR,
        "X": X, "Y": Y, "A": A,
        "pERK_inhibition": pERK_inh,
        "apoptosis_fold": apoptosis_fold,
    }


# ---------------------------------------------------------------------------
# Steady-state dose-response (for the PK/PD curve)
# ---------------------------------------------------------------------------

def dose_response_curve(
    pk: PKParams,
    pd: PDParams,
    doses_mg: Iterable[float] = (50, 100, 150, 250, 350, 500, 750, 1000),
    tau_h: float = 24.0,
    n_days: float = 14.0,
) -> Dict[str, np.ndarray]:
    """Run a regimen to (near) steady state and return summary metrics."""
    rows = []
    for dose in doses_mg:
        sim = simulate_regimen(pk, pd, dose_mg=dose, tau_h=tau_h,
                               n_days=n_days, t_step=0.5)
        # last 24 h window
        t = sim["t"]
        mask = t >= t[-1] - 24.0
        C_free = sim["C_free"][mask]
        occ = sim["occupancy"][mask]
        perk = sim["pERK_inhibition"][mask]
        apop = sim["apoptosis_fold"][mask]
        rows.append({
            "dose_mg": dose,
            "Cmax_free_mgL": float(np.max(C_free)),
            "Cmax_free_uM": float(np.max(C_free) * 1000.0 / pk.MW),
            "AUC_free_mgL_h": float(np.trapezoid(C_free, t[mask])),
            "occupancy_max": float(np.max(occ)),
            "occupancy_min": float(np.min(occ)),
            "pERK_inh_max": float(np.max(perk)),
            "pERK_inh_min": float(np.min(perk)),
            "apoptosis_max": float(np.max(apop)),
            "apoptosis_AUC": float(np.trapezoid(apop, t[mask])),
        })
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


# ---------------------------------------------------------------------------
# T6: exact lumpability of the two-compartment PK model
# ---------------------------------------------------------------------------

def pk_matrix(pk: PKParams) -> np.ndarray:
    return np.array([
        [-(pk.k10 + pk.k12), pk.k21],
        [pk.k12, -pk.k21],
    ])


def lumping_analysis(pk: PKParams) -> Dict[str, object]:
    """Eigen-decomposition and exact/naive 1D lumpings of the PK model.

    Returns
    -------
    A         2x2 amount matrix
    eig       eigenvalues (1/h)
    V         right eigenvectors (columns)
    l_slow    slow left eigenvector (row)
    lam_slow  slow eigenvalue
    exact     dict with the exact 1D lumping data
    naive     dict with the total-amount 1D lumping data
    """
    A = pk_matrix(pk)
    eig, V = np.linalg.eig(A)
    order = np.argsort(np.abs(eig))  # slow first
    eig = eig[order]
    V = V[:, order]

    # left eigenvectors are right eigenvectors of A^T
    eigT, WT = np.linalg.eig(A.T)
    orderT = np.argsort(np.abs(eigT))
    eigT = eigT[orderT]
    WT = WT[:, orderT]

    l_slow = WT[:, 0].real
    v_slow = V[:, 0].real
    # normalize l . v = 1
    l_slow = l_slow / (l_slow @ v_slow)
    lam_slow = eig[0].real
    lam_fast = eig[1].real

    # exact 1D lumping: T = l_slow . x, dT/dt = lam_slow T
    L_exact = l_slow.reshape(1, 2)
    R_exact = L_exact @ A - lam_slow * L_exact

    # naive total-amount lumping: T = A1 + A2, k_eff = CL / Vd
    # The lumped ODE is dT/dt = -k_eff T, so the matrix condition is
    # L A = -k_eff L, i.e. the residual is L A + k_eff L.
    L_naive = np.array([[1.0, 1.0]])
    k_eff = pk.CL / pk.Vd
    R_naive = L_naive @ A + k_eff * L_naive

    # discriminant: complex eigenvalues iff < 0
    tr = np.trace(A)
    det = np.linalg.det(A)
    disc = tr * tr - 4.0 * det

    return {
        "A": A,
        "eig": eig,
        "lam_slow": lam_slow,
        "lam_fast": lam_fast,
        "disc": disc,
        "l_slow": l_slow,
        "v_slow": v_slow,
        "L_exact": L_exact,
        "R_exact": R_exact,
        "L_naive": L_naive,
        "lam_naive": -k_eff,
        "R_naive": R_naive,
        "residual_exact": float(np.linalg.norm(R_exact)),
        "residual_naive": float(np.linalg.norm(R_naive)),
    }


def reduced_pk_profiles(
    pk: PKParams,
    dose_mg: float,
    t: np.ndarray,
    y0: Tuple[float, float] = (0.0, 0.0),
) -> Dict[str, np.ndarray]:
    """Compare full 2-compartment, exact slow-mode, and naive 1-compartment.

    The dose is a bolus into the central compartment at t=0 (IV-like), which
    is the cleanest setting for the lumping comparison.  The exact reduced
    model is exact on the slow manifold; the naive total-amount model is the
    standard 1-compartment approximation with Vd = Vc + Vp.
    """
    A = pk_matrix(pk)
    info = lumping_analysis(pk)
    l = info["l_slow"]
    v = info["v_slow"]
    lam = info["lam_slow"]

    # full model
    sol = solve_ivp(
        rhs_for_lumping, (t[0], t[-1]), np.array(y0, dtype=float),
        t_eval=t, args=(pk,), method="LSODA", rtol=1e-10, atol=1e-12,
    )
    A1_full, A2_full = sol.y[0], sol.y[1]

    # exact slow-mode reduction
    T0 = l @ np.array(y0, dtype=float)
    T_exact = T0 * np.exp(lam * (t - t[0]))
    A1_exact = v[0] * T_exact
    A2_exact = v[1] * T_exact

    # naive 1-compartment: T = A1 + A2, C_naive = T / Vd
    T0_naive = y0[0] + y0[1]
    k_eff = pk.CL / pk.Vd
    T_naive = T0_naive * np.exp(-k_eff * (t - t[0]))
    C_naive = T_naive / pk.Vd

    return {
        "t": t,
        "A1_full": A1_full, "A2_full": A2_full,
        "C_full": A1_full / pk.Vc,
        "A1_exact": A1_exact, "A2_exact": A2_exact,
        "C_exact": A1_exact / pk.Vc,
        "T_naive": T_naive, "C_naive": C_naive,
        "T_exact": T_exact,
        "lam_slow": lam,
        "lam_fast": info["lam_fast"],
        "l_slow": l,
        "v_slow": v,
    }


if __name__ == "__main__":
    from params import GEFITINIB, PD_GEFITINIB

    pk, pd = GEFITINIB, PD_GEFITINIB
    print("baseline PD:", baseline_pd(pd))
    print("Kd_eff:", kd_eff(pd))

    sim = simulate_regimen(pk, pd, n_days=14, t_step=0.25)
    t = sim["t"]
    print("t range:", t[0], t[-1], "n =", len(t))
    mask = t >= t[-1] - 24
    print("day-14 free Cmax (mg/L):", sim["C_free"][mask].max())
    print("day-14 free Cmax (uM):", sim["D_uM"][mask].max())
    print("day-14 occupancy max/min:",
          sim["occupancy"][mask].max(), sim["occupancy"][mask].min())
    print("day-14 pERK inhibition max:",
          sim["pERK_inhibition"][mask].max())
    print("day-14 apoptosis fold max:",
          sim["apoptosis_fold"][mask].max())

    info = lumping_analysis(pk)
    print("\nPK eigenvalues:", info["eig"])
    print("discriminant:", info["disc"])
    print("slow left eigenvector:", info["l_slow"])
    print("exact residual norm:", info["residual_exact"])
    print("naive residual norm:", info["residual_naive"])
    print("naive lambda:", info["lam_naive"])
