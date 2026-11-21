"""
pkpd-composition/params.py

Parameter sets for the compositional PK/PD case study.

Drug: gefitinib (EGFR tyrosine-kinase inhibitor), with erlotinib as a
secondary comparison.

The values are *representative* population values assembled from the
published clinical PK/PD literature and FDA labels.  They are not a
re-estimation of a specific dataset.  The point of the case study is the
compositional workflow and the interface/lumping analysis, not a new
population PK fit.  Every value is documented in README.md with its
provenance and, where it is a round number, labelled as representative.

Units
-----
PK:  mg, L, h; concentrations in mg/L (== ug/mL).
PD:  drug and receptor concentrations in uM; MEK/ERK in normalized units;
     apoptosis in normalized units.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PKParams:
    name: str
    MW: float          # g/mol
    fu: float          # free fraction
    F: float           # oral bioavailability
    ka: float          # absorption rate, 1/h
    Vc: float          # central volume, L
    Vp: float          # peripheral volume, L
    CL: float          # clearance, L/h
    Q: float           # inter-compartmental clearance, L/h
    dose_mg: float     # clinical dose, mg
    tau_h: float       # dosing interval, h

    @property
    def k10(self) -> float:
        return self.CL / self.Vc

    @property
    def k12(self) -> float:
        return self.Q / self.Vc

    @property
    def k21(self) -> float:
        return self.Q / self.Vp

    @property
    def Vd(self) -> float:
        return self.Vc + self.Vp


@dataclass(frozen=True)
class PDParams:
    """Reduced EGFR -> MEK -> ERK -> apoptosis module.

    R   free EGFR (uM)
    DR  drug-bound EGFR (uM)
    X   phospho-MEK (normalized)
    Y   phospho-ERK (normalized)
    A   apoptosis marker (normalized)
    """
    R_total: float     # total EGFR, uM
    Kd: float          # drug-EGFR dissociation constant, uM
    kon: float         # binding rate, 1/(uM h)
    koff: float        # unbinding rate, 1/h
    kint: float        # complex internalization, 1/h
    ksyn: float        # receptor synthesis, uM/h
    kdeg: float        # receptor degradation, 1/h

    KR: float          # half-saturation for MEK activation, uM
    nX: int            # Hill coefficient, MEK tier
    kx: float          # MEK activation rate, 1/h
    kdx: float         # MEK deactivation rate, 1/h

    KX: float          # half-saturation for ERK activation (normalized)
    nY: int            # Hill coefficient, ERK tier
    ky: float          # ERK activation rate, 1/h
    kdy: float         # ERK deactivation rate, 1/h

    KY: float          # half-saturation for apoptosis (normalized)
    nA: int            # Hill coefficient, apoptosis
    ka_apop: float     # apoptosis induction rate, 1/h
    kda: float         # apoptosis clearance rate, 1/h


# ---------------------------------------------------------------------------
# PK parameter sets
# ---------------------------------------------------------------------------

GEFITINIB = PKParams(
    name="gefitinib",
    MW=446.9,
    fu=0.10,       # ~90% protein bound
    F=0.60,        # oral bioavailability
    ka=0.50,       # 1/h, tmax ~3-5 h
    Vc=400.0,      # L
    Vp=600.0,      # L
    CL=15.0,       # L/h, CL/F ~25 L/h with F=0.6
    Q=50.0,        # L/h
    dose_mg=250.0, # standard daily dose
    tau_h=24.0,
)

ERLOTINIB = PKParams(
    name="erlotinib",
    MW=393.4,
    fu=0.07,       # ~93% protein bound
    F=0.60,
    ka=0.50,
    Vc=80.0,
    Vp=120.0,
    CL=4.0,
    Q=10.0,
    dose_mg=150.0,
    tau_h=24.0,
)


# ---------------------------------------------------------------------------
# PD parameter set (gefitinib, EGFR-mutant-like sensitivity)
# ---------------------------------------------------------------------------
#
# The binding module is a simple 1:1 drug-receptor model with kint = kdeg,
# so the total receptor pool is conserved at steady state and the occupancy
# follows the familiar D/(D+Kd) isotherm.  The MEK and ERK tiers are
# cooperative Hill stages; the apoptosis marker is induced when ERK is low.
#
# With these values the baseline is R=0.100 uM, X=0.800, Y=0.719, A=1.025,
# and the pERK IC50 is ~0.1 uM free drug.

PD_GEFITINIB = PDParams(
    R_total=0.100,
    Kd=0.050,
    kon=1.0,
    koff=0.050,
    kint=0.100,
    ksyn=0.010,
    kdeg=0.100,
    KR=0.050,
    nX=2,
    kx=1.0,
    kdx=1.0,
    KX=0.50,
    nY=2,
    ky=1.0,
    kdy=1.0,
    KY=0.50,
    nA=2,
    ka_apop=0.050,
    kda=0.020,
)
