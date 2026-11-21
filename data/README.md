# External PK/PD dataset (LXF A677, gefitinib / erlotinib)

## Source

Eigenmann et al., *Combining Nonclinical Experiments with Translational PKPD
Modeling to Differentiate Erlotinib and Gefitinib*, **Mol Cancer Ther** 2016.

The data are redistributed in the GitHub repository
[`DavAug/ErlotinibGefitinib`](https://github.com/DavAug/ErlotinibGefitinib),
which re-analyses the published study. The repository is licensed under the
**BSD 3-Clause License**. The CSV files here are copied unmodified from the
repository's `data/` directory.

## Files

| File | Rows | Content |
|---|---:|---|
| `PKPD_ErloAndGefi_LXF.csv` | 1106 | Combined tumour-volume observations, dosing records, and individual PK parameters |
| `PK_LXF_gefi.csv` | 216 | Gefitinib plasma concentrations |
| `PK_LXF_erlo.csv` | 472 | Erlotinib plasma concentrations |

## Column dictionary (`PKPD_ErloAndGefi_LXF.csv`)

| Column | Meaning |
|---|---|
| `#ID` | Individual mouse / xenograft identifier |
| `TIME` | Time (days) |
| `DOSE` | Dose amount at a dosing record (μg per mouse) |
| `ADDL` | Number of additional doses (NONMEM-style) |
| `II` | Inter-dose interval (days) |
| `Y` | Observation: tumour volume (mm³) when `YTYPE == 2` |
| `YTYPE` | 2 = tumour-volume observation; missing = dosing record |
| `CELL LINE` | Cell-line / model code (1 = LXF A677) |
| `DOSE GROUP` | 0 (vehicle), 6.25, 25, 100 mg/kg |
| `DOSECAT` | Dose category |
| `DRUG` | 1 = erlotinib, 2 = gefitinib |
| `DRUGCAT` | Drug category (0 = vehicle) |
| `EXPERIMENT` | Experiment identifier |
| `BW` | Body weight (g) |
| `YTV` | Tumour-volume observation used in the model (subset of `Y`) |
| `KA`, `V`, `KE` | Individual PK parameters (absorption rate, volume, elimination rate) |
| `w0` | Initial tumour volume (mm³) |
| `I` | Individual model parameter from the source analysis |

## Structure

- Gefitinib: 4 groups (0, 6.25, 25, 100 mg/kg), 8 mice per group.
- Tumour volume: up to 10 time points per mouse (days 0–30).
- PK: 2 samples per mouse at steady state (≈ day 10 and ≈ day 16), 48 observations.
- Erlotinib: same design with its own dose groups.

## What this dataset can and cannot do

**Can:**

- Validate the dose–exposure relationship (linear vs nonlinear).
- Validate the steep dose–TGI threshold (100 mg/kg effective; 25 mg/kg not).
- Quantify the efficacy / body-weight (tolerability) trade-off.
- Provide individual-level variability for the tumour-growth layer.

**Cannot:**

- It is a mouse xenograft, not a human clinical PK/PD dataset.
- It reports tumour volume, not pERK or apoptosis; the current PD model needs a
  tumour-growth output layer to connect to it.
- The PK is sparse (two samples per mouse), so it is not suitable for a full
  individual PK fit.

## How to use it with the model

Add a tumour-volume output module

```
dV/dt = growth(V) - kill * f(occupancy or pERK)
```

calibrate `growth` on the vehicle group, calibrate `kill` on the 100 mg/kg group,
and use 6.25 / 25 mg/kg as external validation.

## Attribution

If you use these files, cite the original study:

> Eigenmann MJ, et al. Combining Nonclinical Experiments with Translational
> PKPD Modeling to Differentiate Erlotinib and Gefitinib. *Mol Cancer Ther.*
> 2016;15(12):3110–3121.

and the redistribution repository:

> Augustin D. *ErlotinibGefitinib* (GitHub repository), BSD-3-Clause.
