# Data dictionary

The scripts expect the following input files in this `data/` directory. The
individual-level data are **not** distributed with the code. Create these files
from your own dataset, matching the column names below. Time points are days
1, 2, 3, 5 and 7.

Measurement units: MAP and IAP in mmHg; fluid volume (VFR) in litres; age in
years.

---

## `dev_longitudinal_map_iap.csv`
Long-format MAP/IAP for the development cohort (used by **01**). Four columns,
in this order:

| column | description |
|--------|-------------|
| ID   | subject identifier |
| Time | day (1, 2, 3, 5, 7) |
| MAP  | mean arterial pressure |
| IAP  | intra-abdominal pressure |

---

## `dev_cohort_wide.xlsx`
One row per patient, development cohort (used by **02**). Required columns:

| column | description |
|--------|-------------|
| `IAP-MAPgroup` | trajectory group, integer 1–4 |
| `R28day`       | follow-up time to 28 days (days) |
| `R28Death`     | 28-day mortality event (0/1) |

---

## `dev_cohort_full.csv`
One row per patient, development cohort with baseline covariates and the raw
MAP/IAP measurements (used by **03** and **10**). Required columns:

| column | description |
|--------|-------------|
| `IAP-MAPgroup` | trajectory group, integer 1–4 |
| `ID.study.group.` | subject identifier |
| `R28day` | follow-up time to 28 days |
| `R28Death` | 28-day mortality event (0/1) |
| Age, SEX (0 = female, 1 = male), Weight | demographics |
| temperature01, pulse01, respirations01 | day-1 vital signs |
| `Etiology.of.pancreatitis` | aetiology, integer 1–4 (biliary/HTG/alcoholic/other) |
| WBC01, PLT01, BUN01, ALT01, PaO201 | day-1 laboratory values |
| APACHEII1 | day-1 APACHE II score |
| `IAP 01`, `IAP 02`, `IAP 03`, `IAP 05`, `IAP 07` | IAP by day |
| MeanArterialPressure01, meanarterialpressure02, MeanArterialPressure03, MeanArterialPressure05, MeanArterialPressure07 | MAP by day |

---

## `external_validation_315.xlsx`
External validation cohort, one row per patient (used by **03**). Required
columns:

| column | description |
|--------|-------------|
| trajectory_phenotype | assigned subphenotype label (Compensated / Hypodynamic / Decompensated / Pressure-compensated) |
| `28-hospitalday` | follow-up time to 28 days |
| `R-28Death` | 28-day mortality event (0/1) |
| age, gender (0/1), weight | demographics |
| temp_day1, pulse_day1, resp_day1 | day-1 vital signs |
| etiology | aetiology, integer 1–4 |
| WBC, PLT, BUN, ALT | laboratory values |

---

## `external_validation_merged.xlsx`
External cohort in wide MAP/IAP format for re-derivation of trajectories
(used by **05**). Required columns:

| column | description |
|--------|-------------|
| id | subject identifier |
| source | cohort/source label |
| trajectory_group | original assigned group |
| MAP_day1, MAP_day2, MAP_day3, MAP_day5, MAP_day7 | MAP by day |
| IAP_day1, IAP_day2, IAP_day3, IAP_day5, IAP_day7 | IAP by day |

---

## `phase2_lagged_primary.csv` and `phase2_lagged_complete_case.csv`
Long-format, one row per patient-day, for the fluid-response analyses (used by
**06/07/08/11**). The `primary` file carries imputed fluid values
(`VFR_L_imputed`); the `complete_case` file carries observed fluid values only.
Required columns:

| column | description |
|--------|-------------|
| ID | subject identifier |
| Day | day (1–5) |
| Subtype | subphenotype (Compensated / Hypodynamic / Decompensated / Pressure-Compensated) |
| VFR_L | observed fluid volume (L) |
| VFR_L_imputed | imputed fluid volume (L) — `primary` file only |
| dMAP, dIAP | next-day change in MAP / IAP |
| MAP, IAP | current-day MAP / IAP |
| NE, CRRT, PCD | treatment indicators (0/1): norepinephrine, CRRT, percutaneous catheter drainage |
| Age, SEX, APACHEII1, Etiology | baseline covariates |
