# DLBCL gene-expression prognostic score — a pre-registered analysis of public data

**Interactive summary:** <https://gurbuzmericarda0-arch.github.io/dlbcl-prognostic-score/> (built from the summary results in this repository by `site/build_site.py`; research only, not a clinical tool).

> **Research project, not a clinical tool.** Nothing in this repository is intended for diagnosis, prognosis of
> individual patients or treatment decisions. The score is z-scored within each cohort and cannot be applied to a
> single new patient.

A 44-gene elastic-net Cox score for overall survival (OS) in diffuse large B-cell lymphoma (DLBCL), developed on
the training part of GSE181063 and evaluated in a held-out internal test set and two external cohorts on a
different platform. All data-dependent steps used training data only, and the analysis plan was hash-locked before
any model was fitted. The full report is in [`REPORT.md`](REPORT.md).

## Main results

| Set | Patients / deaths | HR per SD (95% CI) | Harrell's C (95% CI) |
|---|---|---|---|
| Training, nested CV (5×5) | 373 / 171 | — | 0.612 |
| Internal test (GSE181063) | 162 / 74 | 1.93 (1.49–2.49) | 0.664 (0.596–0.730) |
| External: GSE10846 | 229 / 58 | 1.69 (1.32–2.17) | 0.659 (0.592–0.732) |
| Supporting external: GSE87371 | 116 / 41 | 1.49 (1.09–2.05) | 0.630 (0.536–0.718) |

Negative controls: label-permutation leakage check mean C 0.499; 100 full-pipeline permutation models scored on the
test set, mean C 0.493 (none reached the observed C; permutation p = 0.010).

**Added value beyond IPI + cell of origin (COO)** — complete cases, LP per SD added to Cox(IPI + COO):

| | Internal test | GSE10846 | GSE87371 |
|---|---|---|---|
| Adjusted HR per SD | 1.63 (1.14–2.34) | 1.29 (0.93–1.80) | 1.17 (0.79–1.73) |
| Likelihood-ratio test | p = 0.007 | p = 0.134 | p = 0.434 |
| ΔHarrell's C (95% bootstrap CI) | +0.020 (−0.004 to +0.053) | +0.015 (−0.006 to +0.051) | +0.008 (−0.008 to +0.045) |

By the pre-specified rules the score is "internally validated", "externally validated" and "adds information beyond
IPI and COO". The last label rests on a weak external condition: **added value was shown statistically only in the
internal test set**; externally it was directionally consistent but not significant, and every ΔC interval includes 0.

## What this work does not show

- Discrimination is modest (C ≈ 0.63–0.66), similar to COO alone and below IPI + COO.
- Gene selection is unstable (4 of 44 genes selected in ≥ 50% of bootstraps); the validated object is the 44-gene
  score as a whole, not individual genes.
- The score mostly re-captures known biology (proliferation/MYC, ABC/IRF4, stromal and T-cell programmes).
- No external calibration or clinical-utility analysis is possible because of within-cohort z-scoring.
- Author-processed GEO values were used (CEL files not re-processed); IPI was missing for some patients and analysed
  as complete cases.

## Pre-registration: what can and cannot be verified

Plan v1.0 (`results/04/analysis_plan_04.md` and `.json`) was locked in `results/04/plan_lock_04.json`
(SHA-256 `8d87cf887e2db4770ae3273755d4901cd0d5a2a010fdc729159e007aa129bc7d`) before model fitting; the frozen model `results/04/model_04.json` has
SHA-256 `661e2b80be774252a7fd613116ecc06673a8617ed8a2230bdae56170e3908403`. Every later notebook re-verifies these hashes and records what it checked in
its manifest. `results/04/plan_amendments.md` records no amendments. **Limitation:** the lock was kept privately on
Google Drive and was not deposited in an external registry or timestamping service before the test and external
data were opened, so the order of events is documented but cannot be independently proven.

## Data

All data are public and de-identified; **raw data are not redistributed** — the notebooks download them from GEO.

| Cohort | Role | Platform | Source study |
|---|---|---|---|
| [GSE181063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181063) | Development + internal test | GPL14951 (Illumina DASL, FFPE) | Lacy et al., Blood 2020 (PMID 32187361); Painter et al., Br J Haematol 2019 (PMID 30408148) |
| [GSE10846](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846) | External validation | GPL570 | Lenz et al., N Engl J Med 2008 |
| [GSE87371](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87371) | Supporting external validation | GPL570 | Dubois et al., Clin Cancer Res 2017 (PMID 27923841); Dubois et al., EBioMedicine 2019 (PMID 31648986) |

Patient-level score files (`results/05/test_lp_05.csv`, `results/06/external_lp_06.csv`) contain GEO sample accessions and derived scores only; no clinical data are redistributed. Gene sets for step 8: MSigDB v2026.1.Hs (not redistributed).

## Repository layout

```
notebooks/        Colab notebooks, steps 00-09b (outputs cleared)
results/02/       step-2 manifests, cohort summary, patient flow, common gene table
results/03-08/    summary results and manifest of each step (SHA-256 of inputs/outputs)
results/04/       locked analysis plan, lock file, amendment log, frozen model
results/09b/      Supplementary Tables S1-S2, release manifest
figures/          all figures; figures/09b/flow_diagram_09b.png = Figure 1
REPORT.md         full report (TRIPOD+AI and REMARK checklists in the appendices)
SHA256SUMS.txt    checksums of every file in this release
```

## Reproducing

1. Google Colab (free tier; 2 CPU was used) with a Drive folder `MyDrive/dlbcl-prognosis`.
2. Run the notebooks in order: 02a → 02b → 03 → 04a → 04b → 05 → 06 → 07 → 08 (→ 09b). 02a/02b download the GEO
   series matrices and platform annotations (the GPL14951 SOFT file is about 2.4 GB).
3. Each notebook checks the SHA-256 of its inputs against the lock and earlier manifests and stops on a mismatch.
   Long steps: 05 ≈ 85 min (100 full-pipeline permutations).
4. Environment of the recorded runs: Python 3.13.16, scikit-survival 0.28.0, scikit-learn 1.9.1, lifelines 0.30.3.

## Licence and citation

Code: MIT ([`LICENSE`](LICENSE)). Text, figures and summary results: CC BY 4.0. GEO data remain under the terms of
their original submitters. Citation metadata: [`CITATION.cff`](CITATION.cff).

## Use of AI

Analysis code and report drafting were assisted by an AI model (Claude, Anthropic). All numbers come from the
project's run logs and manifests; design decisions and their limitations are documented in `REPORT.md`.
