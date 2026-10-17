# DLBCL Gene-Expression Prognostic Score — Report

Version 1.0 · 2026-10-17 · Meriç Arda Gürbüz · Repository: <https://github.com/gurbuzmericarda0-arch/dlbcl-prognostic-score>

**Working title:** Development and validation of a penalised Cox gene-expression score for overall survival in diffuse large B-cell lymphoma, and its added value beyond IPI and cell of origin: a pre-registered analysis of public data

**Study type:** Prognostic model development with internal (split-sample) and external validation (TRIPOD+AI); prognostic marker study (REMARK). Research study only; not a clinical tool.

## Abstract

**Background.** Gene-expression profiling separates diffuse large B-cell lymphoma (DLBCL) into prognostic groups, but published survival signatures often fail to replicate and rarely show added value beyond the International Prognostic Index (IPI) and cell of origin (COO). We asked whether an overall-survival (OS) score built with a locked, leakage-controlled pipeline on public data validates across platforms and adds information to IPI + COO.

**Methods.** Three public GEO cohorts of R-CHOP-type-treated DLBCL were used: GSE181063 (Illumina DASL, FFPE; development and internal test), GSE10846 (Affymetrix GPL570; primary external) and GSE87371 (GPL570; supporting external). GSE181063 was split 70/30, stratified by OS event (training 380 patients, 173 deaths; test 164, 74 deaths). An analysis plan with three hypotheses and fixed conclusion sentences was written and hash-locked before model fitting. In training only, an elastic-net Cox model (l1\_ratio 0.9) was fitted on the 5,000 most variable of 13,747 cross-platform genes, with penalty chosen by 10-fold cross-validation (λ\_min) and performance estimated by 5×5 nested cross-validation. Positive control: ABC vs GCB survival difference. Negative controls: a label-permutation leakage check and 100 full-pipeline permutations. H1: hazard ratio (HR) per SD > 1 in the internal test set; H2: the same in external cohorts with frozen genes and coefficients; H3: likelihood-ratio test (LRT) of the score added to IPI + COO in the test set.

**Results.** The positive control was reproduced (training HR ABC/GCB 2.21, 95% CI 1.57–3.10). The frozen model contained 44 genes; nested-CV Harrell's C was 0.612. In the test set (162 patients, 74 deaths) the HR per SD was 1.93 (1.49–2.49; H1 met), C 0.664 (0.596–0.730); permutation p = 0.010. Externally, HR per SD was 1.69 (1.32–2.17) in GSE10846 (229 patients, 58 deaths; H2 met) and 1.49 (1.09–2.05) in GSE87371 (116 patients, 41 deaths). Adjusted for IPI + COO, the test-set HR was 1.63 (1.14–2.34; LRT p = 0.007; H3 met), but external adjusted HRs were 1.29 (0.93–1.80) and 1.17 (0.79–1.73), with non-significant LRTs; ΔC was +0.008 to +0.020 with every 95% CI including 0. Gene selection was unstable (4/44 genes selected in ≥50% of bootstraps). The score correlated with COO (ρ 0.46) and IPI (ρ 0.36) and mainly re-captured known proliferation/MYC, ABC/IRF4, stromal and T-cell programmes.

**Conclusions.** By the pre-specified criteria the score is internally validated, externally validated and adds information beyond IPI and COO. However, added value was shown statistically only in the internal test set; in external cohorts it was directionally consistent but not demonstrated, and gains in discrimination were small and uncertain. Because the score is z-scored within each cohort, it cannot be applied to a single new patient. This is a research study, not a clinical tool.

## 1 Introduction

Diffuse large B-cell lymphoma (DLBCL) is clinically and biologically heterogeneous: most patients treated with R-CHOP-type immunochemotherapy are cured, while a substantial minority die of their disease. Two established prognostic tools frame this heterogeneity. The International Prognostic Index (IPI) combines age, stage, LDH, performance status and extranodal involvement \[1\]. Gene-expression profiling defined cell-of-origin (COO) subgroups, with activated B-cell-like (ABC) DLBCL having worse outcome than germinal-centre B-cell-like (GCB) DLBCL \[2,3\]; COO can be assigned on FFPE material with assays such as Lymph2Cx \[4\]. Gene-expression studies also identified proliferation, MYC-related and stromal/microenvironment programmes associated with survival \[3,5\], and a molecular high-grade (MHG) group with poor outcome \[6\].

Many gene-expression survival signatures for DLBCL have been published using public GEO cohorts. Common weaknesses are gene selection on the full dataset before splitting, data-driven cut-offs, absent or same-platform-only validation, and no test of added value beyond IPI and COO. These practices inflate apparent performance and make replication unlikely.

**Objectives.** Using only public, de-identified data, we aimed to (i) reproduce a known prognostic finding (ABC vs GCB) as a positive control; (ii) develop a penalised Cox gene-expression score for overall survival (OS) in a training set only; (iii) validate it in a held-out internal test set and in two external cohorts on a different platform; (iv) test whether it adds prognostic information beyond IPI and COO; and (v) describe its biology. The study combines model development and validation (TRIPOD+AI type: development with split-sample internal validation plus external validation) \[7\] and is reported also against REMARK \[8,9\].

**Intended use.** None in patient care. The score is a research instrument for studying which biological programmes carry prognostic information in DLBCL. It is not intended for diagnosis, prognosis of individual patients or treatment decisions.

## 2 Methods

### 2.1 Data sources and design

All data are public and de-identified, downloaded from NCBI GEO \[12\] as series-matrix files with author-processed expression values. Raw CEL files were not re-processed.

| Cohort | Role | Platform | Tissue | Primary reference |
| --- | --- | --- | --- | --- |
| GSE181063 | Development (training) and internal test | GPL14951, Illumina HumanHT-12 WG-DASL V4.0 | FFPE | HMRN population-based cohort \[10\] |
| GSE10846 | Primary external validation | GPL570, Affymetrix HG-U133 Plus 2.0 | Frozen | Lenz et al. \[5\] |
| GSE87371 | Supporting external validation | GPL570 | Frozen | LYSA cohort, Dubois et al. \[11\] |

The study had seven sequential steps, each in its own Colab notebook with a manifest of SHA-256 hashes of inputs and outputs: clinical cleaning (2a), expression processing (2b), positive control (3), analysis-plan lock (4a), model development (4b), internal validation (5), external validation (6), clinical comparison (7) and exploratory biology (8). The analysis plan (version 1.0) was written and hash-locked (lock file SHA-256 `8d87cf887e2d…`) before any model was fitted; no amendments were made (`plan_amendments.md` empty). The internal test set was opened once, after the model was frozen (model SHA-256 `661e2b80be77…`); external expression files were opened once, after internal validation.

### 2.2 Participants and eligibility

- **GSE181063:** diagnosis DLBCL NOS; `qc_fail == 0`; first-line R-CHOP with curative intent; OS time > 0; pre-treatment biopsy only (`EXCLUDE_NON_PRETREATMENT = True`). 633 arrays corresponded to 596 patients; duplicates (33 patients) were resolved by a clinical fingerprint (age, sex, OS time, OS status), with 0 conflicts against published patient IDs; one array per patient was kept (pre-treatment, highest probe detection). Final: 544 patients, 247 deaths, median follow-up 97.5 months.
- **GSE10846:** reference CD19 samples removed; R-CHOP-like treatment; OS in years converted to months; LDH ratio 0 set to missing. Final: 232 patients, 59 deaths, median follow-up 33.7 months.
- **GSE87371:** R-CHOP14/21 or R-mini-CHOP21 (13 patients, sensitivity flag); PMBL excluded. Final: 116 patients, 41 deaths, median follow-up 49.2 months.

Sample-level QC (2b) flagged correlation outliers (robust z < −4) and sex mismatches between clinical sex and expression (GSE181063: four Y-chromosome genes, XIST absent on GPL14951; GPL570: Y genes − XIST). The sex rule (threshold midway between clinical male and female medians; ambiguous band ± 25% of the distance) was set after the first QC run, using training/reference samples only and no survival data. Per the locked plan, `sex_mismatch_confident` samples were excluded in all sets; correlation outliers were kept in primary analyses and removed in sensitivity analyses.

### 2.3 Outcome

The primary outcome was overall survival (OS), death from any cause. Event coding was verified from the data: GSE87371 OS event = 1 − `cens_os` (0 = death), confirmed by censored PFS with no death and PFS ≤ OS; GSE181063 `os_status` 1 = death. PFS was a secondary endpoint, not used for model training; it is not reported here. No outcome assessment was blinded (retrospective public data).

### 2.4 Expression processing (training-only references)

All three sets were already log2-scaled (p99 ≤ 50) and had no missing values. Each sample was quantile-mapped to a reference distribution; for GSE181063 the reference used training samples only. Probes were mapped to genes by Entrez ID (coverage 99.7% GPL14951, 100% GPL570); multi-gene probes were dropped; per gene, the probe with the highest mean was kept (GSE181063: training mean). Probe annotation: NCBI `GPL570.annot.gz`; for GPL14951 no `.annot` file existed, so `GPL14951_family.soft.gz` was used. Genes in the bottom 25% of GSE181063 training means (threshold 8.376) were removed (15,231 genes kept). Intersection with GPL570 gave **13,747 common genes**. Survival data were not read in this step.

### 2.5 Training/test split

GSE181063 was split 70/30 stratified by OS event (seed 20261014): training 380 patients (173 deaths), test 164 (74 deaths). After QC exclusions the primary training population was **373 patients, 171 deaths** and the test population 162 patients, 74 deaths.

### 2.6 Positive control (step 3)

In the training set only, the hypothesis that ABC has worse OS than GCB was written to `plan_03.json` before analysis (success: HR > 1, 95% CI excluding 1, log-rank p < 0.05). Kaplan–Meier, log-rank and univariable Cox (ABC vs GCB) were used, with IPI-adjusted and QC-exclusion sensitivity analyses. Label consistency was checked with a simple Lymph2Cx gene score \[4\] (15/15 genes present), without survival data.

### 2.7 Predictors and model development (step 4b)

Candidate predictors were the 13,747 genes only; clinical variables were not model inputs, and there was no survival-based pre-screening. Within every fit (including every CV fold), the 5,000 most variable genes on that fit's training rows were selected and z-scored (ddof 1). The model was an elastic-net Cox model (scikit-survival `CoxnetSurvivalAnalysis` \[15\]; l1\_ratio 0.9; 100 penalty values; alpha\_min\_ratio 0.01; Breslow ties \[31\]) \[13,30\]. The penalty was chosen by 10-fold CV stratified by event, using the grouped cross-validated partial-likelihood deviance of Verweij and van Houwelingen as implemented in glmnet \[14,13\]. λ\_min was primary; λ\_1se was a sensitivity analysis. Pre-specified rules: a model with no genes would be reported as a null result; λ\_min at the grid edge would be reported without widening the grid. A safeguard against a scikit-survival early-stopping issue with tied times was added; it never fired on real data.

The risk score was the linear predictor LP = Σ coef · (x − train\_mean)/train\_sd over the selected genes. Other training-only outputs: 5×5 nested CV Harrell's C \[17\] (λ\_min and λ\_1se); a cross-validated KM curve (repeat 1); descriptive added value of the out-of-fold LP over IPI + COO; gene-selection stability (100 bootstraps); and sensitivity analyses (λ\_1se, pure LASSO, no variance filter, correlation outliers excluded, five CV seeds). None of these were carried to the test set.

### 2.8 Internal validation and negative controls (step 5)

Before opening the test set: (i) the frozen model was re-derived from training data (identical genes, λ and coefficients); (ii) a **leakage check** ran nested CV (5×1) on 10 label permutations of the training data, with a stop rule if mean C > 0.55; (iii) 100 **negative-control** models were built, each re-running the whole pipeline (including CV) on permuted training labels (time and event permuted together). Test expression was scaled with training means and SDs. Pre-specified measures: **H1** HR per SD of LP > 1 with 95% CI excluding 1 (single primary hypothesis); Harrell's C with 1,000-bootstrap CI; Uno's C \[18\] at τ = 60 months (shortened if < 20% at risk); KM and log-rank using the training-median cut-off (−0.0048; LP > cut-off = high risk); calibration slope; proportional-hazards checks \[20\]. The 100 permutation models were scored on the real test labels to give a permutation p-value.

### 2.9 External validation (step 6)

The 44 frozen genes and coefficients were applied unchanged. Because platforms differ, each gene was z-scored within each external cohort's primary population (ddof 1); training means/SDs, baseline hazard and the numeric cut-off were not transported. Risk groups used each cohort's own median. Calibration was not assessed (plan §10.4). **H2:** HR per SD > 1 with CI excluding 1 in GSE10846; "externally validated" also required HR > 1 in GSE87371. Sensitivity: GSE87371 without R-mini-CHOP21 (13 patients); correlation outliers excluded; no re-scaling in either.

### 2.10 Comparison with clinical predictors (step 7)

In each validation set (complete cases for IPI and COO), seven unpenalised Cox models (lifelines \[16\], Efron ties \[32\]) were fitted: IPI (0–5, linear), COO (GCB reference; ABC; other = UNC/MHG/unclassified), LP (per SD), IPI + COO, IPI + LP, COO + LP and IPI + COO + LP. **H3:** LRT (1 df) of Cox(IPI + COO) vs Cox(IPI + COO + LP) in the test set, p < 0.05. Also reported: adjusted HR per SD, apparent ΔHarrell's C with 1,000 bootstraps (both models refitted per resample; seed 20261515; percentile CI), missingness diagnostics (event rate and LP by IPI availability) and a sensitivity analysis with MHG as a separate category. "Adds information beyond IPI and COO" required H3 in the test set and adjusted HR > 1 (point estimate) in both external cohorts.

### 2.11 Exploratory biology (step 8)

In the training population only, without survival data: over-representation analysis (one-sided hypergeometric, BH-FDR \[23\]) of all 44 genes and of risk-increasing (26) and risk-decreasing (18) genes, with the 5,000-gene variance-filtered pool as background (13,747 genes as sensitivity); a score-based gene-set test using Spearman ρ of every gene with the LP (model genes removed; 4,956 genes), assessed both by LP permutation (2,000) and a PAGE-type test \[22\], with "robust" sets requiring FDR < 0.05 in both; and correlations of the LP with signature scores (mean z-score) — Lymph2Cx COO \[4\], Hallmark proliferation/MYC sets, Hallmark EMT and angiogenesis as stromal proxies (not the Lenz 2008 gene lists), and 30 pre-named DLBCL/lymphoma C2 CGP sets — from MSigDB v2026.1.Hs \[21\]. Settings were saved in `analysis_config_08.json` before results were seen.

### 2.12 Sample size, multiplicity and software

No formal sample-size calculation for development was made; the cohort size was fixed by data availability (training 171 events for up to 13,747 candidate genes, handled by penalisation). Schoenfeld power \[19\] (80%, per SD) gave minimum detectable HRs of ≈ 1.38 in the test set (74 deaths), ≈ 1.44 in GSE10846 (59) and ≈ 1.55 in GSE87371 (41). H1 was the only primary hypothesis; H2 and H3 were secondary and not adjusted for multiplicity (plan §15). Master seed 20261015. Software: Google Colab (2 CPU), Python 3.13.16, scikit-survival 0.28.0, scikit-learn 1.9.1 \[33\], lifelines 0.30.3. Notebooks were tested on synthetic data before each real run. Red flags pre-specified: any C > 0.90; permutation-model mean C > 0.55; validation C exceeding nested-CV mean by > 0.10.

### 2.13 Ethics

Only publicly available, de-identified data were used; no new ethical approval was sought. There was no patient or public involvement. No fairness analysis was planned (race/ethnicity were not available).

## 3 Results

### 3.1 Cohorts and patient flow

**Table 1. Cohorts after eligibility filtering (before sample QC).**

| Cohort | Platform | Patients | Deaths | Median follow-up (months) | IPI missing |
| --- | --- | --- | --- | --- | --- |
| GSE181063 | GPL14951 (DASL, FFPE) | 544 | 247 | 97.5 | 129 |
| GSE10846 | GPL570 | 232 | 59 | 33.7 | 76 |
| GSE87371 | GPL570 | 116 | 41 | 49.2 | 0 |

GSE181063 was split into training (380 patients, 173 deaths) and test (164, 74). Sample QC flags are in Table 2. After excluding confident sex mismatches, analysis populations were: training 373/171, test 162/74, GSE10846 229/58, GSE87371 116/41 (patients/deaths). Median follow-up (reverse KM) was 101.3 months (test), 34.2 (GSE10846) and 49.2 (GSE87371). In the training population, IPI was missing in 92/373; COO was GCB 185, ABC 104, UNC 61, MHG 23. 13,747 genes were common to both platforms, including all 44 model genes (matched by Entrez ID).

**Table 2. Sample QC flags.**

| Cohort | Set | Correlation outliers | Sex mismatch (all) | Confident sex mismatch |
| --- | --- | --- | --- | --- |
| GSE181063 | Training (380) | 19 | 18 | 7 |
| GSE181063 | Test (164) | 6 | 9 | 2 |
| GSE10846 | External (232) | 8 | 3 | 3 |
| GSE87371 | External (116) | 2 | 0 | 0 |

No batch clustering was visible on PCA (PC1 12.1%, PC2 8.8%).

### 3.2 Positive control (training set)

The pre-specified success criterion was met. ABC had worse OS than GCB: log-rank χ² = 21.94, p = 2.8 × 10⁻⁶; univariable HR 2.21 (1.57–3.10; n = 295, 134 deaths), C 0.599; proportional-hazards test p = 0.744.

**Table 3. OS by COO in the training set.**

| Group | Patients | Deaths | Median OS, months (95% CI) | 5-year OS (95% CI) |
| --- | --- | --- | --- | --- |
| ABC | 105 | 67 | 64.8 (36.6–91.8) | 51.4% (41.5–60.5) |
| GCB | 190 | 67 | Not reached | 77.3% (70.7–82.7) |
| MHG | 24 | 13 | 50.9 (9.5–NA) | 49.4% (28.3–67.4) |
| UNC | 61 | 26 | 124.9 (77.5–NA) | 70.2% (57.0–80.1) |

Four-group log-rank (descriptive) p = 2.2 × 10⁻⁵. Restricted to IPI-complete cases (n = 223, 92 deaths) the unadjusted HR was 2.50 (1.66–3.78), C 0.614, and the IPI-adjusted HR 2.67 (1.76–4.03), C 0.736 (IPI per point HR 1.77, 1.47–2.13). Excluding confident sex mismatches (n = 289) gave HR 2.25 (1.60–3.16); excluding correlation outliers (n = 282) gave 2.21 (1.57–3.12). The Lymph2Cx gene score separated the authors' ABC and GCB labels with AUC 0.952.

### 3.3 Model development (training set)

λ\_min = 0.08626 (grid position 18/100, not at the edge) selected **44 genes** (26 risk-increasing, 18 risk-decreasing); λ\_1se (position 5/100) selected 3. The largest absolute coefficients (log HR per SD; negative = higher expression, lower risk) were JAML −0.164, PXDNL −0.137, PLAC8L1 +0.080, FRRS1 +0.071, ZBTB7C −0.069, PALD1 −0.069, CXCL2 +0.068 and GPR141 −0.066. No path events, numerical errors or convergence warnings occurred in the final fit, nested CV or bootstraps. Training-median cut-off: −0.0048.

Nested CV (5×5) gave Harrell's C **0.612** (SD 0.051; range 0.539–0.726; median 19 genes; 0/25 null models) for λ\_min and 0.539 (SD 0.049; median 2 genes; 12/25 null models) for λ\_1se. In the cross-validated KM (repeat 1, descriptive), the high-risk group (195 patients, 114 deaths) had 5-year OS 54.2% (46.9–60.9) and median 74.4 months; the low-risk group (178, 57) had 80.9% (74.3–85.9), median not reached.

**Descriptive added value in training** (281 patients with IPI, 118 deaths): out-of-fold LP HR per SD 1.22 (0.98–1.51), LRT χ² 3.23, p = 0.072; C 0.738 → 0.742 (ΔC +0.003).

**Gene-selection stability** (100 bootstraps): only 4 of 44 genes were selected in ≥ 50% of resamples — PXDNL 84%, JAML 68%, PLAC8L1 60%, AGBL4 50% (next: ATAD3B 49%, CXCL2 44%). Bootstrap model size: median 58.5 genes (IQR 54–63).

**Table 4. Training-only sensitivity analyses.**

| Analysis | Genes | Shared with primary | Jaccard | Nested-CV C |
| --- | --- | --- | --- | --- |
| λ\_1se | 3 | 3 | 0.07 | 0.539 (5×5) |
| Pure LASSO | 44 | 44 | 1.00 | 0.614 (5×1) |
| No variance filter (13,747 genes) | 14 | 9 | 0.18 | 0.606 (5×1) |
| Correlation outliers excluded (n = 356) | 20 | 16 | 0.33 | — |
| Five other CV seeds | 44 / 11 / 35 / 40 / 40 | all | 1.00 / 0.25 / 0.80 / 0.91 / 0.91 | — |

Smaller models from other CV seeds were subsets of the 44-gene model.

### 3.4 Internal validation and negative controls (test set)

The frozen model was reproduced exactly from training data (44 genes, λ = 0.086259, coefficient difference 0). The leakage check passed: mean C 0.499 (folds 0.422–0.534; 70% null models). Across the 100 negative-control models, test C averaged 0.493 (2.5–97.5th percentile 0.395–0.559; maximum 0.580); none reached the observed C (permutation p = 0.010, the smallest possible with 100 permutations being 0.0099).

**Table 5. Internal validation (test set, 162 patients, 74 deaths).**

| Measure | Value |
| --- | --- |
| **H1: HR per SD** | **1.93 (1.49–2.49), p = 5.2 × 10⁻⁷ — met** |
| Harrell's C (1,000 bootstraps) | 0.664 (0.596–0.730) |
| Uno's C (τ = 60 months) | 0.670 |
| Log-rank (training cut-off) | p = 2.4 × 10⁻⁴ |
| HR high vs low | 2.54 (1.52–4.24) |
| Calibration slope | 1.56 (0.95–2.16) |
| Proportional hazards | p = 0.218 (rank), 0.267 (KM) |

With the training cut-off, 90 patients were high risk (54 deaths; median OS 82.1 months, 50.1–132.9; 5-year OS 56.7%, 45.8–66.2) and 72 low risk (20 deaths; median not reached; 5-year OS 81.9%, 70.8–89.0). Tertiles were monotone (descriptive): T1 39 patients, 6 deaths, 5-year OS 86.9%; T2 73, 34, 68.5%; T3 50, 34, 52.0%. Sensitivity: within-test z-scoring HR 1.91 (1.48–2.47), C 0.662; correlation outliers excluded (n = 156) HR 1.97 (1.51–2.58), C 0.668. Test C minus nested-CV C = +0.052; no red flag. **Per plan §12: "internally validated".**

### 3.5 External validation

**Table 6. External validation (frozen genes and coefficients, within-cohort z-scores).**

| Measure | GSE10846 (primary) | GSE87371 (supporting) |
| --- | --- | --- |
| Patients / deaths (after QC) | 229 / 58 (3 excluded) | 116 / 41 (0 excluded) |
| **H2: HR per SD** | **1.69 (1.32–2.17), p = 3.7 × 10⁻⁵ — met** | **1.49 (1.09–2.05), p = 0.013 — met** |
| Harrell's C (1,000 bootstraps) | 0.659 (0.592–0.732) | 0.630 (0.536–0.718) |
| Uno's C | 0.645 (τ = 47.9 months) | 0.626 (τ = 55.2 months) |
| Log-rank (own median) | p = 9.2 × 10⁻⁴ | p = 0.016 |
| HR high vs low | 2.46 (1.42–4.25) | 2.17 (1.14–4.14) |
| Proportional hazards (rank) | p = 0.946 | p = 0.441 |

In both cohorts fewer than 20% of patients were at risk at 60 months, so τ was shortened by the plan rule. Five-year OS: GSE10846 high 61.9% (51.5–70.7) vs low 75.5% (61.2–85.1); GSE87371 high 45.6% (30.4–59.6) vs low 75.2% (58.8–85.7); medians were 50.3 vs 68.6 months in GSE87371 and not reached in GSE10846. The GSE87371 low-risk median rests on a single late death when only two low-risk patients remained at risk, and the two curves cross after about 68 months; this median should not be interpreted. Tertiles were monotone in both (GSE10846 deaths 10/18/30, 5-year OS 82.7/67.3/56.6%; GSE87371 deaths 9/13/19, 76.8/60.2/45.0%).

Sensitivity analyses: GSE87371 without R-mini-CHOP21 (13 patients, 5 deaths removed; n = 103, 36 deaths) HR 1.39 (0.98–1.96), C 0.606, log-rank p = 0.078 — the CI now includes 1; GSE10846 without correlation outliers (n = 221, 55 deaths) HR 1.86 (1.39–2.49), C 0.667; GSE87371 without outliers (n = 114, 40 deaths) HR 1.46 (1.06–2.02), C 0.623. C minus nested-CV mean was +0.047 (GSE10846) and +0.018 (GSE87371); no red flag. **Per plan §12: "externally validated".**

### 3.6 Comparison with IPI and COO

Complete cases were 126/162 patients (54 deaths) in the test set, 154/229 (39) in GSE10846 and 116/116 (41) in GSE87371. COO was never missing. In the test set, patients with missing IPI had a higher death proportion (55.6% vs 42.9%; log-rank p = 0.074); in GSE10846 it was the same (25.3% vs 25.3%; p = 0.541). LP distributions did not differ by IPI availability (Mann–Whitney p = 0.309 and 0.522).

**Table 7. Harrell's C (complete cases, apparent) and ΔC for adding LP to IPI + COO.**

| Model | Test | GSE10846 | GSE87371 |
| --- | --- | --- | --- |
| IPI | 0.679 | 0.706 | 0.731 |
| COO | 0.647 | 0.695 | 0.641 |
| LP | 0.688 | 0.722 | 0.630 |
| IPI + COO | 0.727 | 0.765 | 0.752 |
| IPI + COO + LP | 0.746 | 0.780 | 0.760 |
| **ΔC (95% bootstrap CI)** | **+0.020 (−0.004 to +0.053)** | **+0.015 (−0.006 to +0.051)** | **+0.008 (−0.008 to +0.045)** |

**Table 8. LP hazard ratio per SD and likelihood-ratio tests.**

|  | Test | GSE10846 | GSE87371 |
| --- | --- | --- | --- |
| Unadjusted HR | 2.18 (1.59–3.00) | 1.88 (1.41–2.51) | 1.49 (1.09–2.05) |
| Adjusted for IPI | 1.93 (1.39–2.67) | 1.55 (1.15–2.09) | 1.35 (0.96–1.89) |
| Adjusted for COO | 1.94 (1.38–2.73) | 1.52 (1.09–2.12) | 1.26 (0.88–1.80) |
| **Adjusted for IPI + COO** | **1.63 (1.14–2.34)** | **1.29 (0.93–1.80)** | **1.17 (0.79–1.73)** |
| **LRT, LP added to IPI + COO** | **χ² 7.35, p = 0.007** | χ² 2.25, p = 0.134 | χ² 0.61, p = 0.434 |
| LRT, LP added to IPI | p = 4.4 × 10⁻⁵ | p = 0.004 | p = 0.079 |
| LRT, LP added to COO | p = 1.0 × 10⁻⁴ | p = 0.015 | p = 0.210 |

In the full model, HRs for IPI (per point) and ABC were 1.51 and 2.51 (test), 1.63 and 3.01 (GSE10846) and 1.98 and 1.62 (GSE87371). Spearman correlation of LP with IPI was +0.40 (test), +0.38 (GSE10846) and +0.28 (GSE87371). Treating MHG as a separate category (test, n = 7 MHG) gave LRT p = 0.007, adjusted HR 1.65 (1.13–2.39), ΔC +0.017 (−0.005 to +0.051). All 1,000 bootstrap resamples were valid in every analysis.

**Per plan §12: "adds information beyond IPI and COO"** (H3 met in the test set; adjusted HR point estimates 1.29 and 1.17 > 1 in the external cohorts). This label rests on a weak external condition: in both external cohorts the adjusted HR CIs include 1 and the LRTs are non-significant, and in all three cohorts the ΔC CI includes 0. Added value beyond IPI + COO was therefore **demonstrated statistically only in the internal test set**.

### 3.7 Exploratory biology (training set)

ORA found no gene set with FDR < 0.05 against either background, as expected for 44 unstable genes. The score-based test found robust sets (FDR < 0.05 in both tests) in Hallmark 11/50, Reactome 198/580 and GO BP 306/2,802; because sets overlap heavily, themes rather than individual sets are interpreted.

- **Higher with higher risk score:** E2F targets (z +12.3), G2M checkpoint (+11.4), MYC targets V1 (+8.1) and V2 (+5.6), mTORC1 (+6.7), unfolded protein response (+3.2); Reactome cell cycle, DNA replication and homologous-recombination repair. Leading genes: CDC25A, MAD2L1, TRIP13, MYC, EXO1, POLQ, CENPA (ρ ≈ +0.32 to +0.38).
- **Higher with lower risk score:** extracellular matrix/stroma — Hallmark EMT (z −5.5), myogenesis (−5.7), apical junction, apical surface, coagulation; Reactome ECM organisation, elastic fibres, collagen formation, proteoglycans (HTRA1, TGFB3, MFAP4, LAMB2, SPARC, LUM); and T-cell genes (CD2, LCK, ICOS, CD6, TRAT1, PRKCQ, SLAMF1).

**Table 9. Correlation of LP with signature scores (Spearman ρ, training, n = 373 unless stated).**

| Signature | ρ (95% CI) |
| --- | --- |
| Lymph2Cx COO score | +0.46 (0.37–0.55) |
| Hallmark MYC targets | +0.42 |
| Hallmark G2M | +0.41 |
| Hallmark E2F | +0.39 |
| MORI Eμ-MYC lymphoma (C2 CGP) | +0.39 |
| SHIPP\_DLBCL\_CURED\_VS\_FATAL\_DN (high in fatal disease) | +0.37 (+0.31 without model genes) |
| IPI (n = 281) | +0.36 |
| DAWSON methylation set (C2 CGP) | −0.36 |
| SHAFFER IRF4 targets (C2 CGP) | +0.31 to +0.32 |
| Hallmark EMT (stromal-1 proxy) | −0.23 |
| SHIPP\_DLBCL\_CURED\_VS\_FATAL\_UP (high in cured) | −0.20 |
| Hallmark angiogenesis (stromal-2 proxy) | −0.20 |
| IFN-γ, inflammatory, TNFα/NF-κB, IL6/JAK/STAT3 | \|ρ\| ≤ 0.11 |

The Lymph2Cx score separated the authors' ABC/GCB labels with AUC 0.959, and the LP differed across COO groups (Kruskal–Wallis p = 5.0 × 10⁻¹⁹). Of 30 pre-named C2 CGP sets, 17 had FDR < 0.05. Lymph2Cx alone explained 18.6% of LP variance and six signatures together 35.2% (descriptive, model genes removed); about 65% was not explained by these signatures.

## 4 Discussion

### 4.1 Main findings

A 44-gene penalised Cox score, developed under a hash-locked plan with every data-dependent step confined to training data, was associated with overall survival in a held-out internal test set (HR per SD 1.93) and, with frozen genes and coefficients, in two external cohorts on a different platform and tissue type (HR per SD 1.69 and 1.49). Leakage and negative controls behaved as expected (permuted-label C ≈ 0.50). By the pre-specified conclusion rules the score is "internally validated" and "externally validated". Discrimination was modest throughout (C 0.63–0.66), comparable to COO alone and below IPI + COO.

### 4.2 Added value beyond IPI and COO: weaker than the label

The pre-specified sentence "adds information beyond IPI and COO" is formally met, but the evidence is weaker than this wording suggests. The independent contribution was statistically supported only in the internal test set (adjusted HR 1.63, LRT p = 0.007). In GSE10846 and GSE87371 the adjusted HRs (1.29 and 1.17) were above 1, but their CIs included 1 and the LRTs were not significant (p = 0.134 and 0.434). The external condition in the plan required only a point estimate > 1, which is a weak criterion. ΔC was small in every cohort (+0.008 to +0.020) and its 95% CI included 0 each time; because these ΔC values are apparent (the same data fit and evaluate the extended model), true gains are probably smaller. Adjustment for IPI + COO shrank the HR in every cohort (test 2.18 → 1.63; GSE10846 1.88 → 1.29; GSE87371 1.49 → 1.17), and the LP correlated with IPI (ρ 0.28–0.40) and COO (ρ 0.46 with Lymph2Cx). The fair summary is: **independent contribution is directionally consistent across cohorts but was not demonstrated externally, and no clinically meaningful improvement in discrimination was shown.**

Why the test set alone showed added value is unclear. Chance is plausible (54 IPI-complete deaths). The test set also shares platform, laboratory and COO classifier with the training set, so the score may capture DASL-specific residual information about COO or case mix that does not transfer. In GSE87371, IPI alone discriminated better than the score (C 0.731 vs 0.630), and removing the 13 R-mini-CHOP21 patients lowered the score's HR from 1.49 to 1.39 (CI including 1); part of the score may overlap with age or frailty, which IPI partly captures. With 36–41 deaths these differences cannot be resolved.

### 4.3 Validation performance above the nested-CV estimate

In all three validation sets, discrimination exceeded the nested-CV estimate from training (C 0.612): test +0.052, GSE10846 +0.047, GSE87371 +0.018; the test-set HR (1.93) and the test-set adjusted HR (1.63) also exceeded what training suggested (≈ 1.4–1.5 and 1.22). None crossed the pre-specified red-flag threshold (+0.10), and the leakage and negative controls were clean. Possible explanations, none proven:

1. **Chance.** Each validation C has a wide CI (test 0.596–0.730; GSE10846 0.592–0.732; GSE87371 0.536–0.718) that includes 0.612. However, the excess has the same sign in three sets, which makes pure chance a less complete explanation.
2. **Pessimistic nested CV.** Nested-CV fold models were fitted on \~80% of the training data (\~298 patients) and were sparser (median 19 genes) than the final model (373 patients, 44 genes). Learning-curve effects would make nested CV understate the final model's performance.
3. **Case mix and follow-up.** Harrell's C depends on the follow-up and risk distribution of each cohort. GSE10846 has short follow-up (median 34.2 months), so its C is dominated by early deaths, which may be easier to rank. Within-cohort z-scoring also adapts the score to each cohort's case mix.
4. **A known bias in the opposite direction.** Three unsupervised steps (quantile reference, probe selection, low-expression filter) used all 380 training samples outside CV, which would make nested CV slightly optimistic, not pessimistic. It therefore does not explain the gap.

The validation estimates should be reported with their intervals, not as the "true" performance; the most defensible single expectation for new data remains C ≈ 0.6–0.65.

### 4.4 Biology: mostly known programmes

The score largely re-discovers established DLBCL biology. Its high-risk pole is enriched for proliferation (E2F, G2M), MYC targets and an ABC/IRF4 programme; its low-risk pole for extracellular matrix/stroma and T-cell genes. This matches the proliferation (adverse) and lymph-node/stromal (favourable) signatures of Rosenwald et al. \[3\], the favourable stromal-1 (ECM, histiocytic) signature after R-CHOP of Lenz et al. \[5\], the adverse ABC phenotype \[2,3,4\], IRF4 target genes described by Shaffer et al. \[25\], and the MYC-related high-grade biology of the MHG group \[6\]; MYC and BCL2 co-expression at protein level is also adverse after R-CHOP \[26\]. The score correlated in the expected direction with an outcome signature derived in an independent cohort by Shipp et al. \[24\] (ρ +0.37), but because that signature was itself derived from outcome, this is concordance, not independent biological evidence. Together, six signatures explained about 35% of LP variance; the remaining \~65% is unexplained and may contain noise as well as additional biology.

**Unresolved contradiction (stromal-2).** Lenz et al. reported an angiogenesis-related stromal-2 signature associated with worse survival \[5\]. Our angiogenesis proxy (Hallmark angiogenesis, not the original gene list) instead correlated with the low-risk pole (ρ −0.20). Because the proxy overlaps with ECM genes, stromal-1 and stromal-2 cannot be separated here. We report this as a contradiction, not as explained.

**T-cell signal (hypothesis).** T-cell genes (CD2, LCK, ICOS) were higher at the low-risk pole, and JAML — one of the four stable genes — had a negative coefficient. JAML is a costimulatory receptor on epithelial γδ T cells \[27\] and is expressed by activated/tumour-infiltrating CD8 and tissue-resident memory CD8 T cells \[28,29\]. These observations point in the same direction, but Hallmark IFN-γ correlated only weakly (ρ −0.11) and we found no direct evidence for JAML in DLBCL. The functions of PXDNL and PLAC8L1, two of the other stable genes, are essentially unknown, and no biological interpretation is offered.

### 4.5 Unstable gene list, stable signal

Prediction was robust to tuning (C ≈ 0.61 across pure LASSO and no variance filter), but the gene list was not (Jaccard 0.18–0.33; only 4/44 genes selected in ≥ 50% of bootstraps; other CV seeds chose 11–44 genes). This is the expected behaviour of sparse regression when many correlated genes carry small, overlapping effects. What was validated is the 44-gene score as a whole, not individual genes; individual-gene claims should not be made.

### 4.6 Comparison with common practice

The study followed reporting and design recommendations for prediction models \[7–9\]: pre-registration in a locked internal plan, training-only feature selection, nested CV, frozen-model external validation, label-permutation controls, and explicit comparison with IPI + COO. These choices produced modest but replicable associations rather than the high apparent performance often reported for signatures selected on full datasets.

## 5 Limitations

**Data and processing**

- Author-processed GEO values were used; CEL files were not re-processed, and pre-processing differs between sets (GSE87371 appears GCRMA-like, median 2.55).
- Development on FFPE/DASL and external validation on frozen/Affymetrix tissue; lower external performance was expected. The platform difference is also a strength for generalisability but prevents transporting absolute scores.
- Three unsupervised pre-processing steps used all 380 training samples outside CV, so nested CV may be slightly optimistic (plan §19); test results are unaffected.
- 1,518 of 13,747 genes had different symbols across platforms (same Entrez ID), probably due to annotation dates; current GPL570 symbols were used for interpretation.
- The GSE181063 sex check relied on four Y genes (no XIST), and the sex rule was set after the first QC run (without survival data). The meaning of the GEO field `sample_pre_active_treatment` is not fully certain, and a sensitivity analysis including non-pre-treatment biopsies was not done.

**Design and analysis**

- The internal test set comes from the same cohort, platform and laboratory as training; batch, platform and centre effects are tested only externally.
- External cohorts are small (58 and 41 deaths after QC); GSE10846 has short follow-up (median 34.2 months), so 5-year estimates are unstable; in GSE87371 the HR CI included 1 after removing R-mini-CHOP21 patients.
- Complete-case analysis for IPI excluded 36 test and 75 GSE10846 patients; in the test set missingness may not be random (death proportion 55.6% vs 42.9%, p = 0.074), with unknown direction of bias.
- IPI was modelled as linear 0–5; components were not modelled separately (unavailable in GSE87371).
- COO classifiers differ between cohorts, and MHG exists only in GSE181063 (24 training patients). COO labels come from the authors' classifier applied to the same arrays, so label-consistency checks are not independent validation.
- ΔC values and their bootstrap CIs are apparent and share the same optimism.
- Calibration was assessed only internally (slope 1.56, 0.95–2.16); within-cohort z-scoring makes absolute risk non-transportable, so external calibration and clinical-utility measures (e.g., decision curves) were not possible.
- H2 and H3 were secondary and unadjusted for multiplicity; the permutation p-value is bounded at 0.0099 by 100 permutations.
- No fairness or subgroup-performance analysis was done; race/ethnicity were unavailable. All cohorts were treated with R-CHOP-family regimens; no generalisation to other therapies is possible.
- No formal sample-size calculation was made for development.

**Model and interpretation**

- Gene selection is unstable; the 44-gene model depends on the planned CV seed, and other seeds select sparser nested models.
- The biological analysis is exploratory, training-only and in-sample; signature scores are simple mean z-scores, not the original weighted scores; stromal proxies are not the Lenz 2008 gene lists; |ρ| < 0.2 should be treated as weak.
- **The score cannot be applied to a single new patient.** Because each cohort is z-scored internally, the score is a relative ranking within a cohort, depends on the cohort's case mix (e.g., ABC proportion), and has no transportable threshold or baseline hazard. It is not a clinical tool.

## 6 Conclusion

Under a pre-registered, leakage-controlled pipeline, a 44-gene elastic-net Cox score for overall survival in DLBCL was **internally validated** (test HR per SD 1.93) and **externally validated** across platforms (GSE10846 HR 1.69; GSE87371 HR 1.49), meeting H1 and H2. By the plan's wording it also **adds information beyond IPI and COO**, but this was shown statistically only in the internal test set: externally the adjusted effect was directionally consistent but not significant, and discrimination gains were small with confidence intervals including zero. The score mostly re-captures known proliferation/MYC, ABC/IRF4, stromal and T-cell biology, and its individual genes are unstable. The main contribution of this work is methodological transparency rather than a new prognostic tool. This is a research study; it makes no diagnostic, prognostic or therapeutic claim for individual patients.

### Declarations

- **Ethics:** public, de-identified GEO data only; no new approval sought.
- **Funding:** none.
- **Conflicts of interest:** none declared.
- **Patient and public involvement:** none.
- **Protocol/registration:** internal hash-locked analysis plan v1.0 (`results/04/analysis_plan_04.md`, lock SHA-256 `8d87cf887e2db4770ae3273755d4901cd0d5a2a010fdc729159e007aa129bc7d`); not registered in an external registry.
- **Data availability:** raw data are available from GEO (GSE181063, GSE10846, GSE87371); raw data are not redistributed.
- **Code availability:** notebooks, the locked analysis plan, manifests, summary results and figures are public at <https://github.com/gurbuzmericarda0-arch/dlbcl-prognostic-score> (MIT licence for code; CC BY 4.0 for text, figures and summary results).
- **Use of AI:** analysis code and report drafting were assisted by an AI model (Claude, Anthropic); all numbers come from the project's run logs.

## References

1. International Non-Hodgkin's Lymphoma Prognostic Factors Project. A predictive model for aggressive non-Hodgkin's lymphoma. N Engl J Med. 1993;329(14):987–994.
2. Alizadeh AA, Eisen MB, Davis RE, et al. Distinct types of diffuse large B-cell lymphoma identified by gene expression profiling. Nature. 2000;403(6769):503–511.
3. Rosenwald A, Wright G, Chan WC, et al. The use of molecular profiling to predict survival after chemotherapy for diffuse large-B-cell lymphoma. N Engl J Med. 2002;346(25):1937–1947.
4. Scott DW, Wright GW, Williams PM, et al. Determining cell-of-origin subtypes of diffuse large B-cell lymphoma using gene expression in formalin-fixed paraffin-embedded tissue. Blood. 2014;123(8):1214–1217.
5. Lenz G, Wright G, Dave SS, et al. Stromal gene signatures in large-B-cell lymphomas. N Engl J Med. 2008;359(22):2313–2323. (Source study of GSE10846.)
6. Sha C, Barrans S, Cucco F, et al. Molecular high-grade B-cell lymphoma: defining a poor-risk group that requires different approaches to therapy. J Clin Oncol. 2019;37(3):202–212.
7. Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. BMJ. 2024;385:e078378. [doi:10.1136/bmj-2023-078378](https://doi.org/10.1136/bmj-2023-078378)
8. McShane LM, Altman DG, Sauerbrei W, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). J Natl Cancer Inst. 2005;97(16):1180–1184.
9. Altman DG, McShane LM, Sauerbrei W, Taube SE. Reporting Recommendations for Tumor Marker Prognostic Studies (REMARK): explanation and elaboration. PLoS Med. 2012;9(5):e1001216.
10. Lacy SE, Barrans SL, Beer PA, et al. Targeted sequencing in DLBCL, molecular subtypes, and outcomes: a Haematological Malignancy Research Network report. Blood. 2020;135(20):1759–1771. [doi:10.1182/blood.2019003535](https://doi.org/10.1182/blood.2019003535) PMID 32187361. GEO-listed citation for GSE181063, together with: Painter D, Barrans S, Lacy S, Smith A, et al. Cell-of-origin in diffuse large B-cell lymphoma: findings from the UK's population-based Haematological Malignancy Research Network. Br J Haematol. 2019;185(4):781–784. PMID 30408148.
11. Dubois S, Viailly P-J, Mareschal S, et al. Biological and clinical relevance of associated genomic alterations in MYD88 L265P and non-L265P-mutated diffuse large B-cell lymphoma: analysis of 361 cases. Clin Cancer Res. 2017;23(9):2232–2244. ([HAL record](https://hal.archives-ouvertes.fr/hal-01526582); PMID 27923841.) GEO-listed citation for GSE87371, together with: Dubois S, Tesson B, Mareschal S, Viailly PJ, et al. Refining diffuse large B-cell lymphoma subgroups using integrated analysis of molecular profiles. EBioMedicine. 2019;48:58–69. PMID 31648986. GSE87371 is a sub-series of GSE136971.
12. Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for functional genomics data sets—update. Nucleic Acids Res. 2013;41(D1):D991–D995.
13. Simon N, Friedman J, Hastie T, Tibshirani R. Regularization paths for Cox's proportional hazards model via coordinate descent. J Stat Softw. 2011;39(5):1–13.
14. Verweij PJM, van Houwelingen HC. Cross-validation in survival analysis. Stat Med. 1993;12(24):2305–2314.
15. Pölsterl S. scikit-survival: a library for time-to-event analysis built on top of scikit-learn. J Mach Learn Res. 2020;21(212):1–6.
16. Davidson-Pilon C. lifelines: survival analysis in Python. J Open Source Softw. 2019;4(40):1317.
17. Harrell FE Jr, Lee KL, Mark DB. Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. Stat Med. 1996;15(4):361–387.
18. Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. Stat Med. 2011;30(10):1105–1117.
19. Schoenfeld DA. Sample-size formula for the proportional-hazards regression model. Biometrics. 1983;39(2):499–503.
20. Grambsch PM, Therneau TM. Proportional hazards tests and diagnostics based on weighted residuals. Biometrika. 1994;81(3):515–526.
21. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database hallmark gene set collection. Cell Syst. 2015;1(6):417–425.
22. Kim SY, Volsky DJ. PAGE: parametric analysis of gene set enrichment. BMC Bioinformatics. 2005;6:144.
23. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57(1):289–300.
24. Shipp MA, Ross KN, Tamayo P, et al. Diffuse large B-cell lymphoma outcome prediction by gene-expression profiling and supervised machine learning. Nat Med. 2002;8(1):68–74.
25. Shaffer AL, Emre NCT, Lamy L, et al. IRF4 addiction in multiple myeloma. Nature. 2008;454(7201):226–231. (IRF4 target genes were defined in myeloma; their use in DLBCL is by analogy.)
26. Johnson NA, Slack GW, Savage KJ, et al. Concurrent expression of MYC and BCL2 in diffuse large B-cell lymphoma treated with rituximab plus cyclophosphamide, doxorubicin, vincristine, and prednisone. J Clin Oncol. 2012;30(28):3452–3459.
27. Witherden DA, Verdino P, Rieder SE, et al. The junctional adhesion molecule JAML is a costimulatory receptor for epithelial γδ T cell activation. Science. 2010;329(5996):1205–1210. ([PMC2943937](https://pmc.ncbi.nlm.nih.gov/articles/PMC2943937))
28. McGraw JM, et al. JAML promotes CD8 and γδ T cell antitumor immunity and is a novel target for cancer immunotherapy. J Exp Med. 2021;218(10):e20202644. ([journal page](https://rupress.org/jem/article/218/10/e20202644/212597/JAML-promotes-CD8-and-T-cell-antitumor-immunity))
29. Eschweiler S, Wang A, Ramírez-Suástegui C, et al. JAML immunotherapy targets recently activated tumor-infiltrating CD8+ T cells. Cell Rep. 2023;42(2):112040. [doi:10.1016/j.celrep.2023.112040](https://doi.org/10.1016/j.celrep.2023.112040)
30. Tibshirani R. The lasso method for variable selection in the Cox model. Stat Med. 1997;16(4):385–395.
31. Breslow N. Covariance analysis of censored survival data. Biometrics. 1974;30(1):89–99.
32. Efron B. The efficiency of Cox's likelihood function for censored data. J Am Stat Assoc. 1977;72(359):557–565.
33. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. J Mach Learn Res. 2011;12:2825–2830.

*Verification note:* references 7, 10, 11, 27, 28 and 29 were checked online on 2026-10-16; the GEO-listed citations of GSE181063 and GSE87371 (refs 10, 11) were confirmed on the GEO records on 2026-10-17. The others are standard citations entered from knowledge; volume/page details should be checked against PubMed before submission.

## List of figures and tables

All paths are relative to the repository root. Step-2 aggregate outputs are in `results/02/`; raw and processed expression data and patient-level clinical data are not redistributed.

### Main figures

| No. | Content | File |
| --- | --- | --- |
| Fig. 1 | Study design and patient flow (cohorts → eligibility → split → QC → analysis populations) | `figures/09b/flow_diagram_09b.png` |
| Fig. 2 | Positive control: KM, ABC vs GCB (training) | `figures/03/km_ABC_vs_GCB.png` |
| Fig. 3 | (a) CV deviance curve with λ\_min/λ\_1se; (b) model coefficients | `figures/04/cv_deviance_04b.png`; `figures/04/coefficients_04b.png` |
| Fig. 4 | Internal validation: KM by training-median cut-off (test) | `figures/05/km_test_median_05.png` |
| Fig. 5 | Negative control: test C of 100 permutation models vs observed | `figures/05/negative_control_05.png` |
| Fig. 6 | External validation: KM by cohort median | `figures/06/km_GSE10846_median_06.png`; `figures/06/km_GSE87371_median_06.png` |
| Fig. 7 | Forest plot: HR per SD in test, external cohorts and sensitivity analyses | `figures/06/forest_hr_per_sd_06.png` |
| Fig. 8 | Forest plot: LP HR per SD unadjusted and adjusted for IPI, COO, IPI + COO | `figures/07/forest_adjusted_hr_07.png` |
| Fig. 9 | (a) Hallmark score-based gene-set test; (b) signature–LP correlations and LP vs Lymph2Cx | `figures/08/hallmark_rank_test_08.png`; `figures/08/signatures_vs_score_08.png` |

### Supplementary figures

| No. | Content | File |
| --- | --- | --- |
| S1 | Sample QC: normalisation boxplots, PCA (training), sex check | `figures/02b/boxplots_normalisation.png`; `figures/02b/pca_train.png`; `figures/02b/sex_check.png` |
| S2 | KM by all four COO groups (training) | `figures/03/km_all_COO.png` |
| S3 | Lymph2Cx score by COO label | `figures/03/coo_score_by_label.png` |
| S4 | Nested-CV C distribution (λ\_min, λ\_1se) | `figures/04/nested_cv_c_04b.png` |
| S5 | Cross-validated KM (repeat 1) | `figures/04/cv_km_repeat1_04b.png` |
| S6 | Gene-selection stability (30 most often selected) | `figures/04/stability_top30_04b.png` |
| S7 | Leakage check (permuted nested-CV C) | `figures/05/leakage_check_05.png` |
| S8 | KM by tertile and LP distribution (test vs training) | `figures/05/km_test_tertiles_05.png`; `figures/05/lp_distribution_05.png` |
| S9 | KM by tertile and LP distributions (external) | `figures/06/km_GSE10846_tertiles_06.png`; `figures/06/km_GSE87371_tertiles_06.png`; `figures/06/lp_distribution_06.png` |
| S10 | Harrell's C of the clinical/LP models | `figures/07/c_index_models_07.png` |
| S11 | ΔC bootstrap distributions | `figures/07/delta_c_bootstrap_07.png` |
| S12 | 44 model genes with coefficients and bootstrap frequencies | `figures/08/model_genes_coefficients_08.png` |

### Tables

| Report table | Content | Source file(s) |
| --- | --- | --- |
| Table 1 | Cohorts after eligibility | `results/02/cohort_summary.csv`, `results/02/patient_flow.csv` |
| Table 2 | Sample QC flags | sample-level QC flags (not redistributed; regenerate with notebook 02b) |
| Table 3 | OS by COO (positive control) | `results/03/km_summary_all_COO.csv`, `km_summary_ABC_GCB.csv`, `logrank_results.csv`, `cox_results.csv` |
| Table 4 | Training sensitivity analyses | `results/04/sensitivity_04b.csv` (also `nested_cv_summary_04b.json`, `stability_04b.csv`, `added_value_04b.json`) |
| Table 5 | Internal validation | `results/05/internal_validation_05*`, `km_test_05*`, `negative_control_05*`, `sensitivity_05*`, `first_results_05.json` |
| Table 6 | External validation | `results/06/external_validation_06.json`, `km_external_06.csv`, `sensitivity_06.csv`, `first_results_06.json` |
| Table 7 | C-indices and ΔC | `results/07/c_index_models_07.csv`, `delta_c_bootstrap_07.csv` |
| Table 8 | LP HRs and LRTs | `results/07/hr_lp_07.csv`, `lrt_07.csv`, `missingness_07.csv` |
| Table 9 | Signature correlations | `results/08/signature_associations_08.csv`, `rank_test_08.csv`, `ora_08.csv` |
| Suppl. | Final model (44 genes, coefficients, training mean/SD) | `results/04/model_04.json` (SHA-256 `661e2b80be77…`), `model_genes_04b.csv`, `results/08/model_genes_annotated_08.csv` |
| Suppl. | Analysis plan and lock | `results/04/analysis_plan_04.md`, `analysis_plan_04.json`, `plan_lock_04.json`, `plan_amendments.md` |
| Suppl. | Run provenance | `manifest_02a.json`, `manifest_02b.json`, `results/03/manifest_03.json`, `results/04/manifest_04b.json`, `results/05/manifest_05.json`, `results/06/manifest_06.json`, `results/07/manifest_07.json`, `results/08/manifest_08.json` |

Patient-level score files (`results/05/test_lp_05.csv`, `results/06/external_lp_06.csv`) contain GEO sample accessions and derived scores only; no clinical data are redistributed.

## Appendix A. TRIPOD+AI checklist

Item topics are paraphrased from the TRIPOD+AI checklist (version 11 January 2024) \[7\]; check numbering and wording against the official checklist before submission. Status: Met / Partly / Not met / N/A.

| Item | Topic (paraphrased) | Status | Where in this report / gap |
| --- | --- | --- | --- |
| 1 | Title identifies development/validation, population, outcome | Met | Working title |
| 2 | Abstract (TRIPOD+AI for Abstracts) | Partly | Abstract; not checked item by item against the 13-item abstract checklist |
| 3a | Healthcare context, rationale, existing models | Partly | §1; no specific published DLBCL signature compared head-to-head |
| 3b | Target population, intended purpose, intended users | Met | §1 Intended use (none in patient care) |
| 3c | Known health inequities | Not met | Not addressed; no sociodemographic data beyond age/sex |
| 4 | Objectives (development and/or validation) | Met | §1 Objectives |
| 5a | Data sources | Met | §2.1 |
| 5b | Dates of data collection and follow-up | Not met | Accrual periods not extracted from GEO/source papers |
| 6a | Setting and centres | Partly | §2.1 (population-based HMRN, multicentre trials); centre-level detail missing |
| 6b | Eligibility criteria | Met | §2.2 |
| 6c | Treatments received | Met | §2.2 |
| 7 | Data preparation and quality checks | Met | §2.2, §2.4 |
| 8a | Outcome definition and time horizon | Met | §2.3 |
| 8b | Rationale/qualifications for outcome assessment | N/A | OS is objective |
| 8c | Blinding of outcome assessment | Met | §2.3 (not blinded) |
| 9a | Choice of candidate predictors | Met | §2.7 |
| 9b | Definition and measurement of predictors | Met | §2.4, §2.7, §2.10 |
| 9c | Blinding of predictor assessment | Partly | Expression measured without outcome; COO from authors' classifiers (§5) |
| 10 | Sample size | Partly | §2.12: validation power given; no development calculation |
| 11 | Missing data | Met | §2.4 (none in expression), §2.10 (complete case + diagnostics), §3.6 |
| 12a | How data were used (development vs evaluation) | Met | §2.5, §2.7–2.9 |
| 12b | Predictor handling | Met | §2.7 (variance filter, z-score inside folds), §2.9 |
| 12c | Model type, building, tuning, internal validation | Met | §2.7 |
| 12d | Heterogeneity across clusters/centres | Partly | Results per cohort; no formal heterogeneity analysis |
| 12e | Performance measures | Partly | §2.8–2.10; calibration internal only |
| 12f | Model updating/recalibration | N/A | No updating; within-cohort z-scoring pre-specified (§2.9) |
| 12g | How predictions were calculated in evaluation | Met | §2.8, §2.9 |
| 13 | Class imbalance methods | N/A | None used (survival outcome) |
| 14 | Fairness approaches | Not met | No fairness analysis (§2.13, §5) |
| 15 | Model output and thresholds | Met | §2.7–2.9 (LP; median cut-offs, no optimisation) |
| 16 | Differences between development and evaluation data | Met | §2.1, §2.9, §5 |
| 17 | Ethical approval | Met | §2.13 |
| 18a | Funding | Met | Declarations |
| 18b | Conflicts of interest | Met | Declarations |
| 18c | Protocol | Met | Locked plan v1.0, lock file and amendment log public in `results/04/` (lock not externally registered; see 18d) |
| 18d | Registration | Not met | No external registry |
| 18e | Data sharing | Met | Declarations (GEO accessions) |
| 18f | Code sharing | Met | <https://github.com/gurbuzmericarda0-arch/dlbcl-prognostic-score> |
| 19 | Patient and public involvement | Met | Declarations (none) |
| 20a | Participant flow | Met | §3.1; Fig. 1 (`figures/09b/flow_diagram_09b.png`) |
| 20b | Participant characteristics | Met | Table 1; Supplementary Table S1 (age, sex, IPI group, COO by analysis population; stage only within IPI) |
| 20c | Comparison of development and evaluation data | Not met | No side-by-side predictor/outcome distribution table |
| 21 | Participants and events in each analysis | Met | §3.1–3.6 |
| 22 | Full model specification | Met | Supplementary Table S2 (44 genes, coefficients, training means/SDs); `results/04/model_04.json` |
| 23a | Performance with CIs | Met | Tables 5–8 (external calibration not possible) |
| 23b | Heterogeneity in performance | Partly | Per cohort; no pooled/meta estimate |
| 24 | Model updating results | N/A | No updating |
| 25 | Interpretation (incl. fairness) | Partly | §4; fairness not discussed |
| 26 | Limitations | Met | §5 |
| 27a | Handling poor-quality or missing input data in use | Partly | QC rules §2.2; not intended for use |
| 27b | User expertise/interaction | N/A | Not intended for use |
| 27c | Implications and future research | Partly | §4, §6; future-research directions brief |

**Not met:** 3c, 5b, 14, 18d, 20c. Items 18c, 18f, 20a, 20b and 22 were completed in step 9b.

## Appendix B. REMARK checklist

Item topics paraphrased from McShane et al. \[8\] and Altman et al. \[9\].

| Item | Topic (paraphrased) | Status | Where in this report / gap |
| --- | --- | --- | --- |
| 1 | Marker, objectives, pre-specified hypotheses | Met | §1, §2.8–2.10 (H1–H3) |
| 2 | Patient characteristics, source, inclusion/exclusion | Met | §2.2, Table 1, Supplementary Table S1 |
| 3 | Treatments received | Met | §2.2 |
| 4 | Specimen type, preparation, storage | Partly | FFPE vs frozen (§2.1); storage and preparation details not reported |
| 5 | Assay method and quality control | Partly | Platforms, processing, QC (§2.2, §2.4); author-processed data, no re-processing |
| 6 | Study design: case selection, period, follow-up | Partly | §2.2, §3.1; accrual dates missing |
| 7 | Clinical endpoints | Met | §2.3 |
| 8 | All candidate variables | Met | §2.7, §2.10 |
| 9 | Sample-size rationale | Partly | §2.12 |
| 10 | Statistical methods, cut-points, missing data | Met | §2.7–2.10 |
| 11 | Handling of marker values | Met | §2.7, §2.9 |
| 12 | Patient flow | Met | §3.1, Table 2, Fig. 1 |
| 13 | Distribution of demographics, standard prognostic variables, marker | Met | Supplementary Table S1; LP distributions in figures |
| 14 | Relation of marker to standard prognostic variables | Met | §3.6, §3.7 (ρ with IPI, COO) |
| 15 | Univariable analysis with HR and CI | Met | Tables 5, 6, 8 |
| 16 | Multivariable analysis | Met | Table 8 |
| 17 | Further investigations (e.g., interactions, subgroups) | Partly | MHG and treatment sensitivity analyses; no interaction tests |
| 18 | Validation | Met | §3.4, §3.5 |
| 19 | Interpretation in context of hypotheses and literature | Met | §4 |
| 20 | Implications for future research and clinical value | Met | §4, §6 (no clinical value claimed) |

**Not fully met:** items 4, 5, 6, 9 and 17 are partial; none is entirely missing.

## Appendix C. Supplementary tables (added in step 9b)

**Supplementary Table S1. Characteristics of the analysis populations (after QC exclusions).** Percentages use the population size as denominator. Follow-up by reverse Kaplan–Meier.

| Characteristic | GSE181063 training | GSE181063 test | GSE10846 (external) | GSE87371 (supporting external) |
| --- | --- | --- | --- | --- |
| Patients | 373 | 162 | 229 | 116 |
| Deaths | 171 (45.8%) | 74 (45.7%) | 58 (25.3%) | 41 (35.3%) |
| Median follow-up, months (reverse KM) | 96.9 | 101.3 | 34.2 | 49.2 |
| Age, years: median (IQR) | 67 (58–75) | 67 (58–75) | 62 (51–74) | 67 (63–76) |
| Age, years: > 60 | 264 (70.8%) | 114 (70.4%) | 119 (52.0%) | 97 (83.6%) |
| Sex: F | 188 (50.4%) | 84 (51.9%) | 98 (42.8%) | 58 (50.0%) |
| Sex: M | 185 (49.6%) | 78 (48.1%) | 131 (57.2%) | 58 (50.0%) |
| IPI: 0–1 | 85 (22.8%) | 32 (19.8%) | 66 (28.8%) | 12 (10.3%) |
| IPI: 2 | 67 (18.0%) | 41 (25.3%) | 39 (17.0%) | 27 (23.3%) |
| IPI: 3 | 80 (21.4%) | 31 (19.1%) | 28 (12.2%) | 35 (30.2%) |
| IPI: 4–5 | 49 (13.1%) | 22 (13.6%) | 21 (9.2%) | 42 (36.2%) |
| IPI: missing | 92 (24.7%) | 36 (22.2%) | 75 (32.8%) | 0 (0.0%) |
| COO: ABC | 104 (27.9%) | 50 (30.9%) | 91 (39.7%) | 56 (48.3%) |
| COO: GCB | 185 (49.6%) | 73 (45.1%) | 105 (45.9%) | 42 (36.2%) |
| COO: MHG | 23 (6.2%) | 12 (7.4%) | — | — |
| COO: UNC | 61 (16.4%) | 27 (16.7%) | 33 (14.4%) | 18 (15.5%) |

**Supplementary Table S2. Final model: 44 genes ordered by absolute coefficient.** LP = Σ coef × (x − training mean) / training SD, with x the log2 expression after step-2b processing; positive coefficient = higher expression, higher risk. Bootstrap selection = proportion of 100 bootstrap refits selecting the gene. Training-median cut-off -0.0048; λ = 0.086259. Machine-readable: `results/09b/table_s2_model_coefficients_09b.csv`.

| # | Gene (current symbol) | Entrez | Coefficient (log HR per training SD) | HR per SD | Bootstrap selection | Training mean (log2) | Training SD |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | JAML | 120425 | -0.1637 | 0.849 | 68% | 12.199 | 0.953 |
| 2 | PXDNL | 137902 | -0.1367 | 0.872 | 84% | 9.258 | 1.462 |
| 3 | PLAC8L1 | 153770 | +0.0803 | 1.084 | 60% | 9.521 | 0.995 |
| 4 | FRRS1 | 391059 | +0.0706 | 1.073 | 43% | 11.391 | 0.882 |
| 5 | ZBTB7C | 201501 | -0.0694 | 0.933 | 28% | 8.762 | 0.935 |
| 6 | PALD1 | 27143 | -0.0690 | 0.933 | 28% | 11.466 | 0.841 |
| 7 | CXCL2 | 2920 | +0.0683 | 1.071 | 44% | 11.421 | 1.125 |
| 8 | GPR141 | 353345 | -0.0655 | 0.937 | 39% | 10.599 | 1.046 |
| 9 | ATAD3B | 83858 | +0.0546 | 1.056 | 49% | 12.106 | 0.925 |
| 10 | HBD | 3045 | +0.0545 | 1.056 | 39% | 8.745 | 0.878 |
| 11 | HLX | 3142 | +0.0544 | 1.056 | 28% | 10.912 | 0.835 |
| 12 | GPR107 | 57720 | +0.0457 | 1.047 | 43% | 10.217 | 0.972 |
| 13 | JADE3 | 9767 | +0.0456 | 1.047 | 41% | 11.461 | 1.241 |
| 14 | TIMM50 | 92609 | +0.0451 | 1.046 | 32% | 8.652 | 0.827 |
| 15 | CD300LG | 146894 | +0.0365 | 1.037 | 20% | 10.154 | 1.339 |
| 16 | AGBL4 | 84871 | +0.0345 | 1.035 | 50% | 8.996 | 1.087 |
| 17 | SPINT2 | 10653 | -0.0344 | 0.966 | 29% | 11.029 | 1.059 |
| 18 | LDLRAD4 | 753 | -0.0339 | 0.967 | 38% | 9.341 | 1.042 |
| 19 | TMPRSS5 | 80975 | -0.0334 | 0.967 | 34% | 8.446 | 0.852 |
| 20 | SHROOM3 | 57619 | -0.0310 | 0.969 | 35% | 10.152 | 1.519 |
| 21 | HMGN2P46 | 283651 | -0.0280 | 0.972 | 36% | 9.571 | 1.175 |
| 22 | ADGRG1 | 9289 | -0.0278 | 0.973 | 22% | 10.078 | 1.124 |
| 23 | L1TD1 | 54596 | +0.0276 | 1.028 | 41% | 8.980 | 1.082 |
| 24 | GOSR1 | 9527 | -0.0254 | 0.975 | 26% | 9.895 | 0.844 |
| 25 | ENPP5 | 59084 | +0.0249 | 1.025 | 26% | 11.438 | 1.091 |
| 26 | COL13A1 | 1305 | -0.0232 | 0.977 | 17% | 9.799 | 1.021 |
| 27 | FABP4 | 2167 | +0.0222 | 1.022 | 39% | 8.992 | 0.925 |
| 28 | CLEC4C | 170482 | -0.0198 | 0.980 | 32% | 8.723 | 1.040 |
| 29 | BMP8B | 656 | +0.0163 | 1.016 | 23% | 10.231 | 0.859 |
| 30 | BEND4 | 389206 | +0.0151 | 1.015 | 30% | 11.607 | 1.408 |
| 31 | CTTN | 2017 | -0.0110 | 0.989 | 34% | 10.267 | 0.904 |
| 32 | EVA1C | 59271 | -0.0105 | 0.990 | 18% | 10.227 | 0.914 |
| 33 | DSCAML1 | 57453 | +0.0088 | 1.009 | 28% | 9.104 | 0.928 |
| 34 | SCN2A | 6326 | +0.0087 | 1.009 | 33% | 9.476 | 1.737 |
| 35 | CCDC50 | 152137 | +0.0087 | 1.009 | 15% | 11.046 | 0.815 |
| 36 | WNT1 | 7471 | -0.0083 | 0.992 | 17% | 8.601 | 0.883 |
| 37 | MESP1 | 55897 | +0.0082 | 1.008 | 18% | 8.873 | 0.915 |
| 38 | CUX2 | 23316 | -0.0079 | 0.992 | 33% | 9.349 | 1.626 |
| 39 | AHSP | 51327 | +0.0078 | 1.008 | 37% | 9.988 | 1.152 |
| 40 | STK33 | 65975 | +0.0049 | 1.005 | 23% | 12.124 | 1.123 |
| 41 | CIRBP-AS1 | 148046 | +0.0034 | 1.003 | 18% | 11.060 | 0.952 |
| 42 | TMEM63C | 57156 | +0.0030 | 1.003 | 17% | 8.918 | 0.845 |
| 43 | GAB1 | 2549 | +0.0022 | 1.002 | 27% | 11.019 | 0.988 |
| 44 | NIPSNAP3B | 55335 | +0.0007 | 1.001 | 32% | 9.156 | 0.855 |
