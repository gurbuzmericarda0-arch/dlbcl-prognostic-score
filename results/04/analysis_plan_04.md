# Analysis plan — steps 4–7 (version 1.0)

**Project:** Lymphoma AI Research (`dlbcl-prognosis`) · **Drafted:** 2026-10-15 · **Master seed:** 20261015

**Status:** locked when `results/04/plan_lock_04.json` exists. The lock file stores the SHA-256 of this document and of `analysis_plan_04.json` (the parameters the analysis code reads), a timestamp, and the hashes of all input files.

This is a research project. No part of it is a diagnostic or treatment tool, and no result will be presented as one.

## 0. Data seen before this plan

| Data | Seen so far | Outcome data used so far |
|---|---|---|
| GSE181063 training (n = 380) | clinical cleaning, expression QC | descriptive counts (Step 2a); OS by cell of origin, with and without IPI adjustment (Step 3 positive control) |
| GSE181063 test (n = 164) | clinical cleaning, expression QC | descriptive counts only (Step 2a) |
| GSE10846, GSE87371 | clinical cleaning, expression QC | descriptive counts only (Step 2a) |

Descriptive counts: number of deaths and median follow-up per cohort. No association between any gene and survival has been computed in any cohort, and no outcome comparison has been made outside the training set.

## 1. Objectives

- **Primary.** Develop a penalized Cox model of overall survival (OS) from gene expression in the GSE181063 training set, freeze it, and test it once in the held-out GSE181063 test set (Step 5).
- **Secondary.** Test the frozen model in two external cohorts measured on a different platform (Step 6), and test whether it adds prognostic information to the IPI and cell of origin (COO) (Step 7).
- **Exploratory.** Biological interpretation of the selected genes (Step 8). Its methods will be set in Step 8; it is not confirmatory.

## 2. Cohorts and analysis populations

| Role | Cohort | Platform | Patients | Deaths | Expected QC exclusions | First used |
|---|---|---|---|---|---|---|
| Training | GSE181063, train split | GPL14951 (DASL, FFPE) | 380 | 173 | 7 | Step 4 |
| Internal test | GSE181063, test split | GPL14951 (DASL, FFPE) | 164 | 74 | 2 | Step 5 |
| External, main | GSE10846 | GPL570 (HG-U133 Plus 2.0) | 232 | 59 | 3 | Step 6 |
| External, supportive | GSE87371 | GPL570 (HG-U133 Plus 2.0) | 116 | 41 | 0 | Step 6 |

Patient and death counts are before the QC rule. Inclusion criteria are those fixed in Step 2a.

- **QC rule (primary, all cohorts):** exclude samples flagged `sex_mismatch_confident` in Step 2b. Their clinical data may belong to another sample.
- **QC sensitivity analysis:** additionally exclude samples flagged `corr_outlier`.
- All COO groups, including MHG and UNC, are used for model development. COO is not a model input.

## 3. Endpoints

- **Primary:** OS in months; event = death from any cause; censored at last follow-up; time origin as defined by the data providers. No administrative truncation during model development.
- **Secondary:** progression-free survival (PFS), only in cohorts whose GEO records provide it. PFS is not in the Step 2a clean files; if used, it is extracted from the raw GEO tables with the same inclusion rules, and its coding is checked (as in Step 2a) before any PFS analysis. The frozen OS model is applied unchanged. No model is trained on PFS.

## 4. Features and pre-processing

1. Input: the 13,747 genes common to GPL14951 and GPL570 (Step 2b), log2 and quantile-normalized, identified by **Entrez Gene ID**. Gene symbols are for display only.
2. Variance filter: keep the 5,000 genes with the highest variance in the training data of the current fit.
3. Standardization: gene-wise z-score with the mean and SD of the training data of the current fit.
4. "Current fit" means each cross-validation training fold, or the full training set for the final model. Steps 2–3 are re-estimated inside every fold, so held-out samples never influence them.
5. No supervised pre-screening (for example, univariable Cox filtering). Survival data are used only inside the penalized Cox fit.
6. Clinical variables (IPI, COO, age, sex) are not model inputs. They are used only for comparison (Step 7).

## 5. Model

- Cox proportional-hazards model with an elastic-net penalty, `l1_ratio` = 0.9: mostly LASSO; the small ridge part makes selection among correlated genes more stable. Fitted with scikit-survival `CoxnetSurvivalAnalysis`: 100 penalty values λ (called `alpha` in scikit-survival), `alpha_min_ratio` = 0.01, Breslow ties.
- **Penalty selection:** 10-fold cross-validation, folds stratified by OS event, seed 20261015. Criterion: cross-validated partial-likelihood deviance (Verweij & van Houwelingen 1993; the default in glmnet). **λ_min** (lowest mean CV deviance) is primary. λ_1se (largest penalty within one standard error of the minimum) is recorded for sensitivity only.
- **Final model:** fitted on the full QC-filtered training set at λ_min, then frozen and hashed before Step 5: Entrez IDs, coefficients, training means and SDs of the selected genes, penalty value and risk-group cutoff.
- **Risk score:** linear predictor, LP = Σ β_g · z_g over the selected genes.
- **If no gene is selected at λ_min,** this is reported as a null result. Settings are not changed to force genes into the model.
- If λ_min is the smallest penalty on the grid, this is reported and the grid is not extended.
- If the fit fails numerically at some penalty values in any fold, those values are dropped for all folds and this is reported.

## 6. Risk groups

- **Primary:** high vs low risk, split at the median LP of the training set. The same numeric cutoff is used in the GSE181063 test set.
- **External cohorts:** each cohort's own median LP (§10).
- **Descriptive only:** tertiles.
- Data-driven "optimal" cut-points (for example, maximally selected log-rank statistics) are not used.

## 7. Step 4 outputs (training set only)

1. The frozen final model, the CV deviance curve and the coefficients.
2. **Nested cross-validation** (an honest estimate inside the training set): outer 5-fold × 5 repeats, stratified by event. The whole procedure of §4–5 runs inside each outer training fold. Reported: mean and SD of Harrell's C over the 25 outer test folds, for λ_min and λ_1se, and a cross-validated Kaplan–Meier curve (Simon et al. 2011; repeat 1), in which each outer test fold is split at the median LP of its own outer training fold and the folds are then pooled. This estimate is descriptive. It does not decide whether the model goes to the test set: the primary model goes regardless.
3. **Added value inside the training set (descriptive):** the out-of-fold LP (repeat 1, standardized with the LP mean and SD of its outer training fold) added to Cox(IPI + COO) in patients with known IPI; likelihood-ratio test and ΔC.
4. **Gene stability:** 100 bootstrap resamples; §4 pre-processing and the Cox fit at the fixed λ_min of the final model; selection frequency of each gene. Used for interpretation in Step 8, not for selection.
5. **Training-only sensitivity analyses** (never taken to the test set): (a) λ_1se; (b) pure LASSO (`l1_ratio` = 1.0); (c) no variance filter (all 13,747 genes); (d) correlation outliers excluded; (e) five other CV seeds. Reported for each: number of genes and overlap with the primary model (Jaccard index); for (b) and (c) also nested-CV C (5-fold × 1 repeat, same folds as repeat 1 of the primary analysis).

## 8. Internal validation (Step 5)

Step 5 starts only after the SHA-256 of the frozen model is recorded. The model is not changed afterwards.

1. Test-set expression is z-scored with the **training** means and SDs (same platform, same normalization reference). Sensitivity: z-scoring within the test set.
2. Metrics:
   - **Primary:** Cox hazard ratio (HR) per 1 SD of LP (SD of the cohort's own LP), with 95% CI.
   - Harrell's C with 95% bootstrap CI (1,000 resamples, percentile method).
   - Uno's C, truncated at τ = 60 months. If fewer than 20% of a cohort's patients are still at risk (alive and under follow-up) at 60 months, τ is the latest time at which at least 20% are. The censoring distribution is estimated in the same cohort.
   - Kaplan–Meier curves for high vs low risk (training cutoff), log-rank test, HR high vs low.
   - Calibration slope: the Cox coefficient of the raw LP (1 = ideal; below 1 = overfitting).
   - Proportional-hazards check (Schoenfeld residuals), descriptive.
3. **H1 (the primary hypothesis):** HR per SD > 1 with a 95% CI that excludes 1 (two-sided α = 0.05). Supporting evidence: C-index CI excludes 0.5; log-rank p < 0.05; permutation p < 0.05 (§9).
4. QC sensitivity: correlation outliers excluded, without re-scaling.

## 9. Negative control (Step 5)

1. **Label permutation:** survival time and event are permuted together across training patients. This breaks any link with expression but keeps the survival distribution. The full procedure (§4–5, including CV for the penalty) is re-run on each of **100** permutations, and each permuted model is scored on the **true** test labels. Reported: mean and 2.5–97.5th percentiles of the permuted test C, the share of empty models, and the empirical p-value = (1 + number of permuted C ≥ observed C) / (1 + 100). An empty model has C = 0.5 by definition (all scores tied). Expected: permuted C centred at 0.5; observed model p < 0.05.
2. **Leakage check:** nested CV (5-fold × 1 repeat) on 10 permuted training sets. Expected mean C ≈ 0.5. If the mean is above 0.55, the analysis stops and the pipeline is searched for leakage.

## 10. External validation (Step 6)

1. The **frozen genes and frozen coefficients** are applied unchanged. No re-selection, no re-fitting, no re-tuning.
2. Expression is z-scored **within each cohort**, gene by gene, on its primary analysis population, because absolute expression levels are not comparable between platforms. What is not transferred is the absolute scale: the baseline hazard and the numeric LP cutoff.
3. Risk groups: each cohort's own median LP (the same quantile as in training).
4. Metrics as §8.2, except the calibration slope, which cannot be interpreted after within-cohort scaling.
5. **H2 (per cohort):** HR per SD > 1 with a 95% CI that excludes 1. GSE10846 is the main external test; GSE87371 is supportive (§14).
6. Sensitivity analyses, both using the scores of the primary analysis (no re-scaling): GSE87371 without R-mini-CHOP21 patients (n = 13); correlation outliers excluded in each cohort.

## 11. Comparison with IPI and cell of origin (Step 7)

1. Benchmarks, fitted in each validation cohort: IPI (score 0–5, linear, 1 df); COO (GCB as reference; ABC; other = UNC, MHG or unclassified); IPI + COO.
2. **H3 (added value):** Cox(IPI + COO) vs Cox(IPI + COO + LP per SD), likelihood-ratio test with 1 df. Primary in the test set: p < 0.05. Also reported: adjusted HR per SD with 95% CI; ΔHarrell's C with 95% bootstrap CI (1,000 resamples, both models refitted in each resample). Because the combination weights are fitted in the same cohort, ΔC slightly favours the larger model; the likelihood-ratio test accounts for the extra parameter.
3. **Missing IPI:** complete-case analysis. Diagnostics: patients with and without IPI are compared on age, sex, COO, follow-up and deaths, and the univariable HR of LP is reported in both groups.
4. COO sensitivity in GSE181063: MHG as a separate category.
5. Descriptive: C-index of IPI alone, COO alone, LP alone and their combinations, all in the same complete-case subset.
6. Reference point from Step 3 (training set, ABC/GCB patients with known IPI, n = 223): COO alone C = 0.614; IPI + COO C = 0.736. A gene score is not expected to beat IPI + COO on its own; the question is whether it adds to it (H3).

## 12. Conclusion rules (wording fixed in advance)

- "**Internally validated**": H1 met.
- "**Externally validated**": H2 met in GSE10846, and HR > 1 in GSE87371 (its CI may include 1).
- "**Adds information beyond IPI and COO**": H3 met in the test set, and adjusted HR > 1 in both external cohorts.
- Results that do not meet a rule are reported as negative or inconclusive, with the same detail as positive ones.

## 13. Red flags: stop and look for errors before reporting

- Any C-index > 0.90.
- Mean C of the permuted-label models > 0.55 (§9.2).
- Test-set C higher than the nested-CV mean C by more than 0.10.

## 14. Power (approximate)

Schoenfeld approximation for a continuous score per SD, two-sided α = 0.05: power = Φ(√d · ln HR − 1.96), where d = number of deaths.

| Cohort | Deaths | Smallest HR per SD detectable with 80% power | Power at HR 1.3 | Power at HR 1.5 |
|---|---|---|---|---|
| GSE181063 test | 74 | 1.38 | 0.62 | 0.94 |
| GSE10846 | 59 | 1.44 | 0.52 | 0.88 |
| GSE87371 | 41 | 1.55 | 0.39 | 0.74 |

GSE87371 has little power for moderate effects; a non-significant result there is not evidence against the model.

## 15. Multiplicity

H1 is the only primary hypothesis. All other tests are secondary, reported without adjustment and interpreted as supportive.

## 16. Data-access discipline

- Step 4 uses only: the training rows of the clinical table and of the QC table (other rows are dropped right after loading), the split file, the training expression file and the gene table.
- Step 4 never opens the test expression file or the external cohorts. They are only hashed, to prove later that they were not changed.
- The test set is opened in Step 5, after the frozen model's hash is recorded. The external cohorts are opened in Step 6.

## 17. Deviations

- The locked files are never edited. Any change is appended to `results/04/plan_amendments.md` with the date, the reason and which data had been seen at that time.
- After the test set is opened, changes are allowed only to correct errors, and both the original and the corrected results are then reported.

## 18. Reproducibility and reporting

- All random seeds derive from the master seed 20261015; the scheme is in `analysis_plan_04.json`. Software versions and SHA-256 hashes of inputs and outputs are recorded in each step's manifest.
- The report will follow the TRIPOD+AI and REMARK reporting guidelines.

## 19. Known limitations at the planning stage

- Three unsupervised steps of Step 2b used all 380 training samples, outside cross-validation: the quantile-normalization reference, the choice of one probe per gene (highest training mean) and the low-expression filter. None of them uses outcome data. They can make the nested-CV estimate (§7) very slightly optimistic; they do not bias the test-set results, because their parameters came from training samples only.
- The training data are FFPE tissue on a DASL array; the external cohorts are Affymetrix arrays. Lower external performance is expected.
- Within-cohort z-scoring depends on each cohort's case mix (for example, its share of ABC cases).
- COO classifiers differ between cohorts; MHG exists only in GSE181063.
- IPI is missing for part of GSE181063 and GSE10846; complete-case analyses can be biased if missingness is not random.
- The sex-check rule (Step 2b) was set after seeing the training QC results, without outcome data.
- The sensitivity analysis that includes non-pre-treatment biopsies (flag set in Step 2a; the exact meaning of that GEO field is uncertain) is not part of this plan, because it would require re-running Steps 2a–2b.

## References

- Collins GS, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. BMJ 2024;385:e078378.
- Harrell FE, Lee KL, Mark DB. Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. Stat Med 1996;15:361–387.
- McShane LM, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). J Natl Cancer Inst 2005;97:1180–1184.
- Pölsterl S. scikit-survival: a library for time-to-event analysis built on top of scikit-learn. J Mach Learn Res 2020;21(212):1–6.
- Schoenfeld DA. Sample-size formula for the proportional-hazards regression model. Biometrics 1983;39:499–503.
- Simon N, Friedman J, Hastie T, Tibshirani R. Regularization paths for Cox's proportional hazards model via coordinate descent. J Stat Softw 2011;39(5):1–13.
- Simon RM, Subramanian J, Li MC, Menezes S. Using cross-validation to evaluate predictive accuracy of survival risk classifiers based on high-dimensional data. Brief Bioinform 2011;12:203–214.
- Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. Stat Med 2011;30:1105–1117.
- Verweij PJM, van Houwelingen HC. Cross-validation in survival analysis. Stat Med 1993;12:2305–2314.
