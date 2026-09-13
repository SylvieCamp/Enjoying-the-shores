# Enjoying the shores — analysis code (revised version)

Code for *Enjoying the shores: a meta-analysis of human impacts on coastal and marine
cultural ecosystem services* (Campagne et al.), as revised for Cell Reports Sustainability.

## What changed from the submitted version, and why

1. **Study-level clustering.** The 223 effect sizes come from 59 studies (median 1, maximum 26
   per study). The submitted ordinal models (`clm`) treated them as independent. The revised
   primary analysis uses cumulative link mixed models with a study-level random intercept
   (`clmm`), marginalised over the random effect so that estimates are population-averaged and
   comparable with the fixed-effect ones. Table S11 reports both, plus a study-aggregated check.
2. **Bayesian estimates within drivers.** Within individual drivers the maximum-likelihood
   random-effect variance is not identifiable (it reaches 104 on the logit scale for
   management). Driver-by-service estimates (Figure 2a, Table S13) are therefore fitted with
   `brms` under a half-normal(0,1) prior on the study SD.
3. **Sampling variance in the meta-analytical models.** The submitted line
   `Data$vi <- max(0,001,abs(Data$yi*0.10))` contained a comma where a decimal point was
   intended; R evaluates it as a single scalar recycled over every observation. Table S12 shows
   that the results are unchanged under four alternative specifications of `vi`.
4. **Neutral / ambiguous labels.** The submitted script labelled `B_null` (neutral) as `AMB`
   and `C_null2` (ambiguous) as `NEU`. Positive and negative probabilities were unaffected;
   Table S8 has been corrected.
5. **Publication-bias tests.** Funnel plots and Egger tests were stated but never run; the
   statement has been removed (see Methods 3.2 for why they are not applicable here).

## Reproduce

Input: `Raw_data_14-02-24.xlsx` exported as `raw.csv` (223 rows). Then, in order:

| script | produces | runtime |
|---|---|---|
| `01_reproduce_and_sensitivity.R` | reproduces the submitted estimates to 1e-13; `S11_final.csv`, `S12_final.csv` | ~3 min |
| `02_bayes_drivers.R` | Bayesian check of the between-study variance under three priors | ~6 min |
| `03_bayes_driver_by_service.R` | `bayes_blocks.csv` = Table S13 | ~10 min |
| `04_ambiguous_excluded.R` | sensitivity to dropping the three ambiguous outcomes | <1 min |
| `00_cells_fixed_effect.R` | per-cell fixed-effect estimates (open diamonds, Figure 2a) | ~1 min |
| `05_figure1.R`, `06_figure2.R` | `Figure1_final.*`, `Figure2_final.*` (icons in `icons/`); `Rscript -e 'COMPARE <- TRUE; source("05_figure1.R")'` gives `FigureS2_final.*` | <1 min |
| `07_table1.R` | `table1_counts.csv` = Table 1 | <1 min |

R 4.3; packages `ordinal`, `metafor`, `MASS`, `brms` (Stan), `ggplot2`, `dplyr`, `tidyr`,
`patchwork`, `png`. Seeds are fixed; Bayesian fits use 4 chains × 2000 post-warm-up draws.

## Verification against the submitted results

`01_reproduce_and_sensitivity.R` reproduces every value in `Results_Vote_counting_15-03-24.xlsx`
(ordinal probabilities, to machine precision) and `Resuts_Metaanalyse_14-02-24.Final.xlsx`
(22 meta-analytical fits, |Δβ| < 1e-13, identical k and significance) before any change of
specification is applied.
