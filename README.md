# Log-Adjusted and Calibration-Assisted Correlation Inference

This repository contains the code, analysis extracts, simulation designs, numerical summaries, figures, and validation checks for the manuscript:

**Log-Adjusted and Calibration-Assisted Correlation Inference under Covariate Distortion and Measurement Error**  
Yisa Abolade and Yichuan Zhao

The project addresses latent Pearson correlation under two complementary measurement regimes.

1. **Reference-free LAR/ME-LAR.** Positive measurements are subject to covariate-dependent multiplicative distortion. On the logarithmic scale, the distortion is an additive location shift. A kernel conditional exponential-mean normalization defines the log-adjusted residual (LAR) estimator. ME-LAR then subtracts a known or externally estimated log-scale measurement-error covariance from the adjusted second moments.
2. **Collocation-assisted correlation estimation.** An independent collocation sample contains reference values together with multichannel measurements. Channel-specific calibration functions and the full post-calibration residual covariance are estimated from the collocation data, transported to the target sample, and used to correct the target proxy covariance.

One-sample and two-sample jackknife empirical likelihood (JEL) and adjusted JEL (AJEL) procedures are provided. The two-sample implementation separately jackknifes collocation and target observations and variance-matches their pseudo-values before empirical-likelihood inference.

## Software

```bash
python -m pip install -r requirements.txt
```

Run all commands from the repository root.

## Repository structure

- `R/methods.R` - LAR/ME-LAR, calibration, covariance correction, jackknife, JEL/AJEL, and influence-bootstrap functions.
- `python/lar_melar.py` - reference-free LAR/ME-LAR and one-sample JEL/AJEL implementation used by the reported Monte Carlo study.
- `R/epa_data.R` - EPA data loading, platform/channel definitions, and EPA-informed parameter extraction.
- `R/simulation.R` - collocation-assisted latent-variable and measurement-error simulation mechanisms.
- `data/` - analysis extracts from the public EPA/USFS burn-integrated chamber study.
- `simulation/lar_melar_scenarios.csv` - nine reference-free simulation configurations.
- `simulation/run_lar_melar.py` - reported reference-free Monte Carlo study.
- `simulation/scenarios.csv` - EPA-informed collocation-assisted simulation configurations.
- `simulation/run_scenario.R`, `run_all.R`, `summarize_results.R` - collocation simulation workflow.
- `analysis/epa_application.R` - full-sample and 2,000 repeated held-out-burn EPA analyses.
- `analysis/make_lar_figures.py` - reference-free figures.
- `analysis/make_figures.R` - collocation/EPA figures.
- `results/` - tracked numerical summaries reported in the manuscript and supplement.
- `figures/` - tracked manuscript and supplementary figures.
- `tests/validate_lar.py`, `tests/validate.R` - deterministic numerical checks.

## Data source

The EPA analysis uses the U.S. Environmental Protection Agency and U.S. Forest Service 2019 and 2021 Missoula Fire Sciences Laboratory burn-integrated sensor-performance dataset:

**2019_2021_chamber_study_burn_integrated_data**  
DOI: `10.23719/1532438`

The files under `data/` retain the variables required for the reported CO--NO2 dependence analysis. See `data/README.md`.

## Validation

Run the two deterministic validation suites:

```bash
python tests/validate_lar.py
Rscript tests/validate.R
```

## Reference-free LAR/ME-LAR study

A short computational check can be run with:

```bash
python simulation/run_lar_melar.py --reps 20 --output results/work/lar_melar_check.csv
```

The reported study uses 1,000 Monte Carlo replications per configuration:

```bash
python simulation/run_lar_melar.py --reps 1000
python analysis/make_lar_figures.py
```

Tracked outputs:

```text
results/simulation/lar_melar_summary.csv
figures/simulation_lar_melar_bias.png
figures/simulation_lar_melar_coverage.png
```

## EPA application

```bash
Rscript analysis/epa_application.R
Rscript analysis/make_figures.R
```

Tracked outputs:

```text
results/epa/full_sample_results.csv
results/epa/heldout_validation_summary.csv
figures/epa_correlation_validation.png
```

## EPA-informed collocation simulation

A short scenario check:

```bash
Rscript simulation/run_scenario.R P01 10 1 49
```

The reported study uses 1,000 Monte Carlo replications per scenario and 499 influence-bootstrap draws:

```bash
Rscript simulation/run_all.R 1000 499
Rscript analysis/make_figures.R
```

Tracked summaries:

```text
results/simulation/primary_summary.csv
results/simulation/robustness_summary.csv
```
