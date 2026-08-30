# Simulation studies

The repository contains two simulation families.

## Reference-free LAR/ME-LAR

`lar_melar_scenarios.csv` defines nine log-scale distortion settings covering linear, quadratic, and periodic distortion at sample sizes 40, 80, and 160. Run:

```bash
python simulation/run_lar_melar.py --reps 1000
```

## EPA-informed collocation-assisted estimator

`scenarios.csv` defines the primary and additional configurations used for the multichannel collocation estimator.

- `EPA_primary`: A10 sample-size settings P00--P04.
- `EPA_replication`: Air Pro setting P05.
- `rho_grid`: target-correlation configurations.
- `error_strength`: measurement-contamination configurations.
- `error_structure`: covariance-structure configurations.
- `distribution`: heavy-tailed and skewed latent distributions.
- `misspecification`: varying-gain and nonlinear-additive calibration configurations.
- `transport`: target mechanism shifts toward the 2019 A10 mechanism.

Run one scenario:

```bash
Rscript simulation/run_scenario.R P01 10 1 49
```

Run the reported collocation simulation:

```bash
Rscript simulation/run_all.R 1000 499
```
