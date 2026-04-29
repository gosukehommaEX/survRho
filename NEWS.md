# survRho 0.1.0

Initial release accompanying the article Homma and Komukai (20XX),
*"When Does Unequal Randomization Shorten Survival Trials? Closed-Form
Bounds on the Required Accrual Rate"*.

## Features

- `calc_trial_duration()`: computes the required accrual duration and
  total trial duration under the Lachin-Foulkes (1986) sample size
  framework with Schoenfeld's (1981) required number of events.
- `calc_relative_efficiency()` / `calc_relative_efficiency_grid()`:
  evaluate the relative efficiency of trial duration under r:1
  randomization across a single scenario or a grid of scenarios.
- `calc_rho_bounds()` / `calc_rho_bounds_grid()`: compute the three
  closed-form bounds on the accrual rate multiplier rho* derived in the
  article, namely the universal lower bound, the tight lower bound, and
  the upper bound.
- `find_rho_for_target_re()` / `find_rho_for_target_re_grid()`: solve
  numerically for the rho* that achieves a target relative efficiency.
- S3 print methods for all main output classes.

## Reproducibility

- The script `inst/table_and_figure_paper.R` regenerates all figures and
  the table reported in the article.
