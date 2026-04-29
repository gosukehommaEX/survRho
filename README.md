# survRho

`survRho` provides closed-form expressions and analytical bounds for the
accrual rate multiplier $\rho^{*}$ that preserves the overall trial
duration when switching from $1:1$ to $r:1$ randomization in time-to-event
clinical trials. The methodology is developed under the Lachin-Foulkes
(1986) sample size framework combined with Schoenfeld's (1981) required
number of events, and is described in Homma and Komukai (20XX), *"When
Does Unequal Randomization Shorten Survival Trials? Closed-Form Bounds on
the Required Accrual Rate"*.

The package implements:

- the **universal lower bound** $\rho_{L}^{\rm univ} = (1 + r)^{2} / (4r)$,
  which depends only on the randomization ratio $r$;
- the **tight lower bound** $\rho_{L}^{\rm tight}$, which additionally
  incorporates the hazard ratio HR, the control hazard rate, and the
  dropout hazard;
- the **upper bound** $\rho_{U}$, which depends on $r$ and HR;
- numerical computation of the exact $\rho^{*}$ via root-finding.

## Installation

The package is currently available from GitHub:

```r
# install.packages("remotes")
remotes::install_github("gosukehommaEX/survRho")
```

## Quick start

### Example 1: Compute the three closed-form bounds on $\rho^{*}$

```r
library(survRho)

# Phase 3 oncology trial: HR = 0.7, control median survival 18 months,
# annual dropout probability 0.10, with r:1 randomization for r in {2, 3, 4}.
bounds <- calc_rho_bounds(
  HR       = 0.7,
  hazard_c = log(2) / 18,
  d        = 0.10,
  r        = c(2, 3, 4),
  time_unit = "months"
)
print(bounds)
```

### Example 2: Find the exact $\rho^{*}$ that preserves the trial duration

```r
library(survRho)

# Same trial as above, with omega1 = 50 patients per month under 1:1.
# Solve for the accrual rate multiplier rho* that gives the same total
# trial duration as 1:1 randomization, when switching to r = 4.
result <- find_rho_for_target_re(
  HR       = 0.7,
  hazard_c = log(2) / 18,
  a_k      = 0,
  omega_k  = 50,
  f        = 12,
  d_t      = 0.10,
  d_c      = 0.10,
  r        = 4,
  alpha    = 0.025,
  power    = 0.9,
  metric   = "tau",
  target_re = 1
)
print(result)
```

## Reproducing the figures and table in the article

The script that generates all figures and the table reported in Homma and
Komukai (20XX) is included in the package and can be located via:

```r
system.file("paper", "table_and_figure_paper.R", package = "survRho")
```

To regenerate the outputs locally, set the desired output directory and
source the script:

```r
library(survRho)
OUT_DIR <- "~/survRho_outputs"   # any path you prefer
source(system.file("paper", "table_and_figure_paper.R", package = "survRho"))
```

The script writes four figure files (EPS) and the application-section
table (TeX and CSV) to `OUT_DIR`. If `OUT_DIR` is not set, the script
defaults to a subdirectory `outputs/` of the current working directory.

## License

MIT License. See the `LICENSE` file for details.
