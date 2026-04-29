# Tests for calc_trial_duration()
# Verifies that the trial duration calculation reproduces values reported
# in Section 4 of Homma & Komukai (20XX).

test_that("Trial duration reproduces tau_1 in Table 1 (Phase 3 oncology)", {
  # Phase 3 oncology baseline (r = 1): tau_1 = 27.3 months
  td <- calc_trial_duration(
    HR = 0.7, hazard_c = log(2) / 18,
    a_k = 0, omega_k = 50,
    f = 12, d_t = 0.10, d_c = 0.10,
    r = 1, alpha = 0.025, power = 0.9,
    time_unit = "months", formula = "Schoenfeld"
  )
  expect_equal(round(td$tau, 1), 27.3)
})

test_that("Trial duration reproduces tau_1 in Table 1 (CV outcome)", {
  # CV outcome baseline (r = 1): tau_1 = 53.8 months
  td <- calc_trial_duration(
    HR = 0.8, hazard_c = log(2) / 60,
    a_k = 0, omega_k = 100,
    f = 20, d_t = 0.15, d_c = 0.15,
    r = 1, alpha = 0.025, power = 0.9,
    time_unit = "months", formula = "Schoenfeld"
  )
  expect_equal(round(td$tau, 1), 53.8)
})

test_that("Output has the documented columns (T/C naming)", {
  td <- calc_trial_duration(
    HR = 0.7, hazard_c = log(2) / 18,
    a_k = 0, omega_k = 50,
    f = 12, d_t = 0.10, d_c = 0.10,
    r = 2, alpha = 0.025, power = 0.9,
    time_unit = "months", formula = "Schoenfeld"
  )
  expected_cols <- c("HR", "mst_t", "mst_c", "hazard_t", "hazard_c",
                     "d_t", "d_c", "f", "r", "alpha", "power",
                     "E", "n_t", "n_c", "N", "a_K", "tau", "time_unit")
  expect_true(all(expected_cols %in% names(td)))
})

test_that("Sample size relationship n_t = r * n_c holds", {
  # Under r:1 randomization with continuous values, n_t / n_c should equal r.
  for (rval in c(2, 3, 4)) {
    td <- calc_trial_duration(
      HR = 0.7, hazard_c = log(2) / 18,
      a_k = 0, omega_k = 50,
      f = 12, d_t = 0.10, d_c = 0.10,
      r = rval, alpha = 0.025, power = 0.9,
      time_unit = "months", formula = "Schoenfeld"
    )
    expect_equal(td$n_t / td$n_c, rval, tolerance = 1e-8)
  }
})

test_that("Schoenfeld required events E_r reproduces Table 1 (Phase 3, r = 4)", {
  # Phase 3 oncology at r = 4: paper reports E_r = 517 (rounded up)
  td <- calc_trial_duration(
    HR = 0.7, hazard_c = log(2) / 18,
    a_k = 0, omega_k = 50,
    f = 12, d_t = 0.10, d_c = 0.10,
    r = 4, alpha = 0.025, power = 0.9,
    time_unit = "months", formula = "Schoenfeld"
  )
  expect_equal(ceiling(td$E), 517)
})
