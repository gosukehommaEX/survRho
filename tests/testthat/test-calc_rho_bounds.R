# Tests for calc_rho_bounds()
# Verifies that the closed-form bounds match the numerical values reported
# in Table 1 of Homma & Komukai (20XX).

test_that("Universal lower bound rho_L^univ depends only on r", {
  # rho_L^univ = (1 + r)^2 / (4 * r)
  for (rval in c(2, 3, 4)) {
    bnd <- calc_rho_bounds(
      HR = 0.7, hazard_c = log(2) / 18, d = 0.10,
      r = rval, time_unit = "months"
    )
    expected <- (1 + rval)^2 / (4 * rval)
    expect_equal(bnd$bounds_df$rho_lower_loose, expected, tolerance = 1e-10)
  }
})

test_that("Closed-form bounds reproduce Table 1 (Phase 3 oncology, r = 4)", {
  # Phase 3 oncology: HR = 0.7, MST_C = 18, d = 0.10
  # Paper values at r = 4: rho_L^univ = 1.562, rho_L^tight = 1.599,
  # rho_U = 1.748
  bnd <- calc_rho_bounds(
    HR = 0.7, hazard_c = log(2) / 18, d = 0.10,
    r = 4, time_unit = "months"
  )
  expect_equal(round(bnd$bounds_df$rho_lower_loose, 3), 1.562)
  expect_equal(round(bnd$bounds_df$rho_lower_tight, 3), 1.599)
  expect_equal(round(bnd$bounds_df$rho_upper, 3), 1.748)
})

test_that("Closed-form bounds reproduce Table 1 (CV outcome, r = 2 and 4)", {
  # CV outcome trial: HR = 0.8, MST_C = 60, d = 0.15
  # Paper values at r = 2: rho_L^tight = 1.149, rho_U = 1.168
  # Paper values at r = 4: rho_L^tight = 1.624, rho_U = 1.674
  bnd_r2 <- calc_rho_bounds(
    HR = 0.8, hazard_c = log(2) / 60, d = 0.15,
    r = 2, time_unit = "months"
  )
  expect_equal(round(bnd_r2$bounds_df$rho_lower_tight, 3), 1.149)
  expect_equal(round(bnd_r2$bounds_df$rho_upper, 3), 1.168)

  bnd_r4 <- calc_rho_bounds(
    HR = 0.8, hazard_c = log(2) / 60, d = 0.15,
    r = 4, time_unit = "months"
  )
  expect_equal(round(bnd_r4$bounds_df$rho_lower_tight, 3), 1.624)
  expect_equal(round(bnd_r4$bounds_df$rho_upper, 3), 1.674)
})

test_that("Three-level hierarchy holds: rho_L^univ <= rho_L^tight <= rho_U", {
  # Test across a small grid spanning the realistic design space.
  test_grid <- expand.grid(
    HR = c(0.5, 0.7, 0.9),
    mst_c = c(12, 24),
    d = c(0, 0.15, 0.3),
    r = c(2, 3, 4)
  )
  for (i in seq_len(nrow(test_grid))) {
    row <- test_grid[i, ]
    bnd <- calc_rho_bounds(
      HR = row$HR, hazard_c = log(2) / row$mst_c, d = row$d,
      r = row$r, time_unit = "months"
    )
    expect_lte(bnd$bounds_df$rho_lower_loose,
               bnd$bounds_df$rho_lower_tight + 1e-10)
    expect_lte(bnd$bounds_df$rho_lower_tight,
               bnd$bounds_df$rho_upper + 1e-10)
  }
})
