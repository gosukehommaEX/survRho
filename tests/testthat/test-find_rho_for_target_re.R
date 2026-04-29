# Tests for find_rho_for_target_re()
# Verifies that the numerical solver reproduces the rho* values reported
# in Section 4 of Homma & Komukai (20XX) and that rho* always falls within
# the analytical bounds.

test_that("rho* reproduces Table 1 (Phase 2 oncology, r = 4)", {
  # Phase 2 oncology at r = 4: paper reports rho* = 1.760
  result <- find_rho_for_target_re(
    HR = 0.6, hazard_c = log(2) / 9,
    a_k = 0, omega_k = 25,
    f = 6, d_t = 0.05, d_c = 0.05,
    r = 4, alpha = 0.10, power = 0.8,
    metric = "tau", target_re = 1
  )
  expect_equal(round(result$rho, 3), 1.760)
})

test_that("rho* reproduces Table 1 (Phase 3 oncology, r = 4)", {
  # Phase 3 oncology at r = 4: paper reports rho* = 1.692
  result <- find_rho_for_target_re(
    HR = 0.7, hazard_c = log(2) / 18,
    a_k = 0, omega_k = 50,
    f = 12, d_t = 0.10, d_c = 0.10,
    r = 4, alpha = 0.025, power = 0.9,
    metric = "tau", target_re = 1
  )
  expect_equal(round(result$rho, 3), 1.692)
})

test_that("rho* reproduces Table 1 (CV outcome, r = 2 and r = 4)", {
  # CV outcome at r = 2: paper reports rho* = 1.161
  res_r2 <- find_rho_for_target_re(
    HR = 0.8, hazard_c = log(2) / 60,
    a_k = 0, omega_k = 100,
    f = 20, d_t = 0.15, d_c = 0.15,
    r = 2, alpha = 0.025, power = 0.9,
    metric = "tau", target_re = 1
  )
  expect_equal(round(res_r2$rho, 3), 1.161)

  # CV outcome at r = 4: paper reports rho* = 1.654
  res_r4 <- find_rho_for_target_re(
    HR = 0.8, hazard_c = log(2) / 60,
    a_k = 0, omega_k = 100,
    f = 20, d_t = 0.15, d_c = 0.15,
    r = 4, alpha = 0.025, power = 0.9,
    metric = "tau", target_re = 1
  )
  expect_equal(round(res_r4$rho, 3), 1.654)
})

test_that("rho* always lies within [rho_L^tight, rho_U]", {
  # Verify the three-level hierarchy at the numerical level.
  test_cases <- list(
    list(HR = 0.6, mst_c = 9, f = 6, d = 0.05, omega1 = 25,
         alpha = 0.10, power = 0.8, r = 4),
    list(HR = 0.7, mst_c = 18, f = 12, d = 0.10, omega1 = 50,
         alpha = 0.025, power = 0.9, r = 4),
    list(HR = 0.8, mst_c = 60, f = 20, d = 0.15, omega1 = 100,
         alpha = 0.025, power = 0.9, r = 4)
  )
  for (tc in test_cases) {
    rho_star <- find_rho_for_target_re(
      HR = tc$HR, hazard_c = log(2) / tc$mst_c,
      a_k = 0, omega_k = tc$omega1,
      f = tc$f, d_t = tc$d, d_c = tc$d,
      r = tc$r, alpha = tc$alpha, power = tc$power,
      metric = "tau", target_re = 1
    )$rho
    bnd <- calc_rho_bounds(
      HR = tc$HR, hazard_c = log(2) / tc$mst_c, d = tc$d,
      r = tc$r, time_unit = "months"
    )
    expect_lte(bnd$bounds_df$rho_lower_tight, rho_star + 1e-6)
    expect_lte(rho_star, bnd$bounds_df$rho_upper + 1e-6)
  }
})
