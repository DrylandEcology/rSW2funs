
# Inputs
Nd <- 365
Nl <- 10
Nl2 <- round(Nl / 2)
Nl3 <- round(Nl / 3)
cond1 <- matrix(FALSE, nrow = Nd, ncol = Nl)
cond2 <- matrix(TRUE, nrow = Nd, ncol = Nl)
cond3 <- cbind(matrix(TRUE, nrow = Nd, ncol = Nl2),
              matrix(FALSE, nrow = Nd, ncol = Nl2))
cond4 <- cbind(matrix(TRUE, nrow = Nd, ncol = Nl3),
              matrix(FALSE, nrow = Nd, ncol = Nl3),
              matrix(TRUE, nrow = Nd, ncol = Nl3))



test_that("GISSM_get_KilledBySoilLayers", {

  expect_identical(GISSM_get_KilledBySoilLayers(NA, cond1), NA)
  expect_false(GISSM_get_KilledBySoilLayers(Nl, cond1))
  expect_true(GISSM_get_KilledBySoilLayers(Nl, cond2))
  expect_true(GISSM_get_KilledBySoilLayers(Nl2, cond3))
  expect_false(GISSM_get_KilledBySoilLayers(2 * Nl2, cond3))
  expect_true(GISSM_get_KilledBySoilLayers(Nl3, cond4))
  expect_false(GISSM_get_KilledBySoilLayers(2 * Nl3, cond4))
  expect_false(GISSM_get_KilledBySoilLayers(3 * Nl3, cond4))

  expect_identical(
    GISSM_get_KilledBySoilLayers(rep(NA, Nd), cond1), rep(NA, Nd)
  )
  expect_identical(
    GISSM_get_KilledBySoilLayers(rep(10, Nd), cond1), rep(FALSE, Nd)
  )
  expect_identical(
    GISSM_get_KilledBySoilLayers(rep(10, Nd), cond2), rep(TRUE, Nd)
  )
  expect_identical(
    GISSM_get_KilledBySoilLayers(rep(Nl2, Nd), cond3), rep(TRUE, Nd)
  )
  expect_identical(
    GISSM_get_KilledBySoilLayers(rep(2 * Nl2, Nd), cond3), rep(FALSE, Nd)
  )
  expect_identical(
    GISSM_get_KilledBySoilLayers(rep(Nl3, Nd), cond4), rep(TRUE, Nd)
  )
  expect_identical(
    GISSM_get_KilledBySoilLayers(rep(2 * Nl3, Nd), cond4), rep(FALSE, Nd)
  )
  expect_identical(
    GISSM_get_KilledBySoilLayers(rep(3 * Nl3, Nd), cond4), rep(FALSE, Nd)
  )

  #--- Errors
  # relevantLayers: too long
  expect_error(GISSM_get_KilledBySoilLayers(rep(NA, Nd + 1), cond1))

  # relevantLayers: too large values
  expect_error(GISSM_get_KilledBySoilLayers(Nl + 1, cond1))

  # relevantLayers: negative values
  expect_error(GISSM_get_KilledBySoilLayers(-1, cond1))
})
