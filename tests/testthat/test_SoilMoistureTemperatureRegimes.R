
test_that("SMTR", {
  #--- Expectations
  expected_STR_name <- "Cryic"
  tmp <- STR_names()
  expected_STR <- rep(0, length(tmp))
  names(expected_STR) <- tmp
  expected_STR[expected_STR_name] <- 1

  expected_SMR_name <- c("Ustic", "Typic-Tempustic")
  tmp <- c(SMR_names(), SMRq_names())
  expected_SMR <- rep(0, length(tmp))
  names(expected_SMR) <- tmp
  expected_SMR[expected_SMR_name] <- 1


  #--- Check with rSOILWAT2 example data
  sw_in <- rSOILWAT2::sw_exampleData
  sw_out <- rSOILWAT2::sw_exec(inputData = sw_in)

  SMTR1 <- calc_SMTRs(sim_in = sw_in, sim_out = sw_out)

  expect_true(SMTR1[["regimes_done"]])

  expect_true(all(colnames(SMTR1[["STR"]]) %in% STR_names()))

  expect_true(
    all(colnames(SMTR1[["SMR"]]) %in% c(SMR_names(), SMRq_names()))
  )

  expect_equal(SMTR1[["STR"]], expected_STR, ignore_attr = TRUE)

  expect_equal(SMTR1[["SMR"]], expected_SMR, ignore_attr = TRUE)


  #--- Check with insufficient soil layers (add and recalculate)
  sw_in2 <- sw_in

  # Dissolve soil layers
  xsoils <- rSOILWAT2::swSoils_Layers(sw_in2)
  depths <- xsoils[, "depth_cm"]
  colnames(xsoils) <- sapply(
    strsplit(colnames(xsoils), split = "_", fixed = TRUE),
    function(x) x[1]
  )

  tmp <- stats::reshape(
    data = data.frame(id = 1, xsoils),
    timevar = "depth",
    idvar = "id",
    direction = "wide"
  )

  xsoils_wide <- tmp[, -1]

  colnames(xsoils_wide) <- paste0(
    rep(colnames(xsoils[, -1]), length(depths)),
    "_L",
    rep(seq_along(depths), each = ncol(xsoils) - 1)
  )

  xsoils2_wide <- rSW2data::update_soil_profile(
    soil_layers = t(depths),
    requested_soil_layers = c(5, 20, 60, 80, 100),
    soil_data = xsoils_wide,
    vars_exhaust = c(
      "EvapBareSoil",
      "transpGrass", "transpShrub", "transpTree", "transpForb",
      "impermeability"
    ),
    keep_prev_soildepth = TRUE,
    keep_prev_soillayers = FALSE
  )

  tmp <- stats::reshape(
    data = data.frame(id = 1, xsoils2_wide[["soil_data"]]),
    varying = 1 + seq_len(ncol(xsoils2_wide[["soil_data"]])),
    sep = "_L",
    direction = "long"
  )

  tmp <- cbind(depth = xsoils2_wide[["soil_layers"]][1, ], tmp[, - (1:2)])
  xsoils2 <- tmp[complete.cases(tmp), ]

  colnames(xsoils2) <- colnames(rSOILWAT2::swSoils_Layers(sw_in2))
  rSOILWAT2::swSoils_Layers(sw_in2) <- data.matrix(xsoils2)

  #---
  sw_out2 <- rSOILWAT2::sw_exec(inputData = sw_in2)

  expect_output(
    SMTR2 <- calc_SMTRs(sim_in = sw_in2, sim_out = sw_out2, verbose = TRUE)
  )

  for (kn in names(SMTR1)) {
    if (kn == "cond_annual") {
      # Small deviations in annual values due to interpolation of soil layers
      expect_equal(SMTR2[[kn]], SMTR1[[kn]], tolerance = 0.1)
    } else {
      # do not affect outcomes
      expect_identical(SMTR2[[kn]], SMTR1[[kn]])
    }
  }
})
