
test_that("SMTR", {
  #--- Expectations
  create_STR_expectation <- function(STRs) {
    tmp <- STR_names()
    stopifnot(STRs %in% tmp)
    expected_STR <- rep(0L, length(tmp))
    names(expected_STR) <- tmp
    expected_STR[STRs] <- 1L
    expected_STR
  }

  create_SMR_expectation <- function(SMRs) {
    tmp <- c(SMR_names(), SMRq_names())
    stopifnot(SMRs %in% tmp)
    expected_SMR <- rep(0L, length(tmp))
    names(expected_SMR) <- tmp
    expected_SMR[SMRs] <- 1L
    expected_SMR
  }


  #--- Check with rSOILWAT2 example data
  sw_in <- rSOILWAT2::sw_exampleData
  sw_out <- rSOILWAT2::sw_exec(inputData = sw_in)

  SMTR1 <- calc_SMTRs(sim_in = sw_in, sim_out = sw_out)

  expect_true(SMTR1[["regimes_done"]])

  expect_true(all(colnames(SMTR1[["STR"]]) %in% STR_names()))

  expect_true(
    all(colnames(SMTR1[["SMR"]]) %in% c(SMR_names(), SMRq_names()))
  )

  expect_true(all(SMTR1[["STR"]]) %in% c(0L, 1L))
  expect_true(all(SMTR1[["SMR"]]) %in% c(0L, 1L))



  # Compare against latest rSOILWAT2 version that changed STMR output
  tmpv <- c("5.0.0", "5.1.0")

  compv <- vapply(
    tmpv,
    function(v) {
      rSOILWAT2::check_version(sw_out, v, level = "minor")
    },
    FUN.VALUE = NA
  )

  eqv <- max(as.numeric_version(names(compv)[compv]))


  if (length(eqv) == 1) {
    if (eqv == "5.1.0") {
      # rSOILWAT2 since v5.1.0
      expected_STR <- create_STR_expectation("Cryic")
      expected_SMR <- create_SMR_expectation(c("Xeric", "Typic-Xeric"))

      } else if (eqv == "5.0.0") {
      # rSOILWAT2 since v5.0.0
      expected_STR <- create_STR_expectation("Cryic")
      expected_SMR <- create_SMR_expectation(c("Ustic", "Typic-Tempustic"))
    }

    expect_identical(SMTR1[["STR"]][1L, , drop = TRUE], expected_STR)
    expect_identical(SMTR1[["SMR"]][1L, , drop = TRUE], expected_SMR)

  } else {
    warning(
      "Test expectations for STR/SMR have not yet been implemented using ",
      "rSOILWAT2 v", rSOILWAT2::get_version(sw_out)
    )
  }


  #--- Check with insufficient soil layers (add and recalculate)
  sw_in2 <- sw_in

  # Dissolve soil layers
  xsoils <- rSOILWAT2::swSoils_Layers(sw_in2)
  depths <- xsoils[, "depth_cm"]
  colnames(xsoils) <- vapply(
    strsplit(colnames(xsoils), split = "_", fixed = TRUE),
    function(x) x[[1]],
    FUN.VALUE = NA_character_
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
