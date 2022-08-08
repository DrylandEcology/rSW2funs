
test_that("SMTR", {
  compare_smtrs <- function(target, current, tol_cond_annual = 0.1) {
    for (kn in names(target)) {
      if (kn == "cond_annual") {
        # Small deviations in annual values
        expect_equal(
          target[[!!kn]],
          current[[!!kn]],
          tolerance = tol_cond_annual
        )
      } else {
        expect_identical(target[[!!kn]], current[[!!kn]])
      }
    }
  }

  #--- Expectations ------
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


  #--- Check with rSOILWAT2 example data ------
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



  # Compare against latest rSOILWAT2 version that changed STMR output ------
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


  #--- Check with insufficient soil layers (add and recalculate) ------
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


  #--- Expect warning about additional soil layers
  expect_output(
    SMTR2 <- calc_SMTRs(
      sim_in = sw_in2,
      sim_out = rSOILWAT2::sw_exec(inputData = sw_in2),
      verbose = TRUE
    )
  )


  # Expect similar values to full set of soil layers
  # * Small deviations in annual values due to interpolation of soil layers
  compare_smtrs(SMTR2, SMTR1)


  #--- Check different SWRC/PDF options (if available) ------
  if (getNamespaceVersion("rSOILWAT2") >= as.numeric_version("6.0.0")) {

    #--- Set PDF from (default) Cosby1984AndOthers to Cosby1984
    # (avoid swc-sat complications for tests with unset `pdf_name`)
    # (default `has_swrcp` is FALSE)
    sw_in3 <- sw_in2
    rSOILWAT2::swSite_SWRCflags(sw_in3)[["pdf_name"]] <- "Cosby1984"
    rSOILWAT2::swSite_hasSWRCp(sw_in3) <- FALSE


    # Reference output (that requires additional soil layers) for next tests
    SMTR3 <- calc_SMTRs(
      sim_in = sw_in3,
      sim_out = rSOILWAT2::sw_exec(inputData = sw_in3)
    )

    # Expect similar values to default PDF of "Cosby1984AndOthers"
    # * Small deviations in annual values due to different swc-sat
    compare_smtrs(SMTR3, SMTR2)


    #--- Prepare test for different PDF options
    sw_in_swrc <- sw_in3
    soils <- rSOILWAT2::swSoils_Layers(sw_in_swrc)

    swrcp <- rSOILWAT2::pdf_estimate(
      sand = soils[, "sand_frac"],
      clay = soils[, "clay_frac"],
      fcoarse = soils[, "gravel_content"],
      bdensity = soils[, "bulkDensity_g/cm^3"],
      pdf_name = rSOILWAT2::swSite_SWRCflags(sw_in_swrc)[["pdf_name"]]
    )


    #--- Check active PDF (default) with swrcp (has_swrcp, non-default) ------
    # (`calc_SMTRs()` re-estimates SWRCp for additional soil layers)
    rSOILWAT2::swSite_hasSWRCp(sw_in_swrc) <- TRUE
    rSOILWAT2::swSoils_SWRCp(sw_in_swrc) <- swrcp

    # Expect identical to `SMTR3`
    expect_equal(
      calc_SMTRs(
        sim_in = sw_in_swrc,
        sim_out = rSOILWAT2::sw_exec(inputData = sw_in_swrc)
      ),
      SMTR3,
      tolerance = sqrt(.Machine[["double.eps"]])
    )


    #--- Check not-implemented PDF but with swrcp/has_swrcp ------
    # (`calc_SMTRs()` cannot re-estimate SWRCp for additional soil layers
    # and interpolates instead)
    rSOILWAT2::swSite_SWRCflags(sw_in_swrc)[["pdf_name"]] <-
      "IsCosby1984ButDontUseIt"
    rSOILWAT2::swSite_hasSWRCp(sw_in_swrc) <- TRUE
    rSOILWAT2::swSoils_SWRCp(sw_in_swrc) <- swrcp


    # Expect similar but not identical to `SMTR3`
    # * Small deviations in annual values due to interpolation of soil layers
    compare_smtrs(
      calc_SMTRs(
        sim_in = sw_in_swrc,
        sim_out = rSOILWAT2::sw_exec(inputData = sw_in_swrc)
      ),
      SMTR3,
      tol_cond_annual = 0.01
    )
  }
})
