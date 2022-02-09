
# Inputs
GISSM_tag <- "debug_GISSM_output"

test_that("GISSM", {
  #--- Prepare simulation inputs ------
  sw_in <- rSOILWAT2::sw_exampleData
  xlat <- rSOILWAT2::swSite_IntrinsicSiteParams(sw_in)[["Latitude"]]
  slyrs <- rSOILWAT2::swSoils_Layers(sw_in)[, 1]
  years <- rSOILWAT2::swYears_StartYear(sw_in):rSOILWAT2::swYears_EndYear(sw_in)
  simTime1 <- rSW2data::setup_time_simulation_run(
    sim_time = list(
      spinup_N = 0,
      startyr = years[1],
      endyr = years[length(years)]
    )
  )
  simTime2 <- rSW2data::simTiming_ForEachUsedTimeUnit(
    useyrs = simTime1[["useyrs"]],
    sim_tscales = "daily",
    latitude = xlat,
    account_NorthSouth = TRUE
  )


  #--- Prepare simulation outputs ------
  res <- rSOILWAT2::sw_exec(inputData = sw_in)
  has_Ts <-
    rSOILWAT2::swSite_SoilTemperatureFlag(sw_in) &&
    !rSOILWAT2::has_soilTemp_failed()

  dt <- "Day"
  tmp_swp <- slot(slot(res, rSW2_glovars[["swof"]]["sw_swp"]), dt)
  tmp_snow <- slot(slot(res, rSW2_glovars[["swof"]]["sw_snow"]), dt)
  tmp_airtemp <- slot(slot(res, rSW2_glovars[["swof"]]["sw_temp"]), dt)
  tmp_soiltemp <- slot(slot(res, rSW2_glovars[["swof"]]["sw_soiltemp"]), dt)
  xsim <- list(
    SWP_MPa = -0.1 * tmp_swp[, - (1:2), drop = FALSE],
    Snowpack_SWE_mm = 10 * tmp_snow[, "snowpackWaterEquivalent_cm"],
    air_Tmin_C = tmp_airtemp[, "min_C"],
    air_Tmean_C = tmp_airtemp[, "avg_C"],
    air_Tmax_C = tmp_airtemp[, "max_C"],
    # Using daily mean soil temperature in the absence of daily min/max
    shallowsoil_Tmin_C = tmp_soiltemp[, "Lyr_1"],
    shallowsoil_Tmean_C = tmp_soiltemp[, "Lyr_1"],
    shallowsoil_Tmax_C = tmp_soiltemp[, "Lyr_1"]
  )


  #--- Example 1: use rSOILWAT2 output directly to run GISSM ------
  GISSM_r1 <- suppressWarnings(calc_GISSM(
    x = res,
    soillayer_depths_cm = slyrs,
    site_latitude = xlat,
    has_soil_temperature = has_Ts
  ))


  #--- Example 2: use list of daily values to run GISSM ------
  #--- Pass time information as `years`
  GISSM_r2a <- calc_GISSM(
    x = xsim,
    soillayer_depths_cm = slyrs,
    site_latitude = xlat,
    years = years
  )

  # Checks
  expect_equal(GISSM_r1, GISSM_r2a)
  expect_named(GISSM_r1, "outcome")
  expect_equal(dim(GISSM_r1[["outcome"]]), c(length(years) - 1, 3))


  #--- Pass time information as `simTime1` instead of `years`
  GISSM_r2b <- calc_GISSM(
    x = xsim,
    soillayer_depths_cm = slyrs,
    site_latitude = xlat,
    simTime1 = simTime1
  )

  # Checks
  expect_equal(GISSM_r2a, GISSM_r2b)


  #--- Pass time information as `simTime1` and `simTime2`
  GISSM_r2c <- calc_GISSM(
    x = xsim,
    soillayer_depths_cm = slyrs,
    simTime1 = simTime1,
    simTime2 = simTime2
  )

  # Checks
  expect_equal(GISSM_r2a, GISSM_r2c)


  #--- Example 3: additional elements in the returned list ------
  GISSM_r3 <- suppressWarnings(calc_GISSM(
    x = res,
    soillayer_depths_cm = slyrs,
    site_latitude = xlat,
    has_soil_temperature = has_Ts,
    debug_output = 1
  ))

  expect_equal(GISSM_r1[["outcome"]], GISSM_r3[["outcome"]])
  expect_length(GISSM_r3, 8)


  #--- Example 4: additional output written to spreadsheet and pdf figure ------
  GISSM_r4 <- suppressWarnings(calc_GISSM(
    x = res,
    soillayer_depths_cm = slyrs,
    site_latitude = xlat,
    has_soil_temperature = has_Ts,
    debug_output = 2,
    path = ".",
    filename_tag = GISSM_tag
  ))

  expect_equal(GISSM_r1[["outcome"]], GISSM_r4[["outcome"]])
  expect_length(GISSM_r4, 8)
  expect_true(file.exists(paste0(GISSM_tag, ".csv")))
  expect_true(file.exists(paste0(GISSM_tag, ".pdf")))

  # clean up
  unlink(paste0(GISSM_tag, ".csv"))
  unlink(paste0(GISSM_tag, ".pdf"))


  #--- Example 5: GISSM start is on a leap or non-leap year ------
  for (sy in 1980:1981) {

    sw_in5 <- rSOILWAT2::sw_exampleData
    rSOILWAT2::swWeather_FirstYearHistorical(sw_in5) <- -1
    rSOILWAT2::swYears_StartYear(sw_in5) <- sy
    rSOILWAT2::swYears_EndYear(sw_in5) <- 2010

    res5 <- rSOILWAT2::sw_exec(inputData = sw_in5)

    GISSM_r5 <- suppressWarnings(calc_GISSM(
      x = res5,
      soillayer_depths_cm = rSOILWAT2::swSoils_Layers(sw_in5)[, 1],
      site_latitude =
        rSOILWAT2::swSite_IntrinsicSiteParams(sw_in5)[["Latitude"]],
      has_soil_temperature =
        rSOILWAT2::swSite_SoilTemperatureFlag(sw_in5) &&
        !rSOILWAT2::has_soilTemp_failed()
    ))

    expect_named(GISSM_r5, "outcome")
    expect_equal(dim(GISSM_r5[["outcome"]]), c(2010 - sy, 3))
  }


  #--- Example 6: Mis-specified time information
  bad_years <- years[-1]

  bad_simTime1 <- rSW2data::setup_time_simulation_run(
    sim_time = list(
      spinup_N = 0,
      startyr = years[10],
      endyr = years[length(years)]
    )
  )

  bad_simTime2 <- rSW2data::simTiming_ForEachUsedTimeUnit(
    useyrs = simTime1[["useyrs"]][-1],
    sim_tscales = "daily",
    latitude = xlat,
    account_NorthSouth = TRUE
  )

  #--- Disagreement between `years` and rSOILWAT2 simulation object
  expect_error(
    calc_GISSM(
      x = res,
      soillayer_depths_cm = slyrs,
      site_latitude = xlat,
      years = bad_years
    ),
    regexp = "Values of `years` disagree with content of `x`"
  )


  #--- Disagreement between `years` and `simTime1`
  expect_error(
    calc_GISSM(
      x = res,
      soillayer_depths_cm = slyrs,
      site_latitude = xlat,
      simTime1 = bad_simTime1
    ),
    regexp = "Values of `simTime1` and `years` are in disagreement"
  )

  expect_error(
    calc_GISSM(
      x = xsim,
      soillayer_depths_cm = slyrs,
      site_latitude = xlat,
      years = years,
      simTime1 = bad_simTime1
    ),
    regexp = "Values of `simTime1` and `years` are in disagreement"
  )

  #--- Disagreement between `simTime1` and `simTime2`
  expect_error(
    calc_GISSM(
      x = res,
      soillayer_depths_cm = slyrs,
      site_latitude = xlat,
      simTime2 = bad_simTime2
    ),
    regexp = "Values of `simTime1` and `simTime2` are in disagreement"
  )

  expect_error(
    calc_GISSM(
      x = xsim,
      soillayer_depths_cm = slyrs,
      site_latitude = xlat,
      years = years,
      simTime2 = bad_simTime2
    ),
    regexp = "Values of `simTime1` and `simTime2` are in disagreement"
  )

})
