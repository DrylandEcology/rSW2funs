context("GISSM")

# Inputs
GISSM_tag <- "debug_GISSM_output"

test_that("GISSM", {
  sw_in <- rSOILWAT2::sw_exampleData
  res <- rSOILWAT2::sw_exec(inputData = sw_in)

  # Example 1: use rSOILWAT2 output directly to run GISSM
  GISSM_r1 <- suppressWarnings(calc_GISSM(
    x = res,
    soillayer_depths_cm = rSOILWAT2::swSoils_Layers(sw_in)[, 1],
    site_latitude = rSOILWAT2::swSite_IntrinsicSiteParams(sw_in)[["Latitude"]],
    has_soil_temperature =
      rSOILWAT2::swSite_SoilTemperatureFlag(sw_in) &&
      !rSOILWAT2::has_soilTemp_failed()
  ))

  # Example 2: use list of daily values to run GISSM
  dt <- "Day"
  tmp_swp <- slot(slot(res, rSW2_glovars[["swof"]]["sw_swp"]), dt)
  tmp_snow <- slot(slot(res, rSW2_glovars[["swof"]]["sw_snow"]), dt)
  tmp_airtemp <- slot(slot(res, rSW2_glovars[["swof"]]["sw_temp"]), dt)
  tmp_soiltemp <- slot(slot(res, rSW2_glovars[["swof"]]["sw_soiltemp"]), dt)
  years <- rSOILWAT2::swYears_StartYear(sw_in):rSOILWAT2::swYears_EndYear(sw_in)

  GISSM_r2 <- calc_GISSM(
    x = list(
      SWP_MPa = -1 / 10 * tmp_swp[, -(1:2), drop = FALSE],
      Snowpack_SWE_mm = 10 * tmp_snow[, "snowpackWaterEquivalent_cm"],
      air_Tmin_C = tmp_airtemp[, "min_C"],
      air_Tmean_C = tmp_airtemp[, "avg_C"],
      air_Tmax_C = tmp_airtemp[, "max_C"],
      # Using daily mean soil temperature in the absence of daily min/max
      shallowsoil_Tmin_C = tmp_soiltemp[, "Lyr_1"],
      shallowsoil_Tmean_C = tmp_soiltemp[, "Lyr_1"],
      shallowsoil_Tmax_C = tmp_soiltemp[, "Lyr_1"]
    ),
    soillayer_depths_cm = rSOILWAT2::swSoils_Layers(sw_in)[, 1],
    site_latitude = rSOILWAT2::swSite_IntrinsicSiteParams(sw_in)[["Latitude"]],
    years = years
  )

  # Checks
  expect_equal(GISSM_r1, GISSM_r2)
  expect_named(GISSM_r1, "outcome")
  expect_equal(dim(GISSM_r1[["outcome"]]), c(length(years) - 1, 3))


  # Example 3: additional elements in the returned list
  GISSM_r3 <- suppressWarnings(calc_GISSM(
    x = res,
    soillayer_depths_cm = rSOILWAT2::swSoils_Layers(sw_in)[, 1],
    site_latitude = rSOILWAT2::swSite_IntrinsicSiteParams(sw_in)[["Latitude"]],
    has_soil_temperature =
      rSOILWAT2::swSite_SoilTemperatureFlag(sw_in) &&
      !rSOILWAT2::has_soilTemp_failed(),
    debug_output = 1
  ))

  expect_equal(GISSM_r1[["outcome"]], GISSM_r3[["outcome"]])
  expect_length(GISSM_r3, 8)


  # Example 4: additional output written to spreadsheet and pdf figure
  GISSM_r4 <- suppressWarnings(calc_GISSM(
    x = res,
    soillayer_depths_cm = rSOILWAT2::swSoils_Layers(sw_in)[, 1],
    site_latitude = rSOILWAT2::swSite_IntrinsicSiteParams(sw_in)[["Latitude"]],
    has_soil_temperature =
      rSOILWAT2::swSite_SoilTemperatureFlag(sw_in) &&
      !rSOILWAT2::has_soilTemp_failed(),
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
})
