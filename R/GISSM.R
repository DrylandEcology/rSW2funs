# ----------------------------
#------ GISSM functions ------
# Schlaepfer, D.R., Lauenroth, W.K. & Bradford, J.B. (2014). Modeling
# regeneration responses of big sagebrush (Artemisia tridentata) to abiotic
# conditions. Ecol Model, 286, 66-77.

# Function to convert soil depth to soil layer
SoilLayer_at_SoilDepth <- function(depth_cm, layers_depth) {
  pmax(
    1,
    pmin(
      length(layers_depth),
      1 + findInterval(depth_cm - 0.01, layers_depth)
    )
  )
}


# Function to calculate for each day of the year, duration in days of
# upcoming favorable conditions accounting for consequences_unfavorable = 0
# (if conditions become unfavorable, then restart the count), =1 (resume)
calc_DurationFavorableConds <- function(RYyear, consequences_unfavorable,
  Germination_WhileFavorable, RYyear_ForEachUsedDay) {

  index.year <- RYyear_ForEachUsedDay == RYyear
  conditions <- Germination_WhileFavorable[index.year]
  doys <- seq_len(sum(index.year))
  doys[!conditions] <- NA  #calculate only for favorable days
  out <- rep(NA, times = sum(index.year))

  if (consequences_unfavorable == 0) {
    # if conditions become unfavorable, then restart the count afterwards
    temp.rle <- rle(conditions)
    if (sum(!temp.rle$values) > 0) {
      temp.unfavorable_startdoy <- c((1 + c(0,
        # add starts for odd- and even-lengthed rle
        cumsum(temp.rle$lengths)))[!temp.rle$values], 1 + sum(index.year))

      temp.rle$values <- if (temp.rle$values[1]) {
        # first rle period is favorable
        rep(temp.unfavorable_startdoy, each = 2)
      } else {
        # first rle period is unfavorable
        rep(temp.unfavorable_startdoy[-1], each = 2)
      }
      temp.rle$values <- temp.rle$values[seq_along(temp.rle$lengths)]

    } else {
      # every day is favorable
      temp.rle$values <- length(conditions) + 1
    }

    # difference to next following start of a period of unfavorable conditions
    out <- inverse.rle(temp.rle) - doys

  } else if (consequences_unfavorable == 1) {
    # if conditions become unfavorable, then resume the count afterwards
    temp <- sum(conditions)
    count <- if (temp > 0) {
      temp:1
    } else {
      # every day is unfavorable
      vector("numeric", length = 0)
    }

    # sum of following favorable conditions in this year
    out <- stats::napredict(stats::na.action(stats::na.exclude(doys)), count)
  }

  out
}

# Based on the \var{NLR} model (equation 5) in Hardegree (2006) and modified
# by Schlaepfer et al. (2014) by making time to germinate dependent on
# mean January temperature and soil water potential
#
# @references Hardegree SP (2006) Predicting Germination Response to
#   Temperature. I. Cardinal-temperature Models and Subpopulation-specific
#   Regression. Annals of Botany, 97, 1115-1125.
get_modifiedHardegree2006NLR <- function(RYdoy, Estimate_TimeToGerminate,
  TmeanJan, a, b, c, d, k1_meanJanTemp, k2_meanJanTempXIncubationTemp,
  k3_IncubationSWP, Tgerm.year, SWPgerm.year, durations, rec.delta = 1,
  nrec.max = 10L) {

  for (nrec in seq_len(nrec.max)) {
    Estimate_TimeToGerminate <- prev_est_TimeToGerminate <-
      max(0, round(Estimate_TimeToGerminate))

    ids <- RYdoy:(RYdoy + Estimate_TimeToGerminate - 1)
    Tgerm <- mean(Tgerm.year[ids], na.rm = TRUE)
    SWPgerm <- mean(SWPgerm.year[ids], na.rm = TRUE)

    temp.c.lim <- - (Tgerm - b) * (d ^ 2 - 1) / d
    c <- if (c > 0) {
      if (c > temp.c.lim) c else {
        temp.c.lim + rSW2_glovars[["tol"]]
      }
    } else if (c < 0) {
      if (c < temp.c.lim) c else {
        temp.c.lim - rSW2_glovars[["tol"]]
      }
    }

    # NLR model (eq.5) in Hardegree SP (2006)
    temp <- a * exp(-0.693147181 / log(d) ^ 2 * log(1 + (Tgerm - b) *
        (d ^ 2 - 1) / (c * d)) ^ 2) # 0.693147181 is equal to log(2)

    # drs addition to time to germinate dependent on mean January temperature
    # and soil water potential
    temp <- 1 / temp +
      k1_meanJanTemp * TmeanJan +
      k2_meanJanTempXIncubationTemp * TmeanJan * Tgerm +
      k3_IncubationSWP * SWPgerm
    Estimate_TimeToGerminate <- max(1, round(temp))

    # break if convergence or not enough time in this year
    if (
      abs(Estimate_TimeToGerminate - prev_est_TimeToGerminate) <= rec.delta ||
      RYdoy + Estimate_TimeToGerminate - 1 > 365
    ) {
      break
    }
  }

  out <- if (nrec >= nrec.max) {
    round(mean(c(Estimate_TimeToGerminate, prev_est_TimeToGerminate)), 0)
  } else {
    Estimate_TimeToGerminate
  }

  # test whether enough time to germinate
  if (out <= durations[RYdoy] && RYdoy + out <= 365) out else NA
}

# Function to estimate time to germinate for each day of a given year and
# conditions (temperature, top soil \var{SWP})
#
# @param seed A seed set, \code{NULL}, or \code{NA}. \code{NA} will not affect
#  the state of the \var{RNG}; \code{NULL} will re-initialize the \var{RNG};
#  and all other values are passed to \code{\link{set.seed}}.
calc_TimeToGerminate <- function(RYyear, Germination_WhileFavorable,
  LengthDays_FavorableConditions, RYyear_ForEachUsedDay, soilTmeanSnow,
  swp_shallow, TmeanJan, params, seed = NA) {

  if (!is.na(seed)) set.seed(seed)
  runifs <- stats::runif(2)

  #values for current year
  index.year <- RYyear_ForEachUsedDay == RYyear
  conditions <- Germination_WhileFavorable[index.year]

  # determining time to germinate for every day
  a <- max(rSW2_glovars[["tol"]], params[["Hardegree_a"]])
  b <- params[["Hardegree_b"]]

  tmp <- if (params[["Hardegree_d"]] == 1) {
    if (runifs[1] > 0.5) {
      1 + rSW2_glovars[["tol"]]
    } else {
      1 - rSW2_glovars[["toln"]]
    }
  } else {
    params[["Hardegree_d"]]
  }
  d <- max(rSW2_glovars[["tol"]], tmp)

  tmp.c <- if (params[["Hardegree_c"]] != 0) {
    params[["Hardegree_c"]]
  } else {
    sign(runifs[2] - 0.5) * rSW2_glovars[["tol"]]
  }

  # consequences of unfavorable conditions coded in here
  TimeToGerminate.favorable <- unlist(lapply(which(conditions),
    get_modifiedHardegree2006NLR,
    Estimate_TimeToGerminate = 1,
    TmeanJan = TmeanJan,
    a = a,
    b = b,
    c = tmp.c,
    d = d,
    k1_meanJanTemp = params[["TimeToGerminate_k1_meanJanTemp"]],
    k2_meanJanTempXIncubationTemp =
      params[["TimeToGerminate_k2_meanJanTempXIncubationTemp"]],
    k3_IncubationSWP = params[["TimeToGerminate_k3_IncubationSWP"]],
    Tgerm.year = soilTmeanSnow[index.year],
    SWPgerm.year = swp_shallow[index.year],
    durations = LengthDays_FavorableConditions[index.year]
  ))

  res <- rep(NA, length(conditions))
  if (length(TimeToGerminate.favorable) > 0) {
    res[conditions] <- TimeToGerminate.favorable
  }

  res
}

do.vector <- function(kill.vector, max_time_to_kill) {
  doys <- seq_along(kill.vector)
  doys[!kill.vector] <- NA  #calculate only for kill days
  temp.rle <- rle(kill.vector)

  if (sum(!temp.rle$values) > 0) {
    temp.startdoy <- (1 + c(0, cumsum(temp.rle$lengths)))[!temp.rle$values]
    temp.rle$values <- if (temp.rle$values[1]) {
      rep(temp.startdoy, each = 2)
    } else {
      rep(temp.startdoy[-1], each = 2)
    }
    temp.rle$values <- temp.rle$values[seq_along(temp.rle$lengths)]

  } else {
    # every day is kill free
    temp.rle$values <- length(kill.vector) + 1
  }
  kill.durations <- inverse.rle(temp.rle) - doys
  mortality <- rep(FALSE, times = length(kill.vector))
  mortality[kill.durations > max_time_to_kill] <- TRUE

  mortality
}

# Function to calculate mortality under conditions and checks survival limit
calc_SeedlingMortality <- function(kill_conds,
  max_time_to_kill) {

  if (length(dim(kill_conds)) > 0) {
    # i.e., is.matrix, columns represent soil layers
    apply(kill_conds, 2, do.vector, max_time_to_kill)
  } else {
    do.vector(kill_conds, max_time_to_kill)
  }
}


# Function to calculate favorable conditions for seedling growth for each day
# of a given year
check_SuitableGrowthThisYear <- function(
  favorable_conditions, consequences_unfavorable) {

  out <- rep(NA, times = length(favorable_conditions))

  if (consequences_unfavorable == 0) {
    # if conditions become unfavorable, then stop growth for rest of season
    temp.rle <- rle(favorable_conditions)
    temp.firstFavorable.index <- which(temp.rle$values)[1]

    if (!is.na(temp.firstFavorable.index) &&
        temp.firstFavorable.index < length(temp.rle$values)) {

      temp <- (temp.firstFavorable.index + 1):length(temp.rle$values)
      temp.rle$values[temp] <- FALSE
      out <- inverse.rle(temp.rle)

    } else {
      # nothing changed, either because all days are either favorable or
      # unfavorable or because first favorable period is also the last in the
      # season
      out <- favorable_conditions
    }

  } else if (consequences_unfavorable == 1) {
    # if conditions become unfavorable, then resume growth afterwards
    out <- favorable_conditions
  }

  out
}


# Function to calculate rooting depth at given age
# Units: [age] = days, [P0, K, r] = mm
# @return A numeric vector of rooting depth in units of centimeters.
SeedlingRootingDepth <- function(age, P0, K, r) {
  depth <- K * P0 * exp(r * age) / (K + P0 * (exp(r * age) - 1))

  pmax(0, depth) / 10
}


get.DoyAtLevel <- function(x, level) {
  which(x == level & x > 0)
}

get.DoyMostFrequentSuccesses <- function(doys, data) {
  # must return one of the values because the quantiles are compared against
  # the values in function 'get.DoyAtLevel'
  res1.max <- sapply(1:2, function(x)
    stats::quantile(doys[doys[, x] > 0, x], probs = c(0.1, 1), type = 3))
  germ.doy <- if (all(!data[, 1])) {
    # no successful germination
    list(NA, NA)
  } else {
    lapply(1:2, function(x) get.DoyAtLevel(doys[, 1], res1.max[x, 1]))
  }
  sling.doy <- if (all(!data[, 2])) {
    # no successful seedlings
    list(NA, NA)
  } else {
    lapply(1:2, function(x) get.DoyAtLevel(doys[, 2], res1.max[x, 2]))
  }
  res1.max <- list(germ.doy, sling.doy)

  unlist(lapply(res1.max, function(x)
    c(min(x[[1]]), stats::median(x[[2]]), max(x[[1]]))))
}



#' Default parameter values of GISSM for big sagebrush
default_parameters_GISSM_bigsagebrush <- function() {
  list(
    Doy_SeedDispersalStart0 = 324.5569743,
    SeedDispersalStart_DependencyOnMeanTempJanuary = 2.039915438,
    GerminationPeriods_0ResetOr1Resume = 1,
    Temp_ExperiencedUnderneathSnowcover = 3.406591694,
    Temp_MaximumForGermination = 43.84200172,
    Temp_MinimumForGermination = 3.620305876,
    SWP_MinimumForGermination = -0.447054684,
    SeedlingGrowth_0StopOr1Resume = 1,
    SWE_MaximumForSeedlingGrowth = 0,
    Days_SnowCover_MaximumForSeedlingSurvival = 31,
    Temp_MinimumForSeedlingGrowth = -2.616505128,
    Temp_MaximumForSeedlingGrowth = 34.46478864,
    Temp_MinimumForSeedlingSurvival = -9.25483659,
    Temp_MaximumForSeedlingSurvival = 34.46478864,
    SWP_ChronicMaximumForSeedlingSurvival = -0.03333,
    Days_ChronicMaximumForSeedlingSurvival = 56,
    SWP_ChronicMinimumForSeedlingSurvival = -2.259034726,
    Days_ChronicMinimumForSeedlingSurvival = 49,
    SWP_AcuteMinimumForSeedlingSurvival = -3.278547773,
    SoilDepth_RelevantToGermination = 3,
    Seedling_SoilDepth.PO = 74,
    Seedling_SoilDepth.K = 1765,
    Seedling_SoilDepth.r = 0.189413231,
    Hardegree_a = 0.649614337,
    Hardegree_b = 13.7107755,
    Hardegree_c = -116.2694843,
    Hardegree_d = 0.365901519,
    TimeToGerminate_k1_meanJanTemp = -0.395946153,
    TimeToGerminate_k2_meanJanTempXIncubationTemp = 0.267351215,
    TimeToGerminate_k3_IncubationSWP = -3.538845403
  )
}

#' Parameter values of GISSM for big sagebrush
#'
#' @param ... Named GISSM parameters and their values.
#' @return A list of default and adjusted parameter values for GISSM.
#'
#' @references Schlaepfer, D.R., Lauenroth, W.K. & Bradford, J.B. (2014).
#'   Modeling regeneration responses of big sagebrush (Artemisia tridentata)
#'   to abiotic conditions. Ecol Model, 286, 66-77.
#'
#' @examples
#' # Default values
#' x <- parameters_GISSM_bigsagebrush()
#'
#' # Fix \var{Doy_SeedDispersalStart0} to 325 and use default values otherwise
#' x <- parameters_GISSM_bigsagebrush(Doy_SeedDispersalStart0 = 325)
#'
#' @export
parameters_GISSM_bigsagebrush <- function(...) {
  dots <- list(...)
  nd <- names(dots)

  x <- default_parameters_GISSM_bigsagebrush()
  nx <- names(x)

  used <- intersect(nd, nx)

  for (k in used) {
    x[[k]] <- dots[[k]]
  }

  badp <- setdiff(nd, nx)
  if (length(badp)) {
    warning(
      "Arguments ", paste(shQuote(badp), collapse = ", "),
      " are not GISSM parameters; they are ignored."
    )
  }

  x
}

#' The germination and individual seedling survival model \var{GISSM}
#'
#' GISSM represents the frequency of years when big sagebrush seeds germinate
#' and seedlings survive in undisturbed natural vegetation
#' (Schlaepfer et al. 2014).
#'
#' @param x A named list or an object of
#'   \pkg{rSOILWAT2} class \code{\linkS4class{swOutput}} with daily output.
#'   If \code{x} is a named list, then it must contain appropriate content for
#'   \var{SWP_MPa}, \var{Snowpack_SWE_mm}, \var{air_Tmin_C}, \var{air_Tmax_C},
#'   \var{air_Tmean_C}, \var{shallowsoil_Tmin_C}, \var{shallowsoil_Tmean_C},
#'   \var{shallowsoil_Tmax_C}. See examples.
#' @param soillayer_depths_cm A numeric vector. The lower bounds of simulated
#'   soil layers.
#' @param params A named list. See \code{\link{parameters_GISSM_bigsagebrush}}.
#' @param site_latitude A numeric value.
#'   Required if \code{simTime2} is \code{NULL}.
#' @param has_soil_temperature A logical value or \code{NULL}. Optional
#'   information whether or not object \code{x} contains (good)
#'   soil temperature values.
#' @param years A numeric vector or \code{NULL}. The sequence of simulated
#'   calendar years; required only
#'   if \code{x} is not of class \code{\linkS4class{swOutput}} and
#'   if \code{simTime1} and/or \code{simTime1} is \code{NULL}.
#' @param simTime1 A named list or \code{NULL}.
#'   See \code{\link[rSOILWAT2]{setup_time_simulation_run}}.
#' @param simTime2 A named list or \code{NULL}.
#'   See \code{\link[rSOILWAT2]{simTiming_ForEachUsedTimeUnit}}.
#' @param debug_output An integer value. Level of additional outputs.
#'   If \code{0}, then standard variables are returned.
#'   If \code{1} or \code{2}, then additional elements are included in the
#'   returned item; see \var{Value}.
#'   If \code{2}, then additional output for debugging purposes is
#'   written to a \var{csv} spreadsheet and a \var{pdf} figure is created.
#' @param path A character string. The path to the directory where additional
#'   output should be written to disk.
#' @param filename_tag A character string. File name without extension used
#'   for writing additional output.
#'
#' @return A named list with one element: \describe{
#'     \item{outcome}{A \var{data.frame} tabulating success/failure
#'     for \var{Germination_Emergence} and for \var{SeedlingSurvival_1stSeason}
#'     for each regeneration year.}
#'   }
#'
#'   If \code{debug_output} is \code{2}, then seven additional elements
#'   are provided for debugging purposes: \describe{
#'     \item{successes_days}{Number of days per year with positive outcomes.}
#'     \item{success_mostfrequent_doy}{
#'       Seasonal timing of most frequent positive outcomes.
#'     }
#'     \item{time_to_germinate_days}{Time to germinate.}
#'     \item{nogermination_days}{Days without germination.}
#'     \item{nogermination_periods_yrs}{Consecutive years without germination.}
#'     \item{noseedlings_periods_yrs}{
#'       Consecutive years without seedling recruitment.
#'     }
#'     \item{mortality_causes}{Causes of seedling mortality.}
#'   }
#'
#' @references Schlaepfer, D.R., Lauenroth, W.K. & Bradford, J.B. (2014).
#'   Modeling regeneration responses of big sagebrush (Artemisia tridentata)
#'   to abiotic conditions. Ecol Model, 286, 66-77.
#'
#' @examples
#' sw_in <- rSOILWAT2::sw_exampleData
#' res <- rSOILWAT2::sw_exec(inputData = sw_in)
#'
#' # Example 1: use rSOILWAT2 output directly to run GISSM
#' GISSM_r1 <- calc_GISSM(
#'   x = res,
#'   soillayer_depths_cm = rSOILWAT2::swSoils_Layers(sw_in)[, 1],
#'   site_latitude = rSOILWAT2::swSite_IntrinsicSiteParams(sw_in)[["Latitude"]],
#'   has_soil_temperature =
#'     rSOILWAT2::swSite_SoilTemperatureFlag(sw_in) &&
#'     !rSOILWAT2::has_soilTemp_failed()
#' )
#'
#'
#' # Example 2: use list of daily values to run GISSM
#' #   populate daily values from rSOILWAT2 output as here or from other model
#' dt <- "Day"
#' tmp_swp <- slot(slot(res, "SWPMATRIC"), dt)
#' tmp_snow <- slot(slot(res, "SNOWPACK"), dt)
#' tmp_airtemp <- slot(slot(res, "TEMP"), dt)
#' tmp_soiltemp <- slot(slot(res, "SOILTEMP"), dt)
#'
#' GISSM_r2 <- calc_GISSM(
#'   x = list(
#'     SWP_MPa = -1 / 10 * tmp_swp[, -(1:2), drop = FALSE],
#'     Snowpack_SWE_mm = 10 * tmp_snow[, "snowpackWaterEquivalent_cm"],
#'     air_Tmin_C = tmp_airtemp[, "min_C"],
#'     air_Tmean_C = tmp_airtemp[, "avg_C"],
#'     air_Tmax_C = tmp_airtemp[, "max_C"],
#'     # Using daily mean soil temperature in the absence of daily min/max
#'     shallowsoil_Tmin_C = tmp_soiltemp[, "Lyr_1"],
#'     shallowsoil_Tmean_C = tmp_soiltemp[, "Lyr_1"],
#'     shallowsoil_Tmax_C = tmp_soiltemp[, "Lyr_1"]
#'   ),
#'   soillayer_depths_cm = rSOILWAT2::swSoils_Layers(sw_in)[, 1],
#'   site_latitude = rSOILWAT2::swSite_IntrinsicSiteParams(sw_in)[["Latitude"]],
#'   years =
#'     rSOILWAT2::swYears_StartYear(sw_in):rSOILWAT2::swYears_EndYear(sw_in)
#' )
#'
#' all.equal(GISSM_r1, GISSM_r2)
#'
#' # Calculate the frequency of years when big sagebrush seeds germinate
#' # and seedlings survive in undisturbed natural vegetation
#' # (Schlaepfer et al. 2014)
#' colMeans(GISSM_r1[["outcome"]][, -1])
#'
#' @export
calc_GISSM <- function(
  x,
  soillayer_depths_cm,
  params = parameters_GISSM_bigsagebrush(),
  site_latitude = NULL,
  has_soil_temperature = NULL,
  years = NULL,
  simTime1 = NULL,
  simTime2 = NULL,
  debug_output = 0L,
  path = NULL,
  filename_tag = "GISSM"
) {

  #------ Prepare inputs
  req_vars <- c(
    "SWP_MPa", "Snowpack_SWE_mm",
    "air_Tmin_C", "air_Tmean_C", "air_Tmax_C",
    "shallowsoil_Tmin_C", "shallowsoil_Tmean_C", "shallowsoil_Tmax_C"
  )

  is_sw2 <- inherits(x, "swOutput")

  if (is_sw2) {
    sim_vals <- list()

    has_sw2_daily <- slot(x, "dy_nrow") > 0
    if (!has_sw2_daily) {
      stop(
        "This function requires that `x` has daily output if ",
        "it is a rSOILWAT2 output object."
      )
    }

    tmp <- slot(slot(x, rSW2_glovars[["swof"]]["sw_temp"]), "Day")[, 1]
    years <- unique(tmp)

  } else {
    sim_vals <- x

    has_req_vars <- all(!sapply(req_vars, function(v) is.null(sim_vals[[v]])))

    if (!has_req_vars) {
      stop(
        "This function requires either that",
        "\n\t* `x` is a rSOILWAT2 output object with daily output, or that",
        "\n\t* `x` is a list with complete forcing data, i.e., ",
        paste(shQuote(req_vars), collapse = ", ")
      )
    }

    years <- sort(unique(years))
  }

  has_bad_years <- is.null(years) || length(years) == 0


  #--- Get time sequence information
  st1_elem_names <- c(
    "useyrs", "no.usedy", "startyr", "simstartyr", "index.usedy"
  )


  is_simTime1_good <-
    !is.null(simTime1) &&
    all(!sapply(st1_elem_names, function(en) is.null(simTime1[[en]])))

  if (is_simTime1_good) {
    st1 <- simTime1

  } else {
    if (has_bad_years) {
      stop(
        "Insufficient information on time: 'years' is needed to calculate ",
        "'simTime1', but they couldn't be determined from inputs."
      )
    }

    st1 <- rSOILWAT2::setup_time_simulation_run(
      sim_time = list(
        spinup_N = 0,
        startyr = years[1],
        endyr = years[length(years)]
      )
    )
  }

  st2_elem_names <- c("year_ForEachUsedDay", "doy_ForEachUsedDay")


  is_simTime2_good <-
    !is.null(simTime2) &&
    all(!sapply(st2_elem_names, function(en) is.null(simTime2[[en]])))

  if (is_simTime2_good) {
    st2 <- simTime2

  } else {
    if (has_bad_years) {
      stop(
        "Insufficient information on time: 'years' is needed to calculate ",
        "'simTime2', but they couldn't be determined from inputs."
      )
    }

    st2 <- rSOILWAT2::simTiming_ForEachUsedTimeUnit(
      useyrs = years,
      sim_tscales = "daily",
      latitude = site_latitude,
      account_NorthSouth = TRUE
    )
  }


  #--- Extract (missing) inputs from daily rSOILWAT2 output object
  if (is_sw2) {
    if (!exists("SWP_MPa", where = sim_vals)) {
      sim_vals[["SWP_MPa"]] <- -1 / 10 * slot(
        slot(x, rSW2_glovars[["swof"]]["sw_swp"]),
        "Day"
      )[, - (1:2), drop = FALSE]
    }

    if (!exists("Snowpack_SWE_mm", where = sim_vals)) {
      sim_vals[["Snowpack_SWE_mm"]] <- 10 * slot(
        slot(x, rSW2_glovars[["swof"]]["sw_snow"]),
        "Day"
      )[, "snowpackWaterEquivalent_cm"]
    }

    if (!exists("air_Tmin_C", where = sim_vals)) {
      sim_vals[["air_Tmin_C"]] <- slot(
        slot(x, rSW2_glovars[["swof"]]["sw_temp"]),
        "Day"
      )[, "min_C"]
    }

    if (!exists("air_Tmean_C", where = sim_vals)) {
      sim_vals[["air_Tmean_C"]] <- slot(
        slot(x, rSW2_glovars[["swof"]]["sw_temp"]),
        "Day"
      )[, "avg_C"]
    }

    if (!exists("air_Tmax_C", where = sim_vals)) {
      sim_vals[["air_Tmax_C"]] <- slot(
        slot(x, rSW2_glovars[["swof"]]["sw_temp"]),
        "Day"
      )[, "max_C"]
    }

    if (!exists("shallowsoil_Tmin_C", where = sim_vals)) {
      warning("Using daily mean soil temperature instead of daily minimum.")

      sim_vals[["shallowsoil_Tmin_C"]] <- slot(
        slot(x, rSW2_glovars[["swof"]]["sw_soiltemp"]),
        "Day"
      )[, "Lyr_1"]
    }

    if (!exists("shallowsoil_Tmean_C", where = sim_vals)) {
      sim_vals[["shallowsoil_Tmean_C"]] <- slot(
        slot(x, rSW2_glovars[["swof"]]["sw_soiltemp"]),
        "Day"
      )[, "Lyr_1"]
    }

    if (!exists("shallowsoil_Tmax_C", where = sim_vals)) {
      warning("Using daily mean soil temperature instead of daily maximum.")

      sim_vals[["shallowsoil_Tmax_C"]] <- slot(
        slot(x, rSW2_glovars[["swof"]]["sw_soiltemp"]),
        "Day"
      )[, "Lyr_1"]
    }
  }


  # Check that daily forcing values are well-formed
  hasnt_req_vars <- !sapply(
    X = req_vars,
    FUN = function(v) {
      tmp <- !is.null(sim_vals[[v]])

      if (tmp) {
        switch(
          EXPR = v,
          SWP_MPa = nrow(sim_vals[[v]]) >= st1[["no.usedy"]],
          length(sim_vals[[v]]) >= st1[["no.usedy"]]
        )
      } else {
        tmp
      }
    }
  )

  if (any(hasnt_req_vars)) {
    stop(
      "Daily forcing variable(s) ",
      paste(shQuote(req_vars[hasnt_req_vars]), collapse = ", "),
      " have missing/insufficient values."
    )
  }



  #------ Prepare regeneration time

  # Regeneration year = RY:
  #   RYdoy = 1 == start of seed dispersal = start of 'regeneration year'

  # mean January (N-hemisphere) / July (S-hemisphere) air temperature
  tmp <- st1[["index.usedy"]][st2[["month_ForEachUsedDay_NSadj"]] == 1]
  mean_Jan_airTemp_C <- mean(sim_vals[["air_Tmean_C"]][tmp])

  tmp <-
    params[["Doy_SeedDispersalStart0"]] +
    params[["SeedDispersalStart_DependencyOnMeanTempJanuary"]] *
    mean_Jan_airTemp_C

  Doy_SeedDispersalStart <- as.integer(max(round(tmp, 0) %% 365, 1))

  moveByDays <- if (Doy_SeedDispersalStart > 1) {
    tmp <-
      ISOdate(st1[["useyrs"]][1] - 1, 12, 31, tz = "UTC") -
      ISOdate(st1[["useyrs"]][1] - 1, 1, 1, tz = "UTC") + 1 -
      (Doy_SeedDispersalStart - 1)
    as.integer(max(c(as.numeric(tmp) %% 365, 1)))

  } else {
    1L
  }


  # Determine dates of regeneration year
  st_RY <- list()

  et <- st1[["no.usedy"]]
  itail <- (et - moveByDays + 1):et

  if (st1[["startyr"]] > st1[["simstartyr"]]) {
    # start earlier to complete RY
    st <- st1[["index.usedy"]][1]

    # index indicating which rows of the daily SOILWAT2 output is used
    st_RY[["index.usedy"]] <- c(
      (st - moveByDays):(st - 1),
      st1[["index.usedy"]][-itail]
    )

    # 'regeneration year' for each used day
    st_RY[["year_ForEachUsedDay"]] <- st2[["year_ForEachUsedDay"]]

    # 'doy of the regeneration year' for each used day
    st_RY[["doy_ForEachUsedDay"]] <- st2[["doy_ForEachUsedDay"]]

  } else {
    # start later to get a complete RY
    fyr <- st2[["year_ForEachUsedDay"]][1]
    dy <- if (rSW2utils::isLeapYear(fyr)) 0 else 1

    tmp <- c(seq_len(Doy_SeedDispersalStart - dy), itail)
    st_RY[["index.usedy"]] <- st1[["index.usedy"]][-tmp]

    tmp <- st2[["year_ForEachUsedDay"]] == fyr
    st_RY[["year_ForEachUsedDay"]] <- st2[["year_ForEachUsedDay"]][!tmp]
    st_RY[["doy_ForEachUsedDay"]] <- st2[["doy_ForEachUsedDay"]][!tmp]
  }

  # sequence of 'regeneration years' that are used for aggregation
  st_RY[["useyrs"]] <- unique(st_RY[["year_ForEachUsedDay"]])

  # normal year for each used 'doy of the regeneration year'
  st_RY[["no.usedy"]] <- length(st_RY[["index.usedy"]])
  itail <- (st_RY[["no.usedy"]] - moveByDays + 1):st_RY[["no.usedy"]]
  st_RY[["year_ForEachUsedRYDay"]] <- c(
    rep(st1[["useyrs"]][1] - 1, moveByDays),
    st_RY[["year_ForEachUsedDay"]][-itail]
  )

  # normal doy for each used 'doy of the regeneration year'
  st <- st1[["index.usedy"]][1]
  st_RY[["doy_ForEachUsedRYDay"]] <- c(
    (st - moveByDays):(st - 1),
    st_RY[["doy_ForEachUsedDay"]][-itail]
  )


  #--- Subset daily forcing data to regeneration time
  dyf_swp <- sim_vals[["SWP_MPa"]][st_RY[["index.usedy"]], , drop = FALSE]
  dyf_snow <- sim_vals[["Snowpack_SWE_mm"]][st_RY[["index.usedy"]]]

  dyf_airTmin <- ifelse(
    dyf_snow > 0,
    params[["Temp_ExperiencedUnderneathSnowcover"]],
    sim_vals[["air_Tmin_C"]][st_RY[["index.usedy"]]]
  )

  dyf_airTmax <- sim_vals[["air_Tmax_C"]][st_RY[["index.usedy"]]]


  has_good_soil_temperature <- if (!is.null(has_soil_temperature)) {
    has_soil_temperature
  } else {
    TRUE
  }

  has_good_soil_temperature <- has_good_soil_temperature &&
    !inherits(sim_vals[["shallowsoil_Tmean_C"]], "try-error") &&
    !is.null(sim_vals[["shallowsoil_Tmean_C"]]) &&
    !anyNA(sim_vals[["shallowsoil_Tmean_C"]]) &&
    !all(sim_vals[["shallowsoil_Tmean_C"]] == 0)

  if (has_good_soil_temperature) {
    dyf_soilTmean <- ifelse(
      dyf_snow > 0,
      params[["Temp_ExperiencedUnderneathSnowcover"]],
      sim_vals[["shallowsoil_Tmean_C"]][st_RY[["index.usedy"]]]
    )

    dyf_soilTmin <- ifelse(
      dyf_snow > 0,
      params[["Temp_ExperiencedUnderneathSnowcover"]],
      sim_vals[["shallowsoil_Tmin_C"]][st_RY[["index.usedy"]]]
    )

    dyf_soilTmax <- sim_vals[["shallowsoil_Tmax_C"]][st_RY[["index.usedy"]]]

  } else {
    # air temperature is used instead of unavailable shallow soil temperature
    warning("Soil temperature is unavailable: using air temperature instead.")

    dyf_soilTmean <- ifelse(
      dyf_snow > 0,
      params[["Temp_ExperiencedUnderneathSnowcover"]],
      sim_vals[["air_Tmean_C"]][st_RY[["index.usedy"]]]
    )

    dyf_soilTmin <- dyf_airTmin
    dyf_soilTmax <- dyf_airTmax
  }


  #------ GERMINATION ------

  #--- 1. Germination periods:
  # sequence of days with favorable conditions for germination
  # defined by upper/lower limits

  # Maximal temperature for germination
  Germination_AtBelowTmax <-
    dyf_soilTmax <= params[["Temp_MaximumForGermination"]]

  # Minimal temperature for germination
  Germination_AtAboveTmin <-
    dyf_soilTmin >= params[["Temp_MinimumForGermination"]]

  # Minimum soil water for germination in relevant soil layer
  slyrs_for_germ <- SoilLayer_at_SoilDepth(
    depth_cm = params[["SoilDepth_RelevantToGermination"]],
    layers_depth = soillayer_depths_cm
  )

  if (length(slyrs_for_germ) == 1) {
    Germination_AtMoreThanTopSWPmin <-
      dyf_swp[, slyrs_for_germ] >= params[["SWP_MinimumForGermination"]]

    swp_shallow <- dyf_swp[, slyrs_for_germ]

  } else {
    Germination_AtMoreThanTopSWPmin <- apply(
      X = dyf_swp[, slyrs_for_germ],
      MARGIN = 1,
      FUN = function(x) all(x >= params[["SWP_MinimumForGermination"]])
    )

    swp_shallow <- apply(
      X = dyf_swp[, slyrs_for_germ],
      MARGIN = 1,
      FUN = mean,
      na.rm = TRUE
    )
  }

  # Put all germination limits together
  Germination_WhileFavorable <-
    Germination_AtBelowTmax &
    Germination_AtAboveTmin &
    Germination_AtMoreThanTopSWPmin


  #--- 2. Time to germinate
  # for each day with favorable conditions, determine whether
  # period of favorable conditions (resumed or reset if broken) is long enough
  # for successful completion of germination under current mean conditions

  LengthDays_FavorableConditions <- unlist(lapply(
    X = st_RY[["useyrs"]],
    FUN = calc_DurationFavorableConds,
    consequences_unfavorable = params[["GerminationPeriods_0ResetOr1Resume"]],
    Germination_WhileFavorable = Germination_WhileFavorable,
    RYyear_ForEachUsedDay = st_RY[["year_ForEachUsedDay"]]
  ))


  Germination_TimeToGerminate <- unlist(lapply(
    X = st_RY[["useyrs"]],
    FUN = calc_TimeToGerminate,
    Germination_WhileFavorable = Germination_WhileFavorable,
    LengthDays_FavorableConditions = LengthDays_FavorableConditions,
    RYyear_ForEachUsedDay = st_RY[["year_ForEachUsedDay"]],
    soilTmeanSnow = dyf_soilTmean,
    swp_shallow = swp_shallow,
    TmeanJan = mean_Jan_airTemp_C,
    params = params
  ))


  Germination_RestrictedByTimeToGerminate <- rep(FALSE, st_RY[["no.usedy"]])
  tmp <- Germination_WhileFavorable & is.na(Germination_TimeToGerminate)
  Germination_RestrictedByTimeToGerminate[tmp] <- TRUE


  #--- 3. Successful germination
  GerminationSuccess_Initiated <- !is.na(Germination_TimeToGerminate)
  germ_starts <- which(GerminationSuccess_Initiated)
  germ_durs <- Germination_TimeToGerminate[germ_starts] - 1

  if (params[["GerminationPeriods_0ResetOr1Resume"]] == 1) {
    germ_durs <-
      germ_durs +
      GISSM_germination_wait_times(
        time_to_germinate = Germination_TimeToGerminate,
        duration_fave_cond = LengthDays_FavorableConditions
      )
  }

  # index of start of successful germination + time to germinate
  # (including wait time during unfavorable conditions if 'resume')
  emergence_doys <- germ_starts + germ_durs
  Germination_Emergence <- rep(FALSE, st_RY[["no.usedy"]])
  Germination_Emergence[emergence_doys] <- TRUE
  Germination_emergence_doys <- rep(NA, st_RY[["no.usedy"]])
  Germination_emergence_doys[GerminationSuccess_Initiated] <- emergence_doys



  #------ SEEDLING SURVIVAL ------

  #--- 1. Seedling survival periods:
  #  mortality = !survival: days with conditions which kill a seedling,
  #     defined by upper/lower limits
  #  growth: days with conditions which allows a seedling to grow (here, roots),
  #     defined by upper/lower limits

  SeedlingMortality_UnderneathSnowCover <- calc_SeedlingMortality(
    kill_conds = dyf_snow > params[["SWE_MaximumForSeedlingGrowth"]],
    max_time_to_kill = params[["Days_SnowCover_MaximumForSeedlingSurvival"]]
  )

  SeedlingMortality_ByTmin <- calc_SeedlingMortality(
    kill_conds = dyf_airTmin < params[["Temp_MinimumForSeedlingSurvival"]],
    max_time_to_kill = 0
  )

  SeedlingMortality_ByTmax <- calc_SeedlingMortality(
    kill_conds = dyf_airTmax > params[["Temp_MaximumForSeedlingSurvival"]],
    max_time_to_kill = 0
  )

  SeedlingMortality_ByChronicSWPMax <- calc_SeedlingMortality(
    kill_conds = dyf_swp > params[["SWP_ChronicMaximumForSeedlingSurvival"]],
    max_time_to_kill = params[["Days_ChronicMaximumForSeedlingSurvival"]]
  )

  SeedlingMortality_ByChronicSWPMin <- calc_SeedlingMortality(
    kill_conds = dyf_swp < params[["SWP_ChronicMinimumForSeedlingSurvival"]],
    max_time_to_kill = params[["Days_ChronicMinimumForSeedlingSurvival"]]
  )

  SeedlingMortality_ByAcuteSWPMin <- calc_SeedlingMortality(
    kill_conds = dyf_swp < params[["SWP_AcuteMinimumForSeedlingSurvival"]],
    max_time_to_kill = 0
  )

  SeedlingGrowth_AbsenceOfSnowCover <-
    dyf_snow <= params[["SWE_MaximumForSeedlingGrowth"]]

  SeedlingGrowth_AtAboveTmin <-
    dyf_airTmin >= params[["Temp_MinimumForSeedlingGrowth"]]

  SeedlingGrowth_AtBelowTmax <-
    dyf_airTmax <= params[["Temp_MaximumForSeedlingGrowth"]]


  #--- 2. Grow and kill the seedlings

  # TRUE = seedling that germinate on that day and survives until end of season;
  # FALSE = no germination or seedling dies during the first seaso
  SeedlingSurvival_1stSeason <- Seedling_Starts <- Germination_Emergence

  # deep copy because Rcpp-version of get_KilledBySoilLayers changes in place
  # which would create side effects on Seedling_Starts and Germination_Emergence
  SeedlingSurvival_1stSeason[] <- SeedlingSurvival_1stSeason

  tmp <- paste0(
    "Seedlings1stSeason.Mortality.",
    c(
      "UnderneathSnowCover", "ByTmin", "ByTmax", "ByChronicSWPMax",
      "ByChronicSWPMin", "ByAcuteSWPMin",
      "DuringStoppedGrowth.DueSnowCover", "DuringStoppedGrowth.DueTmin",
      "DuringStoppedGrowth.DueTmax"
    )
  )

  SeedlingMortality_CausesByYear <- matrix(
    0,
    nrow = length(st_RY[["useyrs"]]),
    ncol = length(tmp),
    dimnames = list(NULL, tmp)
  )

  # Loop over regeneration years
  for (y in seq_along(st_RY[["useyrs"]])) {
    ids_year <- st_RY[["year_ForEachUsedDay"]] == st_RY[["useyrs"]][y]

    RYDoys_SeedlingStarts_ThisYear <- which(Seedling_Starts[ids_year])

    if (length(RYDoys_SeedlingStarts_ThisYear) > 0) {
      # if there is day with germination: init values for this year
      days_N <- sum(ids_year)
      thisYear_SeedlingMortality_UnderneathSnowCover <-
        SeedlingMortality_UnderneathSnowCover[ids_year]
      thisYear_SeedlingMortality_ByTmin <-
        SeedlingMortality_ByTmin[ids_year]
      thisYear_SeedlingMortality_ByTmax <-
        SeedlingMortality_ByTmax[ids_year]
      thisYear_SeedlingMortality_ByChronicSWPMax <-
        SeedlingMortality_ByChronicSWPMax[ids_year, , drop = FALSE]
      thisYear_SeedlingMortality_ByChronicSWPMin <-
        SeedlingMortality_ByChronicSWPMin[ids_year, , drop = FALSE]
      thisYear_SeedlingMortality_ByAcuteSWPMin <-
        SeedlingMortality_ByAcuteSWPMin[ids_year, , drop = FALSE]
      thisYear_SeedlingGrowth_AbsenceOfSnowCover <-
        SeedlingGrowth_AbsenceOfSnowCover[ids_year]
      thisYear_SeedlingGrowth_AtAboveTmin <-
        SeedlingGrowth_AtAboveTmin[ids_year]
      thisYear_SeedlingGrowth_AtBelowTmax <-
        SeedlingGrowth_AtBelowTmax[ids_year]

      # Loop over each seedling (cohort) indexed by day of germination
      for (sg_RYdoy in RYDoys_SeedlingStarts_ThisYear) {
        # init values for this seedling and season
        tmp <- seq_len(days_N)
        ids_season <- tmp[tmp > sg_RYdoy]

        # book-keeping of causes of mortality
        killed_byCauses_onRYdoy <- rep(NA, times = 6)
        names(killed_byCauses_onRYdoy) <-
          colnames(SeedlingMortality_CausesByYear)[1:6]

        # book-keeping of causes why growth stopped
        stopped_byCauses_onRYdoy <- rep(NA, times = 3)
        names(stopped_byCauses_onRYdoy) <-
          colnames(SeedlingMortality_CausesByYear)[7:9]

        # Establish days of growth (= TRUE) and surviving &no growth (= FALSE)
        thisSeedlingGrowing <- rep(TRUE, days_N)
        if (sg_RYdoy > 1) {
          # seedling germinated on sg_RYdoy,
          # hence it cannot grow before germination day
          thisSeedlingGrowing[seq_len(sg_RYdoy - 1)] <- FALSE
        }

        #--- Check growth under above-ground conditions
        # Snow cover
        thisSeedlingGrowth_AbsenceOfSnowCover <- check_SuitableGrowthThisYear(
          favorable_conditions =
            thisSeedlingGrowing & thisYear_SeedlingGrowth_AbsenceOfSnowCover,
          consequences_unfavorable = params[["SeedlingGrowth_0StopOr1Resume"]]
        )

        tmp <- !thisSeedlingGrowth_AbsenceOfSnowCover[ids_season]
        if (any(tmp)) {
          stopped_byCauses_onRYdoy["Seedlings1stSeason.Mortality.DuringStoppedGrowth.DueSnowCover"] <- sg_RYdoy + which(tmp)[1] #nolint
        }

        # Minimum temperature
        thisSeedlingGrowth_AtAboveTmin <- check_SuitableGrowthThisYear(
          favorable_conditions =
            thisSeedlingGrowing & thisYear_SeedlingGrowth_AtAboveTmin,
          consequences_unfavorable = params[["SeedlingGrowth_0StopOr1Resume"]]
        )

        tmp <- !thisSeedlingGrowth_AtAboveTmin[ids_season]
        if (any(tmp)) {
          stopped_byCauses_onRYdoy["Seedlings1stSeason.Mortality.DuringStoppedGrowth.DueTmin"] <- sg_RYdoy + which(tmp)[1] #nolint
        }

        # Maximum temperature
        thisSeedlingGrowth_AtBelowTmax <- check_SuitableGrowthThisYear(
          favorable_conditions =
            thisSeedlingGrowing & thisYear_SeedlingGrowth_AtBelowTmax,
          consequences_unfavorable = params[["SeedlingGrowth_0StopOr1Resume"]]
        )

        tmp <- !thisSeedlingGrowth_AtBelowTmax[ids_season]
        if (any(tmp)) {
          stopped_byCauses_onRYdoy["Seedlings1stSeason.Mortality.DuringStoppedGrowth.DueTmax"] <- sg_RYdoy + which(tmp)[1] #nolint
        }

        #--- Update days of growth or surviving
        thisSeedlingGrowing <-
          thisSeedlingGrowing &
          thisSeedlingGrowth_AbsenceOfSnowCover &
          thisSeedlingGrowth_AtAboveTmin &
          thisSeedlingGrowth_AtBelowTmax

        thisSeedlingLivingButNotGrowing <- !thisSeedlingGrowing

        if (sg_RYdoy > 1) {
          # seedling germinated on sg_RYdoy,
          # hence it cannot live before germination day
          thisSeedlingLivingButNotGrowing[seq_len(sg_RYdoy - 1)] <- FALSE
        }

        #--- Book-keeping survival under above-ground conditions
        tmp <- thisYear_SeedlingMortality_UnderneathSnowCover[ids_season]
        if (any(tmp)) {
          killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.UnderneathSnowCover"] <- sg_RYdoy + which(tmp)[1] - 1 #nolint
        }

        tmp <- thisYear_SeedlingMortality_ByTmin[ids_season]
        if (any(tmp)) {
          killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByTmin"] <- sg_RYdoy + which(tmp)[1] - 1 #nolint
        }

        tmp <- thisYear_SeedlingMortality_ByTmax[ids_season]
        if (any(tmp)) {
          killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByTmax"] <- sg_RYdoy + which(tmp)[1] - 1 #nolint
        }

        #--- If not killed (yet) then grow and check survival below-ground
        if (all(is.na(killed_byCauses_onRYdoy))) {
          # Grow: estimate rooting depth for this seedling for each day of year
          thisSeedling_thisYear_RootingDepth <- rep(NA, times = days_N)

          tmp <- sum(thisSeedlingGrowing)
          if (tmp > 0) {
            thisSeedlingGrowing_AgeDays <- seq_len(tmp)
            thisSeedlingGrowing_RootingDepth <- SeedlingRootingDepth(
              age = thisSeedlingGrowing_AgeDays,
              P0 = params[["Seedling_SoilDepth.PO"]],
              K = params[["Seedling_SoilDepth.K"]],
              r = params[["Seedling_SoilDepth.r"]]
            )

            thisSeedling_thisYear_RootingDepth[thisSeedlingGrowing] <-
              thisSeedlingGrowing_RootingDepth

            if (any(thisSeedlingLivingButNotGrowing, na.rm = TRUE)) {
              # for days when growth stopped then copy relevant soil depth
              stopg <- addDepths <- rle(thisSeedlingLivingButNotGrowing)
              RYDoys_stopg <- c(1, cumsum(stopg$lengths))
              for (p in seq_along(stopg$values)[stopg$values]) {
                addDepths$values[p] <- if (
                  is.na(thisSeedling_thisYear_RootingDepth[RYDoys_stopg[p]])
                ) {
                  tmp <-
                    thisSeedling_thisYear_RootingDepth[1 + RYDoys_stopg[p + 1]]

                  if (is.na(tmp)) {
                    params[["Seedling_SoilDepth.K"]]
                  } else {
                    tmp
                  }

                } else {
                  thisSeedling_thisYear_RootingDepth[RYDoys_stopg[p]]
                }
              }

              RYDoys_addDepths <- inverse.rle(addDepths)
              thisSeedling_thisYear_RootingDepth <- ifelse(
                RYDoys_addDepths > 0,
                RYDoys_addDepths,
                thisSeedling_thisYear_RootingDepth
              )
            }

          } else {
            thisSeedling_thisYear_RootingDepth[thisSeedlingLivingButNotGrowing] <- params[["Seedling_SoilDepth.PO"]] / 10 #nolint
          }

          thisSeedling_thisYear_RootingSoilLayers <- SoilLayer_at_SoilDepth(
            depth_cm = thisSeedling_thisYear_RootingDepth,
            layers_depth = soillayer_depths_cm
          )

          #--- Check survival under chronic SWPMax
          thisSeedling_thisYear_SeedlingMortality_ByChronicSWPMax <-
            GISSM_get_KilledBySoilLayers(
              relevantLayers = thisSeedling_thisYear_RootingSoilLayers,
              kill_conditions = thisYear_SeedlingMortality_ByChronicSWPMax
          )

          tmp <- thisSeedling_thisYear_SeedlingMortality_ByChronicSWPMax[ids_season] #nolint
          if (any(tmp)) {
            killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByChronicSWPMax"] <- sg_RYdoy + which(tmp)[1] - 1 #nolint
          }

          #--- Check survival under chronic SWPMin
          thisSeedling_thisYear_SeedlingMortality_ByChronicSWPMin <-
            GISSM_get_KilledBySoilLayers(
              relevantLayers = thisSeedling_thisYear_RootingSoilLayers,
              kill_conditions = thisYear_SeedlingMortality_ByChronicSWPMin
          )

          tmp <- thisSeedling_thisYear_SeedlingMortality_ByChronicSWPMin[ids_season] #nolint
          if (any(tmp)) {
            killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByChronicSWPMin"] <- sg_RYdoy + which(tmp)[1] - 1 #nolint
          }

          #--- Check survival under acute SWPMin
          thisSeedling_thisYear_SeedlingMortality_ByAcuteSWPMin <-
            GISSM_get_KilledBySoilLayers(
              relevantLayers = thisSeedling_thisYear_RootingSoilLayers,
              kill_conditions = thisYear_SeedlingMortality_ByAcuteSWPMin
          )

          tmp <- thisSeedling_thisYear_SeedlingMortality_ByAcuteSWPMin[ids_season] #nolint
          if (any(tmp)) {
            killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByAcuteSWPMin"] <- sg_RYdoy + which(tmp)[1] - 1 #nolint
          }
        }

        #--- If killed then establish which factor killed first and whether
        # and how growth was stopped before kill

        if (any(!is.na(killed_byCauses_onRYdoy))) {
          kill_factor <- which.min(killed_byCauses_onRYdoy)
          SeedlingMortality_CausesByYear[y, kill_factor] <-
            SeedlingMortality_CausesByYear[y, kill_factor] + 1
          stop_factor <- which.min(stopped_byCauses_onRYdoy)

          if (
            any(
              !is.na(stopped_byCauses_onRYdoy)) &&
              killed_byCauses_onRYdoy[kill_factor] >
                stopped_byCauses_onRYdoy[stop_factor]
          ) {
            SeedlingMortality_CausesByYear[y, 6 + stop_factor] <-
              SeedlingMortality_CausesByYear[y, 6 + stop_factor] + 1
          }

          SeedlingSurvival_1stSeason <- GISSM_kill_seedling(
            ss1s = SeedlingSurvival_1stSeason,
            ry_year_day = st_RY[["year_ForEachUsedDay"]],
            ry_useyrs = st_RY[["useyrs"]],
            y = y,
            doy = sg_RYdoy
          )
        }
      }
    } else {
      # no germination during this year -> no seedlings to grow or die
      SeedlingMortality_CausesByYear[y, ] <- NA
    }
  } # end of year loop of seedling growth


  #---Aggregate output
  GISSM <- list()

  index_RYuseyr <- unique(st_RY[["year_ForEachUsedRYDay"]]) %in% st1[["useyrs"]]

  # Number of days per year with success
  dat_gissm1 <- cbind(Germination_Emergence, SeedlingSurvival_1stSeason)
  res1_yr_v0 <- stats::aggregate(
    x = dat_gissm1,
    by = st_RY["year_ForEachUsedRYDay"],
    FUN = sum
  )
  res1_yr <- res1_yr_v0[index_RYuseyr, ]

  # Years with successful germination and seedling survival
  GISSM[["outcome"]] <- data.frame(
    ryear = res1_yr[, 1],
    res1_yr[, -1] > 0
  )

  if (isTRUE(debug_output %in% 1:2)) {
    GISSM[["successes_days"]] <- res1_yr[, -1]
    GISSM[["mortality_causes"]] <- SeedlingMortality_CausesByYear

    # Periods with no successes
    tmp <- rle(GISSM[["outcome"]][, "Germination_Emergence"])
    GISSM[["nogermination_periods_yrs"]] <- if (any(!tmp$values)) {
      tmp$lengths[!tmp$values]
    } else {
      0
    }

    tmp <- rle(GISSM[["outcome"]][, "SeedlingSurvival_1stSeason"])
    GISSM[["noseedlings_periods_yrs"]] <- if (any(!tmp$values)) {
      tmp$lengths[!tmp$values]
    } else {
      0
    }

    # Days of year (in normal count) of most frequent successes among years
    if (FALSE) {
      # convert to normal doys
      toDoy <- function(x) {
        tmp <- x + Doy_SeedDispersalStart - 1
        sort(ifelse(tmp > 365, tmp - 365, tmp))
      }
    }

    res1_dy <- stats::aggregate(
      x = dat_gissm1,
      by = st_RY["doy_ForEachUsedRYDay"],
      FUN = sum
    )

    GISSM[["success_mostfrequent_doy"]] <- get.DoyMostFrequentSuccesses(
      doys = res1_dy,
      data = dat_gissm1
    )

    # Mean number of days when germination is restricted due to conditions
    dat_gissm2 <- cbind(
      !Germination_AtBelowTmax,
      !Germination_AtAboveTmin,
      !Germination_AtMoreThanTopSWPmin,
      !Germination_WhileFavorable,
      Germination_RestrictedByTimeToGerminate
    )

    res2_yr_v0 <- stats::aggregate(
      x = dat_gissm2,
      by = st_RY["year_ForEachUsedRYDay"],
      FUN = sum
    )

    GISSM[["nogermination_days"]] <- res2_yr_v0[index_RYuseyr, -1]

    # Mean time to germinate in days
    res3_yr_v0 <- tapply(
      X = Germination_TimeToGerminate,
      INDEX = st_RY[["year_ForEachUsedRYDay"]],
      FUN = mean,
      na.rm = TRUE
    )

    GISSM[["time_to_germinate_days"]] <- res3_yr_v0[index_RYuseyr]


    # Extra output as side-effects
    if (isTRUE(debug_output == 2L)) {
      write_GISSM_debug(
        dat_gissm1,
        res1_yr_v0,
        res2_yr_v0,
        res3_yr_v0,
        SeedlingMortality_CausesByYear,
        st1,
        index_RYuseyr,
        st_RY[["year_ForEachUsedRYDay"]],
        path,
        filename_tag
      )

      plot_GISSM_debug(
        dyf_snow,
        Germination_TimeToGerminate,
        Doy_SeedDispersalStart,
        SeedlingSurvival_1stSeason,
        GerminationSuccess_Initiated,
        Germination_emergence_doys,
        Germination_RestrictedByTimeToGerminate,
        Germination_AtAboveTmin,
        Germination_AtMoreThanTopSWPmin,
        st1,
        path,
        filename_tag
      )
    }
  }

  GISSM
}



# Write internal GISSM output to spreadsheet
write_GISSM_debug <- function(dat_gissm1,
  res1_yr_v0, res2_yr_v0, res3_yr_v0,
  SeedlingMortality_CausesByYear,
  st1,
  index_RYuseyr,
  year_ForEachUsedRYDay,
  path, filename_tag
) {

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  # Table with data for every year
  res1_yr_doy <- t(simplify2array(
    by(
      dat_gissm1,
      INDICES = year_ForEachUsedRYDay,
      FUN = function(x) get.DoyMostFrequentSuccesses(x, dat_gissm1)
    )
  ))[st1[["index.useyr"]], , drop = FALSE]

  res_yr <- data.frame(
    data.frame(
      res1_yr_v0,
      res2_yr_v0[, -1],
      res3_yr_v0
    )[index_RYuseyr, ],
    SeedlingMortality_CausesByYear,
    res1_yr_doy[index_RYuseyr, ]
  )

  tmp_header2 <- c(
    "DaysWith_GerminationSuccess",
    "DaysWith_SeedlingSurvival1stSeason",
    "Days_GerminationRestrictedByTmax",
    "Days_GerminationRestrictedByTmin",
    "Days_GerminationRestrictedBySWPmin",
    "Days_GerminationRestrictedByAnyCondition",
    "Days_GerminationRestrictedByTimeToGerminate",
    "MeanDays_TimeToGerminate",
    paste(
      "Days",
      colnames(SeedlingMortality_CausesByYear),
      sep = "_"
    ),
    paste(
      rep(c("Start90%", "Median", "End90%"), times = 2),
      rep(
        c(
          "DoyMostFrequent_GerminationSuccess",
          "DoyMostFrequent_SeedlingSurvival1stSeason"
        ),
        each = 3
      ),
      sep = "_"
    )
  )

  colnames(res_yr) <- c("Year", tmp_header2)

  utils::write.csv(
    res_yr,
    file = file.path(path, paste0(filename_tag, ".csv"))
  )
}

# Visualize internal GISSM output
plot_GISSM_debug <- function(
  dyf_snow,
  Germination_TimeToGerminate,
  Doy_SeedDispersalStart,
  SeedlingSurvival_1stSeason,
  GerminationSuccess_Initiated,
  Germination_emergence_doys,
  Germination_RestrictedByTimeToGerminate,
  Germination_AtAboveTmin,
  Germination_AtMoreThanTopSWPmin,
  st1,
  path, filename_tag
) {

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  grDevices::pdf(
    file = file.path(path, paste0(filename_tag, ".pdf")),
    height = 4.5,
    width = max(4, 2 * length(st1[["index.useyr"]]))
  )

  op <- graphics::par(mar = c(1, 3, 0.1, 0.1), mgp = c(2, 0.5, 0), las = 1)
  ylim <- c(
    -17.5,
    max(
      max(dyf_snow, na.rm = TRUE),
      max(Germination_TimeToGerminate, na.rm = TRUE)
    )
  )

  p.cex <- max(0.5, min(1, exp(-0.01 * ylim[2]) + 0.5))
  xp <- seq_along(dyf_snow) + Doy_SeedDispersalStart - 1

  # Snow-water equivalents
  graphics::plot(
    xp,
    dyf_snow,
    type = "l",
    ylim = ylim,
    xlab = "Year",
    ylab = "SWE (mm), Time to germinate (days)",
    axes = FALSE
  )

  graphics::axis(
    side = 1,
    pos = ylim[1],
    at = 365 * seq_along(st1[["index.useyr"]]),
    labels = st1[["useyr"]]
  )
  graphics::axis(
    side = 2,
    pos = graphics::par("usr")[1],
    at = (tmp <- graphics::axTicks(2))[tmp >= 0]
  )

  # Time to germinate
  graphics::lines(
    xp,
    Germination_TimeToGerminate,
    col = "red",
    type = "b",
    pch = 19,
    cex = p.cex / 5
  )

  # Seedling survival
  graphics::points(
    xp,
    ifelse(SeedlingSurvival_1stSeason, 0, NA),
    col = "green",
    pch = 19
  )

  # Emergence
  tmp <- data.frame(xp, ifelse(GerminationSuccess_Initiated, -7.5, NA))
  tmp_x0 <- tmp[stats::complete.cases(tmp), ]

  tmp <- data.frame(
    Germination_emergence_doys + Doy_SeedDispersalStart - 1,
    ifelse(GerminationSuccess_Initiated, -2.5, NA)
  )
  tmp_x1 <- tmp[stats::complete.cases(tmp), ]

  graphics::segments(
    x0 = tmp_x0[, 1],
    y0 = tmp_x0[, 2],
    x1 = tmp_x1[, 1],
    y1 = tmp_x1[, 2],
    col = "blue"
  )

  # Mortality due to too short favorable conditions
  graphics::points(
    xp,
    ifelse(Germination_RestrictedByTimeToGerminate, -10, NA),
    col = "black",
    pch = 4,
    cex = p.cex
  )

  # Mortality due to too cold conditions
  graphics::points(
    xp,
    ifelse(!Germination_AtAboveTmin, -12.5, NA),
    col = grDevices::gray(0.3),
    pch = 4,
    cex = p.cex
  )

  # Mortality due to too dry conditions
  graphics::points(
    xp,
    ifelse(!Germination_AtMoreThanTopSWPmin, -15, NA),
    col = grDevices::gray(0.7),
    pch = 4,
    cex = p.cex
  )

  graphics::legend(
    "topright",
    bty = "n",
    legend = c(
      "SWE (mm)",
      "Time to germinate",
      "Seedling survival",
      "Emergence",
      "Too short favorable conditions",
      "Too cold",
      "Too dry"
    ),
    lty = c(1, 1, -1, 1, -1, -1, -1),
    pch = c(-1, -1, 19, -1, 4, 4, 4),
    col = c(
      "black", "red", "green", "blue", "black",
      grDevices::gray(0.3), grDevices::gray(0.7)
    ),
    merge = TRUE
  )

  graphics::par(op)
  grDevices::dev.off()
}

#------ End of GISSM functions ------
#------------------------------------
