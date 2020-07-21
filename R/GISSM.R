# ----------------------------
#------ GISSM functions ------
# Schlaepfer, D.R., Lauenroth, W.K. & Bradford, J.B. (2014). Modeling
# regeneration responses of big sagebrush (Artemisia tridentata) to abiotic
# conditions. Ecol Model, 286, 66-77.

#' Function to convert soil depth to soil layer
SoilLayer_at_SoilDepth <- function(depth_cm, layers_depth) {
  pmax(1, pmin(length(layers_depth), 1 + findInterval(depth_cm - 0.01,
    layers_depth)))
}


#' Function to calculate for each day of the year, duration in days of
#' upcoming favorable conditions accounting for consequences.unfavorable = 0
#' (if conditions become unfavorable, then restart the count), =1 (resume)
calc_DurationFavorableConds <- function(RYyear, consequences.unfavorable,
  Germination_WhileFavorable, RYyear_ForEachUsedDay) {

  index.year <- RYyear_ForEachUsedDay == RYyear
  conditions <- Germination_WhileFavorable[index.year]
  doys <- seq_len(sum(index.year))
  doys[!conditions] <- NA  #calculate only for favorable days
  out <- rep(NA, times = sum(index.year))

  if (consequences.unfavorable == 0) {
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

  } else if (consequences.unfavorable == 1) {
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

#' Based on the \var{NLR} model (equation 5) in Hardegree (2006) and modified
#' by Schlaepfer et al. (2014) by making time to germinate dependent on
#' mean January temperature and soil water potential
#'
#' @references Hardegree SP (2006) Predicting Germination Response to
#'   Temperature. I. Cardinal-temperature Models and Subpopulation-specific
#'   Regression. Annals of Botany, 97, 1115-1125.
get_modifiedHardegree2006NLR <- function(RYdoy, Estimate_TimeToGerminate,
  TmeanJan, a, b, c, d, k1_meanJanTemp, k2_meanJanTempXIncubationTemp,
  k3_IncubationSWP, Tgerm.year, SWPgerm.year, durations, rec.delta = 1,
  nrec.max = 10L) {

  for (nrec in seq_len(nrec.max)) {
    Estimate_TimeToGerminate <- prev_est_TimeToGerminate <- max(0,
      round(Estimate_TimeToGerminate))

    ids <- RYdoy:(RYdoy + Estimate_TimeToGerminate - 1)
    Tgerm <- mean(Tgerm.year[ids], na.rm = TRUE)
    SWPgerm <- mean(SWPgerm.year[ids], na.rm = TRUE)

    temp.c.lim <- - (Tgerm - b) * (d ^ 2 - 1) / d
    c <- if (c > 0) {
        if (c > temp.c.lim) c else {
          temp.c.lim + SFSW2_glovars[["tol"]]
        }
      } else if (c < 0) {
        if (c < temp.c.lim) c else {
          temp.c.lim - SFSW2_glovars[["tol"]]
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
    temp <- abs(Estimate_TimeToGerminate - prev_est_TimeToGerminate)
    if (temp <= rec.delta | RYdoy + Estimate_TimeToGerminate - 1 > 365)
      break
  }

  out <- if (nrec >= nrec.max) {
      round(mean(c(Estimate_TimeToGerminate, prev_est_TimeToGerminate)), 0)
    } else {
      Estimate_TimeToGerminate
    }

  # test whether enough time to germinate
  if (out <= durations[RYdoy] & RYdoy + out <= 365) out else NA
}

#' Function to estimate time to germinate for each day of a given year and
#' conditions (temperature, top soil \var{SWP})
#'
#' @param seed A seed set, \code{NULL}, or \code{NA}. \code{NA} will not affect
#'  the state of the \var{RNG}; \code{NULL} will re-initialize the \var{RNG};
#'  and all other values are passed to \code{\link{set.seed}}.
calc_TimeToGerminate <- function(RYyear, Germination_WhileFavorable,
  LengthDays_FavorableConditions, RYyear_ForEachUsedDay, soilTmeanSnow,
  swp.TopMean, TmeanJan, param, seed = NA) {

  if (!is.na(seed)) set.seed(seed)
  runifs <- stats::runif(2)

  #values for current year
  index.year <- RYyear_ForEachUsedDay == RYyear
  conditions <- Germination_WhileFavorable[index.year]

  # determining time to germinate for every day
  a <- max(SFSW2_glovars[["tol"]], param$Hardegree_a)
  b <- param$Hardegree_b
  temp <- if (param$Hardegree_d == 1) {
      if (runifs[1] > 0.5) {
        1 + SFSW2_glovars[["tol"]]
      } else {
        1 - SFSW2_glovars[["toln"]]
      }
    } else {
      param$Hardegree_d
    }
  d <- max(SFSW2_glovars[["tol"]], temp)
  temp.c <- if (param$Hardegree_c != 0) param$Hardegree_c else {
     sign(runifs[2] - 0.5) * SFSW2_glovars[["tol"]]
    }

  # consequences of unfavorable conditions coded in here
  TimeToGerminate.favorable <- unlist(lapply(which(conditions),
    get_modifiedHardegree2006NLR,
    Estimate_TimeToGerminate = 1, TmeanJan = TmeanJan,
    a = a, b = b, c = temp.c, d = d,
    k1_meanJanTemp = param$TimeToGerminate_k1_meanJanTemp,
    k2_meanJanTempXIncubationTemp =
      param$TimeToGerminate_k2_meanJanTempXIncubationTemp,
    k3_IncubationSWP = param$TimeToGerminate_k3_IncubationSWP,
    Tgerm.year = soilTmeanSnow[index.year],
    SWPgerm.year = swp.TopMean[index.year],
    durations = LengthDays_FavorableConditions[index.year]))

  res <- rep(NA, length(conditions))
  if (length(TimeToGerminate.favorable) > 0) {
      res[conditions] <- TimeToGerminate.favorable
  }

  res
}

do.vector <- function(kill.vector, max.duration.before.kill) {
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
  mortality[kill.durations > max.duration.before.kill] <- TRUE

  mortality
}

#' Function to calculate mortality under conditions and checks survival limit
calc_SeedlingMortality <- function(kill.conditions,
  max.duration.before.kill) {

  if (length(dim(kill.conditions)) > 0) {
    # i.e., is.matrix, columns represent soil layers
    apply(kill.conditions, 2, do.vector, max.duration.before.kill)
  } else {
    do.vector(kill.conditions, max.duration.before.kill)
  }
}


#' Function to calculate favorable conditions for seedling growth for each day
#' of a given year
check_SuitableGrowthThisYear <- function(
  favorable.conditions, consequences.unfavorable) {

  out <- rep(NA, times = length(favorable.conditions))

  if (consequences.unfavorable == 0) {
    # if conditions become unfavorable, then stop growth for rest of season
    temp.rle <- rle(favorable.conditions)
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
      out <- favorable.conditions
    }

  } else if (consequences.unfavorable == 1) {
    # if conditions become unfavorable, then resume growth afterwards
    out <- favorable.conditions
  }

  out
}


#' Function to calculate rooting depth at given age
#' Units: [age] = days, [P0, K, r] = mm
#' @return A numeric vector of rooting depth in units of centimeters.
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


#' The germination and individual seedling survival model \var{GISSM}
#'
#' GISSM represents the frequency of years when big sagebrush seeds germinate
#' and seedlings survive in undisturbed natural vegetation
#' (Schlaepfer et al. 2014).
#'
#' @export
calc_GISSM <- function(
  x,
  soillayer_depths_cm,
  params = parameters_GISSM_bigsagebrush(),
  site_latitude = NULL,
  has_soil_temperature = NULL,
  simTime1 = NULL,
  simTime2 = NULL,
  plot = FALSE,
  path = NULL,
  filename_tag = "GISSM"
) {

          #---Access daily data, which do not depend on specific species parameters, i.e., start of season

          if (!exists("swpmatric.dy.all")) swpmatric.dy.all <- list(val = -1/10*slot(slot(runDataSC, swof["sw_swp"]), "Day"))  #no vwcdy available!
          temp.snow <- slot(slot(runDataSC, swof["sw_snow"]), "Day")
          temp.temp <- slot(slot(runDataSC, swof["sw_temp"]), "Day")
          TmeanJan <- mean(temp.temp[isim_time[[itime]]$index.usedy, 5][simTime2[[itime]]$month_ForEachUsedDay_NSadj == 1], na.rm = TRUE)  #mean January (N-hemisphere)/July (S-hemisphere) air temperature based on normal 'doy'
          temp.soiltemp <- slot(slot(runDataSC, swof["sw_soiltemp"]), "Day")
          if (inherits(temp.soiltemp, "try-error") || anyNA(temp.soiltemp[, -(1:2)]) || all(temp.soiltemp[, -(1:2)] == 0)) {
            use.soiltemp <- FALSE  #flag whether soil temperature output is available or not (and then air temperature is used instead of top soil temperature)
          } else {
            use.soiltemp <- TRUE  #currently we have only mean daily soil temperatures and not min/max which we need fo the model
          }

            param <- data.frame(t(opt_agg[["GISSM_params"]][, sp]))

            #Regeneration year = RY: RYdoy = 1 == start of seed dispersal = start of 'regeneration year'
            temp <- param$Doy_SeedDispersalStart0 +
              param$SeedDispersalStart_DependencyOnMeanTempJanuary * TmeanJan
            Doy_SeedDispersalStart <- as.integer(max(round(temp, 0) %% 365, 1))

            moveByDays <- if (Doy_SeedDispersalStart > 1) {
                temp <- ISOdate(isim_time[[itime]]$useyrs[1] - 1, 12, 31, tz = "UTC") -
                        ISOdate(isim_time[[itime]]$useyrs[1] - 1, 1, 1, tz = "UTC") + 1 -
                        (Doy_SeedDispersalStart - 1)
                as.integer(max(c(as.numeric(temp) %% 365, 1)))
              } else {
                1L
              }

            #Calculate regeneration year dates
            et <- isim_time[[itime]]$no.usedy
            itail <- (et - moveByDays + 1):et
            if (isim_time[[itime]][["startyr"]] > isim_time[[itime]][["simstartyr"]]) {
              #start earlier to complete RY
              st <- isim_time[[itime]]$index.usedy[1]
              RY.index.usedy <- c((st - moveByDays):(st - 1), isim_time[[itime]]$index.usedy[-itail]) #index indicating which rows of the daily SOILWAT2 output is used
              RYyear_ForEachUsedDay <- simTime2[[itime]]$year_ForEachUsedDay  #'regeneration year' for each used day
              RYdoy_ForEachUsedDay <- simTime2[[itime]]$doy_ForEachUsedDay  #'doy of the regeneration year' for each used day

            } else {
              #start later to get a complete RY
              RY.index.usedy <- isim_time[[itime]]$index.usedy[-c(1:(Doy_SeedDispersalStart - 1), itail)]
              temp <- which(simTime2[[itime]]$year_ForEachUsedDay == simTime2[[itime]]$year_ForEachUsedDay[1])
              RYyear_ForEachUsedDay <- simTime2[[itime]]$year_ForEachUsedDay[-temp]
              RYdoy_ForEachUsedDay <- simTime2[[itime]]$doy_ForEachUsedDay[-temp]
            }
            RY.useyrs <- unique(RYyear_ForEachUsedDay)  #list of 'regeneration years' that are used for aggregation

            # normal year for each used 'doy of the regeneration year'
            RY_N_usedy <- length(RY.index.usedy)
            itail <- (RY_N_usedy - moveByDays + 1):RY_N_usedy
            year_ForEachUsedRYDay <- c(rep(isim_time[[itime]]$useyrs[1] - 1, moveByDays),
                                        RYyear_ForEachUsedDay[-itail])
            # normal doy for each used 'doy of the regeneration year'
            st <- isim_time[[itime]]$index.usedy[1]
            doy_ForEachUsedRYDay <- c((st - moveByDays):(st - 1),
                                      RYdoy_ForEachUsedDay[-itail])

            #Access daily data, the first time and afterwards only if Doy_SeedDispersalStart is different from value of previous species
            if (sp == 1 || Doy_SeedDispersalStart != prev.Doy_SeedDispersalStart) {
              swp <- swpmatric.dy.all$val[RY.index.usedy, 2 + ld, drop = FALSE]
              snow <- temp.snow[RY.index.usedy, 3]*10 #mm swe in snowpack
              airTminSnow <- ifelse(snow > 0, param$Temp_ExperiencedUnderneathSnowcover, temp.temp[RY.index.usedy, 4])
              airTmax <- temp.temp[RY.index.usedy, 3]
              if (use.soiltemp) {
                soilTmeanSnow <- ifelse(snow > 0, param$Temp_ExperiencedUnderneathSnowcover, temp.soiltemp[RY.index.usedy, 3])
                soilTminSnow <- ifelse(snow > 0, param$Temp_ExperiencedUnderneathSnowcover, temp.soiltemp[RY.index.usedy, 3])
                soilTmax <- temp.soiltemp[RY.index.usedy, 3]

              } else {
                soilTmeanSnow <- ifelse(snow > 0, param$Temp_ExperiencedUnderneathSnowcover, temp.temp[RY.index.usedy, 5])
                soilTminSnow <- airTminSnow
                soilTmax <- airTmax
              }
            }

            #----GERMINATION

            #---1. Germination periods: sequence of days with favorable conditions for germination defined by upper/lower limits
            #Maximal temperature for germination
            Germination_AtBelowTmax <- soilTmax <= param$Temp_MaximumForGermination

            #Minimal temperature for germination
            Germination_AtAboveTmin <- soilTminSnow >= param$Temp_MinimumForGermination

            #Minimum soil water for germination in relevant soil layer
            SoilLayers_RelevantToGermination <- SoilLayer_at_SoilDepth(param$SoilDepth_RelevantToGermination, layers_depth)
            if (length(SoilLayers_RelevantToGermination) == 1) {
              Germination_AtMoreThanTopSWPmin <- swp[, SoilLayers_RelevantToGermination] >= param$SWP_MinimumForGermination
              swp.TopMean <- swp[, SoilLayers_RelevantToGermination]
            } else {
              Germination_AtMoreThanTopSWPmin <- apply(swp[, SoilLayers_RelevantToGermination], MARGIN = 1, FUN = function(x) all(x >= param$SWP_MinimumForGermination))
              swp.TopMean <- apply(swp[, SoilLayers_RelevantToGermination], MARGIN = 1, FUN = mean, na.rm = TRUE)
            }

            #Put all limits together
            Germination_WhileFavorable <- Germination_AtBelowTmax & Germination_AtAboveTmin & Germination_AtMoreThanTopSWPmin

            #---2. Time to germinate
            #for each day with favorable conditions, determine whether period of favorable conditions (resumed or reset if broken) is long enough for successful completion of germination under current mean conditions
            LengthDays_FavorableConditions <- unlist(lapply(RY.useyrs, FUN = calc_DurationFavorableConds,
                consequences.unfavorable = param$GerminationPeriods_0ResetOr1Resume,
                Germination_WhileFavorable = Germination_WhileFavorable,
                RYyear_ForEachUsedDay = RYyear_ForEachUsedDay))
            Germination_TimeToGerminate <- unlist(lapply(RY.useyrs, FUN = calc_TimeToGerminate,
                Germination_WhileFavorable = Germination_WhileFavorable,
                LengthDays_FavorableConditions = LengthDays_FavorableConditions,
                RYyear_ForEachUsedDay = RYyear_ForEachUsedDay,
                soilTmeanSnow = soilTmeanSnow,
                swp.TopMean = swp.TopMean,
                TmeanJan = TmeanJan, param = param))

            Germination_RestrictedByTimeToGerminate <- rep(FALSE, RY_N_usedy)
            Germination_RestrictedByTimeToGerminate[Germination_WhileFavorable & is.na(Germination_TimeToGerminate)] <- TRUE

            #---3. Successful germinations
            GerminationSuccess_Initiated <- !is.na(Germination_TimeToGerminate)
            germ.starts <- which(GerminationSuccess_Initiated)
            germ.durs <- Germination_TimeToGerminate[germ.starts] - 1
            if (param$GerminationPeriods_0ResetOr1Resume == 1) {
              germ.durs <- germ.durs + germination_wait_times(Germination_TimeToGerminate,
                LengthDays_FavorableConditions)
            }
            emergence.doys <- germ.starts + germ.durs #index of start of successful germinations + time to germinate (including wait time during unfavorable conditions if 'resume')
            Germination_Emergence <- rep(FALSE, RY_N_usedy)
            Germination_Emergence[emergence.doys] <- TRUE
            Germination_Emergence.doys <- rep(NA, RY_N_usedy)
            Germination_Emergence.doys[GerminationSuccess_Initiated] <- emergence.doys


            #----SEEDLING SURVIVAL

            #---1. Seedling survival periods:
            #  mortality = !survival: days with conditions which kill a seedling, defined by upper/lower limits
            #  growth: days with conditions which allows a seedling to grow (here, roots), defined by upper/lower limits
            SeedlingMortality_UnderneathSnowCover <- calc_SeedlingMortality(kill.conditions = (snow > param$SWE_MaximumForSeedlingGrowth), max.duration.before.kill = param$Days_SnowCover_MaximumForSeedlingSurvival)
            SeedlingMortality_ByTmin <- calc_SeedlingMortality(kill.conditions = (airTminSnow < param$Temp_MinimumForSeedlingSurvival), max.duration.before.kill = 0)
            SeedlingMortality_ByTmax <- calc_SeedlingMortality(kill.conditions = (airTmax > param$Temp_MaximumForSeedlingSurvival), max.duration.before.kill = 0)
            SeedlingMortality_ByChronicSWPMax <- calc_SeedlingMortality(kill.conditions = (swp > param$SWP_ChronicMaximumForSeedlingSurvival), max.duration.before.kill = param$Days_ChronicMaximumForSeedlingSurvival)
            SeedlingMortality_ByChronicSWPMin <- calc_SeedlingMortality(kill.conditions = (swp < param$SWP_ChronicMinimumForSeedlingSurvival), max.duration.before.kill = param$Days_ChronicMinimumForSeedlingSurvival)
            SeedlingMortality_ByAcuteSWPMin <- calc_SeedlingMortality(kill.conditions = (swp < param$SWP_AcuteMinimumForSeedlingSurvival), max.duration.before.kill = 0)

            SeedlingGrowth_AbsenceOfSnowCover <- (snow <= param$SWE_MaximumForSeedlingGrowth)
            SeedlingGrowth_AtAboveTmin <- (airTminSnow >= param$Temp_MinimumForSeedlingGrowth)
            SeedlingGrowth_AtBelowTmax <- (airTmax <= param$Temp_MaximumForSeedlingGrowth)

            #---2. Grow and kill the seedlings
            SeedlingSurvival_1stSeason <- Seedling_Starts <- Germination_Emergence #TRUE = seedling that germinated on that day and survives until end of season; FALSE = no germination or seedling dies during the first season
            SeedlingSurvival_1stSeason[] <- SeedlingSurvival_1stSeason # deep copy because Rcpp-version of get_KilledBySoilLayers changes in place which has otherwise side effects on Seedling_Starts and Germination_Emergence
            SeedlingMortality_CausesByYear <- matrix(0, nrow = length(RY.useyrs), ncol = 9)
            colnames(SeedlingMortality_CausesByYear) <- paste0("Seedlings1stSeason.Mortality.", c("UnderneathSnowCover", "ByTmin", "ByTmax", "ByChronicSWPMax", "ByChronicSWPMin", "ByAcuteSWPMin",
                    "DuringStoppedGrowth.DueSnowCover", "DuringStoppedGrowth.DueTmin", "DuringStoppedGrowth.DueTmax"))
            for (y in seq_along(RY.useyrs)) {#for each year
              index.thisYear <- RYyear_ForEachUsedDay == RY.useyrs[y]
              RYDoys_SeedlingStarts_ThisYear <- which(Seedling_Starts[index.thisYear])
              if (length(RYDoys_SeedlingStarts_ThisYear) > 0) {#if there are any germinations
                #init values for this year
                no.days <- sum(index.thisYear)
                thisYear_SeedlingMortality_UnderneathSnowCover <- SeedlingMortality_UnderneathSnowCover[index.thisYear]
                thisYear_SeedlingMortality_ByTmin <- SeedlingMortality_ByTmin[index.thisYear]
                thisYear_SeedlingMortality_ByTmax <- SeedlingMortality_ByTmax[index.thisYear]
                thisYear_SeedlingMortality_ByChronicSWPMax <- SeedlingMortality_ByChronicSWPMax[index.thisYear, , drop = FALSE]
                thisYear_SeedlingMortality_ByChronicSWPMin <- SeedlingMortality_ByChronicSWPMin[index.thisYear, , drop = FALSE]
                thisYear_SeedlingMortality_ByAcuteSWPMin <- SeedlingMortality_ByAcuteSWPMin[index.thisYear, , drop = FALSE]
                thisYear_SeedlingGrowth_AbsenceOfSnowCover <- SeedlingGrowth_AbsenceOfSnowCover[index.thisYear]
                thisYear_SeedlingGrowth_AtAboveTmin <- SeedlingGrowth_AtAboveTmin[index.thisYear]
                thisYear_SeedlingGrowth_AtBelowTmax <- SeedlingGrowth_AtBelowTmax[index.thisYear]

                for (sg_RYdoy in RYDoys_SeedlingStarts_ThisYear) {#for each seedling indexed by day of germination
                  #init values for this seedling and season
                  temp <- seq_len(no.days)
                  index.thisSeedlingSeason <- temp[temp > sg_RYdoy]
                  killed_byCauses_onRYdoy <- rep(NA, times = 6)  #book-keeping of mortality causes
                  names(killed_byCauses_onRYdoy) <- colnames(SeedlingMortality_CausesByYear)[1:6]
                  stopped_byCauses_onRYdoy <- rep(NA, times = 3)  #book-keeping of causes why growth stopped
                  names(stopped_byCauses_onRYdoy) <- colnames(SeedlingMortality_CausesByYear)[7:9]

                  #Establish days of growth ( = TRUE) and surviving, but no growth ( = FALSE)
                  thisSeedlingGrowing <- rep(TRUE, no.days)
                  if (sg_RYdoy > 1)
                    thisSeedlingGrowing[seq_len(sg_RYdoy - 1)] <- FALSE  #seedling germinated on sg_RYdoy, hence it cannot grow before germination day

                  #Check growth under above-ground conditions
                  #Snow cover
                  thisSeedlingGrowth_AbsenceOfSnowCover <- check_SuitableGrowthThisYear(favorable.conditions = thisSeedlingGrowing & thisYear_SeedlingGrowth_AbsenceOfSnowCover, consequences.unfavorable = param$SeedlingGrowth_0StopOr1Resume)
                  temp <- !thisSeedlingGrowth_AbsenceOfSnowCover[index.thisSeedlingSeason]
                  if (any(temp))
                    stopped_byCauses_onRYdoy["Seedlings1stSeason.Mortality.DuringStoppedGrowth.DueSnowCover"] <- sg_RYdoy + which(temp)[1]
                  #Minimum temperature
                  thisSeedlingGrowth_AtAboveTmin <- check_SuitableGrowthThisYear(favorable.conditions = thisSeedlingGrowing & thisYear_SeedlingGrowth_AtAboveTmin, consequences.unfavorable = param$SeedlingGrowth_0StopOr1Resume)
                  temp <- !thisSeedlingGrowth_AtAboveTmin[index.thisSeedlingSeason]
                  if (any(temp))
                    stopped_byCauses_onRYdoy["Seedlings1stSeason.Mortality.DuringStoppedGrowth.DueTmin"] <- sg_RYdoy + which(temp)[1]
                  #Maximum temperature
                  thisSeedlingGrowth_AtBelowTmax <- check_SuitableGrowthThisYear(favorable.conditions = thisSeedlingGrowing & thisYear_SeedlingGrowth_AtBelowTmax, consequences.unfavorable = param$SeedlingGrowth_0StopOr1Resume)
                  temp <- !thisSeedlingGrowth_AtBelowTmax[index.thisSeedlingSeason]
                  if (any(temp))
                    stopped_byCauses_onRYdoy["Seedlings1stSeason.Mortality.DuringStoppedGrowth.DueTmax"] <- sg_RYdoy + which(temp)[1]
                  #Updated days of growth or surviving
                  thisSeedlingGrowing <- thisSeedlingGrowing & thisSeedlingGrowth_AbsenceOfSnowCover & thisSeedlingGrowth_AtAboveTmin & thisSeedlingGrowth_AtBelowTmax
                  thisSeedlingLivingButNotGrowing <- !thisSeedlingGrowing
                  if (sg_RYdoy > 1)
                    thisSeedlingLivingButNotGrowing[seq_len(sg_RYdoy - 1)] <- FALSE  #seedling germinated on sg_RYdoy, hence it cannot live before germination day

                  #Book-keeping survival under above-ground conditions
                  temp <- thisYear_SeedlingMortality_UnderneathSnowCover[index.thisSeedlingSeason]
                  if (any(temp))
                    killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.UnderneathSnowCover"] <- sg_RYdoy + which(temp)[1] - 1
                  temp <- thisYear_SeedlingMortality_ByTmin[index.thisSeedlingSeason]
                  if (any(temp))
                    killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByTmin"] <- sg_RYdoy + which(temp)[1] - 1
                  temp <- thisYear_SeedlingMortality_ByTmax[index.thisSeedlingSeason]
                  if (any(temp))
                    killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByTmax"] <- sg_RYdoy + which(temp)[1] - 1

                  #If not killed (yet) then grow and check survival below-ground
                  if (all(is.na(killed_byCauses_onRYdoy))) {
                    #Grow: estimate rooting depth for this seedling for each day of this year
                    thisSeedling_thisYear_RootingDepth <- rep(NA, times = no.days)
                    temp <- sum(thisSeedlingGrowing)
                    if (temp > 0) {
                      thisSeedlingGrowing_AgeDays <- seq_len(temp)
                      thisSeedlingGrowing_RootingDepth <- SeedlingRootingDepth(thisSeedlingGrowing_AgeDays, param$Seedling_SoilDepth.PO, param$Seedling_SoilDepth.K, param$Seedling_SoilDepth.r)
                      thisSeedling_thisYear_RootingDepth[thisSeedlingGrowing] <- thisSeedlingGrowing_RootingDepth
                      if (any(thisSeedlingLivingButNotGrowing, na.rm = TRUE)) {
                        #for days when growth stopped then copy relevant soil depth
                        stopg <- addDepths <- rle(thisSeedlingLivingButNotGrowing)
                        RYDoys_stopg <- c(1, cumsum(stopg$lengths))
                        for (p in seq_along(stopg$values)[stopg$values]) {
                          addDepths$values[p] <- if (is.na(thisSeedling_thisYear_RootingDepth[RYDoys_stopg[p]])) {
                              if (is.na(thisSeedling_thisYear_RootingDepth[1 + RYDoys_stopg[p+1]])) {
                                  param$Seedling_SoilDepth.K
                                } else {
                                  thisSeedling_thisYear_RootingDepth[1 + RYDoys_stopg[p+1]]
                                }
                            } else {
                              thisSeedling_thisYear_RootingDepth[RYDoys_stopg[p]]
                            }
                        }
                        RYDoys_addDepths <- inverse.rle(addDepths)
                        thisSeedling_thisYear_RootingDepth <- ifelse(RYDoys_addDepths > 0, RYDoys_addDepths, thisSeedling_thisYear_RootingDepth)
                      }
                    } else {
                      thisSeedling_thisYear_RootingDepth[thisSeedlingLivingButNotGrowing] <- param$Seedling_SoilDepth.PO/10
                    }
                    thisSeedling_thisYear_RootingSoilLayers <- SoilLayer_at_SoilDepth(thisSeedling_thisYear_RootingDepth, layers_depth)

                    #Check survival under chronic SWPMax
                    thisSeedling_thisYear_SeedlingMortality_ByChronicSWPMax <- get_KilledBySoilLayers(thisSeedling_thisYear_RootingSoilLayers, thisYear_SeedlingMortality_ByChronicSWPMax)
                    temp <- thisSeedling_thisYear_SeedlingMortality_ByChronicSWPMax[index.thisSeedlingSeason]
                    if (any(temp))
                      killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByChronicSWPMax"] <- sg_RYdoy + which(temp)[1] - 1
                    #Check survival under chronic SWPMin
                    thisSeedling_thisYear_SeedlingMortality_ByChronicSWPMin <- get_KilledBySoilLayers(thisSeedling_thisYear_RootingSoilLayers, thisYear_SeedlingMortality_ByChronicSWPMin)
                    temp <- thisSeedling_thisYear_SeedlingMortality_ByChronicSWPMin[index.thisSeedlingSeason]
                    if (any(temp))
                      killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByChronicSWPMin"] <- sg_RYdoy + which(temp)[1] - 1
                    #Check survival under acute SWPMin
                    thisSeedling_thisYear_SeedlingMortality_ByAcuteSWPMin <- get_KilledBySoilLayers(thisSeedling_thisYear_RootingSoilLayers, thisYear_SeedlingMortality_ByAcuteSWPMin)
                    temp <- thisSeedling_thisYear_SeedlingMortality_ByAcuteSWPMin[index.thisSeedlingSeason]
                    if (any(temp))
                      killed_byCauses_onRYdoy["Seedlings1stSeason.Mortality.ByAcuteSWPMin"] <- sg_RYdoy + which(temp)[1] - 1
                  }
                  #If killed then establish which factor killed first and if and how growth was stopped before kill
                  if (any(!is.na(killed_byCauses_onRYdoy))) {
                    kill.factor <- which.min(killed_byCauses_onRYdoy)
                    SeedlingMortality_CausesByYear[y, kill.factor] <- SeedlingMortality_CausesByYear[y, kill.factor] + 1
                    stop.factor <- which.min(stopped_byCauses_onRYdoy)
                    if (any(!is.na(stopped_byCauses_onRYdoy)) &&
                        killed_byCauses_onRYdoy[kill.factor] > stopped_byCauses_onRYdoy[stop.factor]) {
                      SeedlingMortality_CausesByYear[y, 6+stop.factor] <- SeedlingMortality_CausesByYear[y, 6+stop.factor] + 1
                    }
                    SeedlingSurvival_1stSeason <- kill_seedling(
                      SeedlingSurvival_1stSeason, RYyear_ForEachUsedDay,
                      RY.useyrs, y, sg_RYdoy)
                  }
                }
              } else {#no germination during this year -> no seedlings to grow or die
                SeedlingMortality_CausesByYear[y, ] <- NA
              }
            }#end of year loop of seedling growth

            #---Aggregate output
            dat_gissm1 <- cbind(Germination_Emergence, SeedlingSurvival_1stSeason)
            dat_gissm2 <- cbind(!Germination_AtBelowTmax, !Germination_AtAboveTmin,
              !Germination_AtMoreThanTopSWPmin, !Germination_WhileFavorable,
              Germination_RestrictedByTimeToGerminate)

            # Number of days per year with success
            index_RYuseyr <- unique(year_ForEachUsedRYDay) %in% isim_time[[itime]]$useyr
            res1.yr_v0 <- stats::aggregate(dat_gissm1, by = list(year_ForEachUsedRYDay), FUN = sum)
            res1.yr <- res1.yr_v0[index_RYuseyr, -1]

            # Frequency of years with success
            stemp <- res1.yr > 0

            # Periods with no successes
            rleGerm <- rle(stemp[, 1])
            rleSling <- rle(stemp[, 2])

            # Days of year (in normal count) of most frequent successes among years
            # toDoy <- function(x) sort(ifelse((temp <- x+Doy_SeedDispersalStart-1) > 365, temp-365, temp)) #convert to normal doys
            res1.dy <- stats::aggregate(dat_gissm1, by = list(doy_ForEachUsedRYDay), FUN = sum)
            mfs_doy <- get.DoyMostFrequentSuccesses(res1.dy, dat_gissm1)

            # Mean number of days when germination is restricted due to conditions
            res2.yr_v0 <- stats::aggregate(dat_gissm2, by = list(year_ForEachUsedRYDay), sum)
            res2.yr <- res2.yr_v0[index_RYuseyr, -1]

            # Mean time to germinate in days
            res3.yr_v0 <- tapply(Germination_TimeToGerminate, year_ForEachUsedRYDay, mean, na.rm = TRUE)
            res3.yr <- res3.yr_v0[index_RYuseyr]

            # Extra output as sideeffects
            if (plot) {
              dir.create(path, recursive = TRUE, showWarnings = FALSE)

              #Table with data for every year
              res1.yr.doy <- t(simplify2array(by(dat_gissm1, INDICES = year_ForEachUsedRYDay,
                FUN = function(x) get.DoyMostFrequentSuccesses(x, dat_gissm1))))[isim_time[[itime]]$index.useyr, , drop = FALSE]

              res.yr <- data.frame(data.frame(res1.yr_v0, res2.yr_v0[, -1], res3.yr_v0)[index_RYuseyr, ], SeedlingMortality_CausesByYear, res1.yr.doy)
              temp.header2 <- c("DaysWith_GerminationSuccess", "DaysWith_SeedlingSurvival1stSeason",
                "Days_GerminationRestrictedByTmax", "Days_GerminationRestrictedByTmin",
                "Days_GerminationRestrictedBySWPmin", "Days_GerminationRestrictedByAnyCondition",
                "Days_GerminationRestrictedByTimeToGerminate", "MeanDays_TimeToGerminate",
                paste("Days", colnames(SeedlingMortality_CausesByYear), sep = "_"),
                paste(rep(c("Start90%", "Median", "End90%"), times = 2),
                  rep(c("DoyMostFrequent_GerminationSuccess", "DoyMostFrequent_SeedlingSurvival1stSeason"),
                    each = 3), sep = "_"))
              colnames(res.yr) <- c("Year", temp.header2)
              utils::write.csv(
                res.yr,
                file = file.path(path, paste0(filename_tag, ".csv"))
              )

              #Plot with data for every day
              grDevices::pdf(
                file = file.path(path, paste0(filename_tag, ".pdf")),
                width = max(4, 2*length(isim_time[[itime]]$index.useyr)),
                height = 4.5
              )

              op <- graphics::par(mar = c(1, 3, 0.1, 0.1), mgp = c(2, 0.5, 0), las = 1)
              ylim <- c(-17.5, max(max(snow, na.rm = TRUE), max(Germination_TimeToGerminate, na.rm = TRUE)))
              p.cex <- max(0.5, min(1, exp(-0.01 * ylim[2]) + 0.5))
              xp <- 1:length(snow) + Doy_SeedDispersalStart-1

              graphics::plot(xp, snow, type = "l", ylim = ylim, xlab = "Year", ylab = "SWE (mm), Time to germinate (days)", axes = FALSE)
              graphics::axis(1, pos = ylim[1], at = 365*(1:(length(isim_time[[itime]]$index.useyr))), labels = isim_time[[itime]]$useyr)
              graphics::axis(2, pos = graphics::par("usr")[1], at = (temp <- graphics::axTicks(2))[temp >= 0])
              graphics::lines(xp, Germination_TimeToGerminate, col = "red", type = "b", pch = 19, cex = p.cex/5)
              graphics::points(xp, ifelse(SeedlingSurvival_1stSeason, 0, NA), col = "green", pch = 19)
              x0.temp <- (temp <- data.frame(xp, ifelse(GerminationSuccess_Initiated, -7.5, NA)))[stats::complete.cases(temp), ]
              x1.temp <- (temp <- data.frame(Germination_Emergence.doys + Doy_SeedDispersalStart-1, ifelse(GerminationSuccess_Initiated, -2.5, NA)))[stats::complete.cases(temp), ]
              graphics::segments(x0 = x0.temp[, 1], y0 = x0.temp[, 2], x1 = x1.temp[, 1], y1 = x1.temp[, 2], col = "blue")
              graphics::points(xp, ifelse(Germination_RestrictedByTimeToGerminate, -10, NA), col = "black", pch = 4, cex = p.cex)
              graphics::points(xp, ifelse(!Germination_AtAboveTmin, -12.5, NA), col = grDevices::gray(0.3), pch = 4, cex = p.cex)
              graphics::points(xp, ifelse(!Germination_AtMoreThanTopSWPmin, -15, NA), col = grDevices::gray(0.7), pch = 4, cex = p.cex)
              graphics::mtext(i_label)
              graphics::legend("topright", legend = c("SWE", "Time to germinate", "Seedling survival", "Emergence", "Too short favorable conditions", "Too cold", "Too dry"),
                bty = "n", lty = c(1, 1, -1, 1, -1, -1, -1), pch = c(-1, -1, 19, -1, 4, 4, 4), col = c("black", "red", "green", "blue", "black", grDevices::gray(0.3), grDevices::gray(0.7)), merge = TRUE)
              graphics::par(op)
              grDevices::dev.off()
            }

  list(
    freq = stemp,
    res1.yr = res1.yr,
    rleGerm = rleGerm,
    rleSling = rleSling,
    mfs_doy = mfs_doy,
    res2.yr = res2.yr,
    res3.yr = res3.yr,
    SeedlingMortality_CausesByYear = SeedlingMortality_CausesByYear
  )
}

#------ End of GISSM functions ------
#------------------------------------
