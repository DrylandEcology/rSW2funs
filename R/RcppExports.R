# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Determine wait times until germination based on information on favorable
#'   conditions and time required to germinate
#'
#' @param time_to_germinate Numerical vector. Elements represent calendar days
#'  and values represent number of days required to germination
#'  if started on that day.
#' @param duration_fave_cond Numerical vector. Elements represent calendar days
#'  and values represent number of days with favorable conditions
#'  for germination.
#'
#'
#' @section Note: The \pkg{Rcpp} version of the function is about 270x faster
#'  for vectors of length 365 and 12,000x faster for vectors of length 11,000
#'  than the R version. The \pkg{Rcpp} version also reduced the memory
#'  footprint by a factor of >> 3080.
#'
#' @references Schlaepfer, D.R., Lauenroth, W.K. & Bradford, J.B. (2014).
#'  Modeling regeneration responses of big sagebrush (Artemisia tridentata)
#'  to abiotic conditions. Ecol Model, 286, 66-77.
#'
#' @examples
#'  # The \pkg{Rcpp} function is equivalent to the following R version
#'    germination_wait_times_R <- function(time_to_germinate, duration_fave_cond) {
#'      N <- length(time_to_germinate)
#'      stats::na.exclude(unlist(lapply(seq_len(N), function(t) {
#'        if (is.finite(time_to_germinate[t])) {
#'          t1 <- duration_fave_cond[t:N]
#'          t2 <- stats::na.exclude(t1)
#'          t3 <- which(t2[time_to_germinate[t]] == t1)[1]
#'          sum(is.na(t1[1:t3]))
#'        } else {
#'          NA
#'        }
#'      })))
#'    }
#'
#' @export
GISSM_germination_wait_times <- function(time_to_germinate, duration_fave_cond) {
    .Call(`_rSW2funs_GISSM_germination_wait_times`, time_to_germinate, duration_fave_cond)
}

#' Determine if all conditions across rooted soil layers are deadly
#'
#' Function that checks whether all relevant (those with roots) soil layers
#'  are under conditions of mortality (kill.conditions) for each day of a
#'  given year
#'
#'  \code{relevantLayers} takes either \code{NA} if no soil layers should be
#'  considered (e.g., because not yet germinated), or an integer number
#'  between 1 and the number of simulated soil layers. The number indicates
#'  the depth to which a seedling has grown roots and over which layers
#'  \code{kill.conditions} will be evaluated.
#'
#' @section Note: The \pkg{Rcpp} version of the function is about 165x
#'  faster than the version previous to commit
#'  \var{6344857a9cdb08acf68fa031c43cf4a596613aad} 'Small speed improvements'
#'  and about 70x faster than the R version. The \pkg{Rcpp} version also
#'  reduced the memory footprint by a factor of 200.
#'
#' @param relevantLayers An integer vector, usually of length 365 or 366
#'  (days).
#' @param kill_conditions A m x p logical matrix with
#'  \code{m >= length(relevantLayers)} and p represents the number of
#'  simulated soil layers, i.e., \code{p >= max(relevantLayers, na.rm = TRUE)}.
#'
#' @references Schlaepfer, D.R., Lauenroth, W.K. & Bradford, J.B. (2014).
#'  Modeling regeneration responses of big sagebrush (Artemisia tridentata)
#'  to abiotic conditions. Ecol Model, 286, 66-77.
#'
#' @return A logical vector of the length of \code{relevantLayers} with
#'  values containing \code{NA} for days when conditions were not evaluated,
#'  \code{TRUE} if all relevant soil layers (columns) of \code{kill.conditions}
#'  were \code{TRUE}, and with \code{FALSE} otherwise
#'
#' @examples
#'  # The \pkg{Rcpp} function is equivalent to the following R version
#'     get_KilledBySoilLayers_R <- function(relevantLayers, kill.conditions) {
#'       vapply(
#'         seq_along(relevantLayers),
#'         function(k) {
#'           if (all(is.finite(relevantLayers[k]))) {
#'             all(as.logical(kill.conditions[k, seq_len(relevantLayers[k])]))
#'           } else NA
#'         },
#'         FUN.VALUE = NA
#'       )
#'     }
#'
#' @export
GISSM_get_KilledBySoilLayers <- function(relevantLayers, kill_conditions) {
    .Call(`_rSW2funs_GISSM_get_KilledBySoilLayers`, relevantLayers, kill_conditions)
}

#' Set seedling as dead for a given day in a given year (\var{\sQuote{ss1s}})
#'
#' @param ss1s Logical vector. Elements represent calendar days included in
#'   all calendar year of \code{ry_useyrs}.
#' @param ry_year_day Numerical vector. Elements represent calendar days
#'   and values represent the calendar year of each day.
#' @param ry_useyrs Numerical vector. The sequence of calendar year.
#' @param y Numerical value. The index of the currently examined calendar year.
#' @param doy Numerical value. The index of the currently examined calendar day.
#'
#' @section Note: The \pkg{Rcpp} version of the function is about 270x faster
#'  for vectors of length 365 and 12,000x faster for vectors of length 11,000
#'  than the R version. The \pkg{Rcpp} version also reduced the memory
#'  footprint by a factor of >> 3080.
#' @section Note: Previous name \code{setFALSE_SeedlingSurvival_1stSeason}.
#'
#' @section C code: \code{ss1s} is a pointer to the data and the original
#'  vector will get altered; one would need for a deep copy:
#'  \code{LogicalVector out = clone(ss1s)}
#'
#' @references Schlaepfer, D.R., Lauenroth, W.K. & Bradford, J.B. (2014).
#'  Modeling regeneration responses of big sagebrush (Artemisia tridentata)
#'  to abiotic conditions. Ecol Model, 286, 66-77.
#'
#' @examples
#'  # The \pkg{Rcpp} function is equivalent to the following R version
#'    kill_seedling_R <- function(ss1s, ry_year_day, ry_useyrs, y,
#'      doy) {
#'      ss1s[ry_year_day == ry_useyrs[y]][doy] <- FALSE
#'      ss1s
#'    }
#'
#' @export
GISSM_kill_seedling <- function(ss1s, ry_year_day, ry_useyrs, y, doy) {
    .Call(`_rSW2funs_GISSM_kill_seedling`, ss1s, ry_year_day, ry_useyrs, y, doy)
}

