#' \pkg{rSW2funs}: Collection of functions that calculate new variables
#' using rSOILWAT2 simulation output.
#'
#' @section LICENSE:
#'    Copyright (C) \Sexpr{format(Sys.Date(), "\%Y")} by
#'    \Sexpr{packageDescription("rSW2funs")[["Maintainer"]]}
#'
#'    This program is free software: you can redistribute it and/or modify
#'    it under the terms of the GNU General Public License as published by
#'    the Free Software Foundation, version 3 of the License.
#'
#' @section DISCLAIMER:
#'    This program is distributed in the hope that it will be useful,
#'    but WITHOUT ANY WARRANTY; without even the implied warranty of
#'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#'    GNU General Public License for more details.
#'
#' @useDynLib rSW2funs, .registration = TRUE
#' @docType package
#' @name rSW2funs
"_PACKAGE"


##------ Package level variables
rSW2_glovars <- new.env()


##------ Import from other packages
#' @import rSW2utils
#' @import rSOILWAT2
## Package uses S3/S4 classes - they are defined in package:methods
#' @import methods
NULL


##------ Support Rcpp
#' @importFrom Rcpp sourceCpp evalCpp
NULL
