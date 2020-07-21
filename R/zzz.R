###############################################################################
# rSW2funs
#    Copyright (C) {2020}  {Daniel Schlaepfer}
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################


.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage(
      "Package ", shQuote(pkgname),
      " v", utils::packageVersion(pkgname),
      " (", utils::packageDate(pkgname), ")",
      " attached/loaded."
    )
  }

  invisible()
}


.onLoad <- function(libname, pkgname) {
  #--- Define package level variables that should be hidden from package user
  # 'rSW2_glovars' is defined in "rSW2funs-package.R"

  assign("swof", rSOILWAT2::sw_out_flags(), envir = rSW2_glovars)
  assign("tol", sqrt(.Machine$double.eps), envir = rSW2_glovars)
  assign("toln", sqrt(.Machine$double.neg.eps), envir = rSW2_glovars)

  invisible()
}


.onUnload <- function(libpath) {
  #--- Clean up C code
  library.dynam.unload("rSW2funs", libpath)

  invisible()
}
