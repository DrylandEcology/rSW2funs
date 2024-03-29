# rSW2funs v0.2.0
* `rSW2funs` now uses new soil water retention curve `SWRC` functionality
  if `rSOILWAT2` v6.0.0 or later is available and uses the old interface
  if `rSOILWAT2` earlier than v6.0.0; affected functions: `calc_SMTRs()`
  (#10; @dschlaep).

# rSW2funs v0.1.3
* Linting updated to `lintr` >= 3 and
  lint workflow switched from package tests to Github Action (#5).
* Github Actions are triggered for `release/**` branches in addition to `main`.
* `calc_GISSM()` and `calc_SMTRs()` now handle `rSOILWAT2` `v3.5.0` output
  with new minimum, average, and maximum soil temperature at surface and
  at layer depth (updated column names) and continue to support output from
  earlier `rSOILWAT2` versions.


# rSW2funs v0.1.2
* `calc_GISSM()` now produces correct output even for a non-leap start year
  (@CaitlinA, #2).
* `calc_GISSM()` now checks more cases of possible inconsistencies
  in time information extracted from
  arguments `x`, `years`, `simTime1`, and `simTime2`.


# rSW2funs v0.1.1
* SOILWAT2-related functionality was re-organized, i.e.,
  a family of rSW2-related packages was created
  (see https://github.com/DrylandEcology/rSOILWAT2#rSW2)
* New dependency on new packages rSW2utils and rSW2data


# rSW2funs v0.1.0
* Initial release

