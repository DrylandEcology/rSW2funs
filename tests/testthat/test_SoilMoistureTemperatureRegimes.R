context("Calculate soil moisture/texture regimes")


test_that("SMTR", {
  #--- Expectations
  expected_STR_name <- "Cryic"
  tmp <- STR_names()
  expected_STR <- rep(0, length(tmp))
  names(expected_STR) <- tmp
  expected_STR[expected_STR_name] <- 1

  expected_SMR_name <- c("Ustic", "Typic-Tempustic")
  tmp <- c(SMR_names(), SMRq_names())
  expected_SMR <- rep(0, length(tmp))
  names(expected_SMR) <- tmp
  expected_SMR[expected_SMR_name] <- 1


  #--- INPUTS
  sw_in <- rSOILWAT2::sw_exampleData

  #--- Run SOILWAT2
  sw_out <- rSOILWAT2::sw_exec(inputData = sw_in)

  #--- Calculate soil moisture/texture regimes
  SMTR <- calc_SMTRs(sim_in = sw_in, sim_out = sw_out)

  expect_true(SMTR[["regimes_done"]])

  expect_true(all(colnames(SMTR[["STR"]]) %in% STR_names()))

  expect_true(
    all(colnames(SMTR[["SMR"]]) %in% c(SMR_names(), SMRq_names()))
  )

  expect_equivalent(SMTR[["STR"]], expected_STR)

  expect_equivalent(SMTR[["SMR"]], expected_SMR)
})
