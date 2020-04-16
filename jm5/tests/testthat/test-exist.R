library(testthat)

expect_error(
  jm_filter(CD4 ~ age, Surv(Time, death) ~ drug, JM::aids, JM::aids$patient),
  "does not exist")
