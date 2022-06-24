test_that("Stan Sampling works", {
  hawkes_fit <- fit_gpmcp(
    t = as.numeric(difftime(train_df$time, lubridate::dmy_hms(010419000000, tz = "Europe/London"), units = "hours")),
    branching_structure = train_df$parent_id,
    a = 0, b = 28 * 24,
    chains = 1, refresh = 0
  )
  expect_true(TRUE)
})
