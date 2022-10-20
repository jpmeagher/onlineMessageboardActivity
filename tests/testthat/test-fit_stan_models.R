library(dplyr)
library(lubridate)

test_that("Stan Sampling works", {
  df <- train_df %>%
    mutate(t = as.numeric(difftime(time, lubridate::dmy_hms(010419000000, tz = "Europe/London"), units = "hours"))) %>%
    filter(as.numeric(tree) < 500) %>%
    arrange(tree) %>%
    group_by(tree) %>%
    mutate(tau = t - min(t)) %>%
    filter(tau < 12) %>%
    mutate(beta = refactor_branching_structure(
      id = id, parent_id = parent_id, is_immigrant = parent_id == 0, S = length(id)
    )) %>%
    ungroup()

  unit <- 24
  omega <- structure(2 * pi * (1:2 / unit), dim = 2)

  # fit_M0 <- fit_branching_point_process(
  #   t = df$t, branching_structure = df$beta, observation_interval = 12,
  #   chains = 1, refresh = 0
  # )
  # fit_M1 <- fit_branching_point_process(
  #   t = df$t, branching_structure = df$beta, observation_interval = 12,
  #   point_type = 1 + (df$beta > 0),
  #   chains = 1, refresh = 0
  # )
  # fit_M2 <- fit_branching_point_process(
  #   t = df$t, branching_structure = df$beta, observation_interval = 12,
  #   point_type = 1 + (df$beta > 0), K = 2, omega = omega,
  #   chains = 1, refresh = 0
  # )
  fit_M3 <- fit_branching_point_process(
    t = df$t, branching_structure = df$beta, observation_interval = 12,
    point_type = 1 + (df$beta > 0),
    is_hetero = c(T, F),
    chains = 1, refresh = 0
  )
  expect_true(TRUE)

  expect_error(
    fit_branching_point_process(
      t = df$t, branching_structure = df$beta, observation_interval = 12,
      point_type = 1 + (df$beta > 0), K = 2, omega = rev(omega),
      chains = 1, refresh = 0
    )
  )
})

