test_that("Stan Sampling works", {
  df <- train_df %>%
    dplyr::filter(as.numeric(discussion) < 500) %>%
    dplyr::arrange(discussion) %>%
    dplyr::group_by(discussion) %>%
    dplyr::mutate(tau = t - min(t)) %>%
    dplyr::filter(tau < 12) %>%
    dplyr::mutate(beta = refactor_branching_structure(
      id = id, parent_id = parent_id, is_immigrant = parent_id == 0, S = length(id)
    )) %>%
    dplyr::ungroup()

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

test_that("Periodic point process sampling works", {
  df <- train_df %>%
    dplyr::filter(parent_id == 0)

  days <- 21
  unit <- 24
  omega <- structure(2 * pi * (1:2 / unit), dim = 2)

  fit_M0 <- fit_periodic_point_process(
    t = df$t, immigrant_observation_interval = days * unit,
    chains = 1, refresh = 0
  )

  expect_equal(
    sapply(
      rstan::extract(fit_M0, "lambda_0")[[1]],
      function(x){
        periodic_poisson_point_process_likelihood(
          t = df$t,
          a = days * unit,
          lambda_0 = x
        )
      }
    ),
    rstan::extract(fit_M0, "loglik")[[1]],
    ignore_attr = T
  )

  expect_true(
    quantile(rstan::extract(fit_M0, "lambda_0")[[1]], probs = c(0.25, 0.75))[1] < nrow(df) / (days * unit)
  )

  expect_true(
    quantile(rstan::extract(fit_M0, "lambda_0")[[1]], probs = c(0.25, 0.75))[2] > nrow(df) / (days * unit)
  )

  fit_M1 <- fit_periodic_point_process(
    t = df$t, immigrant_observation_interval = days * unit,
    K = 2, omega = omega,
    chains = 1, refresh = 0
  )

  expect_equal(
    apply(
      rstan::extract(fit_M1, c("lambda_0", 'alpha')) %>%
        as.data.frame(),
      1,
      function(x){
        periodic_poisson_point_process_likelihood(
          t = df$t,
          a = days * unit,
          lambda_0 = x[1],
          alpha = x[-1],
          omega = omega
        )
      }
    ),
    rstan::extract(fit_M1, "loglik")[[1]],
    ignore_attr = T
  )

  expect_true(
    quantile(rstan::extract(fit_M1, "lambda_0")[[1]], probs = c(0.25, 0.75))[1] < nrow(df) / (days * unit)
  )

  expect_true(
    quantile(rstan::extract(fit_M1, "lambda_0")[[1]], probs = c(0.25, 0.75))[2] > nrow(df) / (days * unit)
  )

})
