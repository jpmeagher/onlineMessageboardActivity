# ##### Set important time-points #####
t_min <- lubridate::ymd_hms(20190401000000, tz = "Europe/London")
tau_threshold_short <- 3
# ## Interesting test case:
tst_tr <- 6
tst_cl <- train_df[train_df$discussion == tst_tr, ] %>%
  dplyr::mutate(
    t = t - min(t)
  ) %>%
  dplyr::mutate(
    parent_id = refactor_branching_structure(
      id = id, parent_id = parent_id, is_immigrant = parent_id == 0, S = nrow(.)
    )
  ) %>%
  dplyr::mutate(
    id = 1:nrow(.)
  ) %>%
  dplyr::filter(t < tau_threshold_short) %>%
  dplyr::select(id, parent_id, t)

day <- 24
f <- 1 / c(day / 2, day, 7 * day)
w <- 2 * pi * f
set.seed(101)
beta <- runif(2 * length(w)) / 10

testthat::test_that("Simulated Hawkes process satisfies sanity checks", {
  set.seed(1011)
  tst <- simulate_hawkes_cluster(
    observed_cluster_df = tst_cl,
    reproduction_number = 0.99,
    gi_exp_decay_rate = 1,
    observation_horizon = tau_threshold_short, simulation_horizon = 2 * tau_threshold_short
  )
  N <- nrow(tst)
  testthat::expect_true(
    all(tst$t[1:(N-1)] < tst$t[2:N])
  )
  testthat::expect_true(
    all(tst$id > tst$parent_id)
  )

  set.seed(1011)
  tst <- simulate_hawkes_cluster(
    observed_cluster_df = tst_cl[1, ],
    reproduction_number = 0.99,
    gi_exp_decay_rate = 1,
    observation_horizon = 0, simulation_horizon = 2
  )
  N <- nrow(tst)
  testthat::expect_true(
    all(tst$t[1:(N-1)] < tst$t[2:N])
  )
  testthat::expect_true(
    all(tst$id > tst$parent_id)
  )
})

testthat::test_that("Simulated DS Hawkes process satisfies sanity checks", {
  set.seed(1011)
  tst <- simulate_ds_hawkes_cluster(
    observed_cluster_df = tst_cl,
    reproduction_number = 0.99, dispersion_parameter = 0.1,
    gi_exp_decay_rate = 1,
    observation_horizon = tau_threshold_short, simulation_horizon = 2 * tau_threshold_short
  )
  N <- nrow(tst)
  testthat::expect_true(
    all(tst$t[1:(N-1)] < tst$t[2:N])
  )
  testthat::expect_true(
    all(tst$id > tst$parent_id)
  )

  set.seed(1012)
  tst <- simulate_ds_hawkes_cluster(
    observed_cluster_df = tst_cl[1,],
    reproduction_number = .99, dispersion_parameter = 100,
    gi_exp_decay_rate = 1,
    observation_horizon = 0, simulation_horizon = 2 * tau_threshold_short
  )
  tst
  N <- nrow(tst)
  testthat::expect_true(
    all(tst$t[1:(N-1)] < tst$t[2:N])
  )
  testthat::expect_true(
    all(tst$id > tst$parent_id)
  )
})

testthat::test_that("Simulated Ti-De Hawkes process satisfies sanity checks", {
  set.seed(1012)
  tst <- simulate_tide_hawkes_cluster(
    observed_cluster_df = tst_cl,
    reproduction_number = 0.99,
    gi_exp_decay_rate = 1,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    observation_horizon = tau_threshold_short, simulation_horizon = 2 * tau_threshold_short
  )
  N <- nrow(tst)
  testthat::expect_true(
    all(tst$t[1:(N-1)] < tst$t[2:N])
  )
  testthat::expect_true(
    all(tst$id > tst$parent_id)
  )

  set.seed(1012)
  tst <- simulate_tide_hawkes_cluster(
    observed_cluster_df = tst_cl[1, ],
    reproduction_number = 0.99,
    gi_exp_decay_rate = 1,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    observation_horizon = 0, simulation_horizon = 2 * tau_threshold_short
  )
  N <- nrow(tst)
  testthat::expect_true(
    all(tst$t[1:(N-1)] < tst$t[2:N])
  )
  testthat::expect_true(
    all(tst$id > tst$parent_id)
  )
})

testthat::test_that("Simulated Ti-De DS Hawkes process satisfies sanity checks", {
  set.seed(1011)
  tst <- simulate_tide_ds_hawkes_cluster(
    observed_cluster_df = tst_cl,
    reproduction_number = 0.99, dispersion_parameter = 0.1,
    gi_exp_decay_rate = 1,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    observation_horizon = tau_threshold_short, simulation_horizon = 2 * tau_threshold_short
  )
  tst
  N <- nrow(tst)
  testthat::expect_true(
    all(tst$t[1:(N-1)] < tst$t[2:N])
  )
  testthat::expect_true(
    all(tst$id > tst$parent_id)
  )

  set.seed(1012)
  tst <- simulate_tide_hawkes_cluster(
    observed_cluster_df = tst_cl[1, ],
    reproduction_number = 0.99,
    gi_exp_decay_rate = 1,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    observation_horizon = 0, simulation_horizon = 2 * tau_threshold_short
  )
  N <- nrow(tst)
  testthat::expect_true(
    all(tst$t[1:(N-1)] < tst$t[2:N])
  )
  testthat::expect_true(
    all(tst$id > tst$parent_id)
  )
})
