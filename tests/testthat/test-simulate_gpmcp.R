test_that("we can simulate Poisson offspring processes", {
  n <- 1000

  R <- 5
  eta <- 0.29

  t0 <- 0
  tau1 <- 0
  tau2 <- Inf

  set.seed(101)
  tst_z <- sapply(
    1:n,
    function(i){
      simulate_poisson_offspring_process(
        t_seed = t0,
        sampling_interval = t0 + c(tau1, tau2),
        reproduction_number = R,
        gi_exp_decay_rate = eta
      ) %>% length()
    })
  expect_equal(
    mean(tst_z), R, tolerance = 2 / sqrt(n)
  )

  tst_gi <- lapply(
    1:n,
    function(i){
      simulate_poisson_offspring_process(
        t_seed = t0,
        sampling_interval = t0 + c(tau1, tau2),
        reproduction_number = R,
        gi_exp_decay_rate = eta
      )
    }) %>%
    unlist()
  expect_equal(
    mean(tst_gi), 1 / eta, tolerance = 2 / sqrt(n)
  )


  t0 <- -10
  tau1 <- qexp(0.25, rate = eta)
  tau2 <- qexp(0.75, rate = eta)

  set.seed(101)
  tst_z <- sapply(
    1:n,
    function(i){
      simulate_poisson_offspring_process(
        t_seed = t0,
        sampling_interval = t0 + c(tau1, tau2),
        reproduction_number = R,
        gi_exp_decay_rate = eta
      ) %>% length()
    })
  expect_equal(
    mean(tst_z), R / 2, tolerance = 2 / sqrt(n)
  )

  tst_gi <- lapply(
    1:n,
    function(i){
      simulate_poisson_offspring_process(
        t_seed = t0,
        sampling_interval = t0 + c(tau1, tau2),
        reproduction_number = R,
        gi_exp_decay_rate = eta
      )
    }) %>%
    unlist()
  expect_true(
    min(tst_gi) > t0 + tau1
  )
  expect_true(
    max(tst_gi) < t0 + tau2
  )

  expect_error(
    simulate_poisson_offspring_process(
      t_seed = "a",
      sampling_interval = t0 + c(tau1, tau2),
      reproduction_number = R,
      gi_exp_decay_rate = eta
    )
  )
  expect_error(
    simulate_poisson_offspring_process(
      t_seed = t0,
      sampling_interval = t0 - c(10, 5),
      reproduction_number = R,
      gi_exp_decay_rate = eta
    )
  )
  expect_error(
    simulate_poisson_offspring_process(
      t_seed = t0,
      sampling_interval = t0 + c(10, 5),
      reproduction_number = R,
      gi_exp_decay_rate = eta
    )
  )
  expect_error(
    simulate_poisson_offspring_process(
      t_seed = t0,
      sampling_interval = t0 + c(tau1, tau2),
      reproduction_number = -1,
      gi_exp_decay_rate = eta
    )
  )
  expect_error(
    simulate_poisson_offspring_process(
      t_seed = t0,
      sampling_interval = t0 + c(tau1, tau2),
      reproduction_number = R,
      gi_exp_decay_rate = -1
    )
  )
})

test_that("we can simulate time-dependant Poisson offspring processes", {
  n <- 1000
  tst_t <- seq(0, 48, length.out = 1001)

  R <- 5
  eta <- 0.29

  t0 <- 0 # Try different values here
  tau1 <- 0
  tau2 <- 100

  day <- 24
  f <- 1 / c(day / 2, day)
  w <- 2 * pi * f
  set.seed(101)
  beta <- c(-0.17, -0.59, -0.25, 0.34)

  set.seed(101)
  tst_z <- sapply(
    1:n,
    function(i){
      simulate_poisson_offspring_process(
        t_seed = t0,
        sampling_interval = t0 + c(tau1, tau2),
        reproduction_number = R,
        gi_exp_decay_rate = eta,
        sinusoid_coefficients = beta, sinusoid_frequencies = w
      ) %>% length()
    })
  sin_exp_int <- pexp(tau2, rate = eta) + sinusoidal_exponential_decay_integral(
    lower_limit = t0 + tau1, upper_limit = t0 + tau2,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    exponential_rate = eta
  )
  expect_equal(
    mean(tst_z), R * sin_exp_int , tolerance = 2 / sqrt(n)
  )

  tst_gi <- lapply(
    1:n,
    function(i){
      simulate_poisson_offspring_process(
        t_seed = t0,
        sampling_interval = t0 + c(tau1, tau2),
        reproduction_number = R,
        gi_exp_decay_rate = eta,
        sinusoid_coefficients = beta, sinusoid_frequencies = w
      )
    }) %>%
    unlist()

  # Visual inspection of the generation interval distribution looks good
  # hist(tst_gi, breaks = 100, freq = F)
  # lines(
  #   tst_t,
  #   (1 + sinusoidal_function(tst_t, sinusoid_coefficients = beta, sinusoid_frequencies = w)) * dexp(tst_t - t0, rate = eta) / sin_exp_int,
  #   col = "red")

  t0 <- -10
  tau1 <- qexp(0.25, rate = eta)
  tau2 <- qexp(0.75, rate = eta)

  tst_gi <- lapply(
    1:n,
    function(i){
      simulate_poisson_offspring_process(
        t_seed = t0,
        sampling_interval = t0 + c(tau1, tau2),
        reproduction_number = R,
        gi_exp_decay_rate = eta,
        sinusoid_coefficients = beta, sinusoid_frequencies = w
      )
    }) %>%
    unlist()
  expect_true(
    min(tst_gi) > t0 + tau1
  )
  expect_true(
    max(tst_gi) < t0 + tau2
  )
})

test_that("simulating hawkes from gpmcp works", {
  # ##### Set important time-points #####
  t_min <- lubridate::ymd_hms(20190401000000, tz = "Europe/London")
  tau1 <- 3
  # ## Interesting test case:
  tst_tr <- 1415
  tst_cl <- train_df[train_df$tree == tst_tr, ] %>%
    dplyr::mutate(
      t = difftime(time, t_min, units = "hours") %>% as.numeric()
    ) %>%
    dplyr::mutate(
      rel_t = t - min(t)
    ) %>%
    dplyr::mutate(
      parent_id = refactor_branching_structure(
        id = id, parent_id = parent_id, is_immigrant = parent_id == 0, S = nrow(.)
      )
    ) %>%
    dplyr::mutate(
      id = 1:nrow(.)
    ) %>%
    dplyr::filter(rel_t < tau1) %>%
    dplyr::select(id, parent_id, t, rel_t)

  R <- 2
  eta <- 0.29

  set.seed(1011)
  tst <- simulate_gpm_cluster_process(
    t_seed = tst_cl$t, branching_structure_seed = tst_cl$parent_id,
    observation_horizon = tst_cl$t[1] + tau1, simulation_horizon = tst_cl$t[1] + 2 * tau1,
    immigrant_reproduction_number = R,
    immigrant_gi_exp_decay_rate = eta
  )
  expect_true(all(tst$nu == R))
  expect_true(all(diff(tst$t) > 0))

  set.seed(1011)
  tst <- simulate_gpm_cluster_process(
    t_seed = tst_cl$t, branching_structure_seed = tst_cl$parent_id,
    observation_horizon = tst_cl$t[1] + tau1, simulation_horizon = tst_cl$t[1] + 2 * tau1,
    immigrant_reproduction_number = R,
    immigrant_gi_exp_decay_rate = eta,
    immigrant_dispersion_parameter = 1
  )
  expect_true(all(tst$nu != R))
  expect_true(all(diff(tst$t) > 0))

  set.seed(1011)
  tst <- simulate_gpm_cluster_process(
    t_seed = tst_cl$t, branching_structure_seed = tst_cl$parent_id,
    observation_horizon = tst_cl$t[1] + tau1, simulation_horizon = tst_cl$t[1] + 2 * tau1,
    immigrant_reproduction_number = R,
    immigrant_gi_exp_decay_rate = eta,
    # immigrant_dispersion_parameter = 1,
    offspring_reproduction_number = R * 0.99,
    offspring_gi_exp_decay_rate = 0.8 * eta
  )
  expect_true(all(tst$nu[-1] == R * 0.99))
  expect_true(all(diff(tst$t) > 0))

  set.seed(1011)
  tst <- simulate_gpm_cluster_process(
    t_seed = tst_cl$t, branching_structure_seed = tst_cl$parent_id,
    observation_horizon = tst_cl$t[1] + tau1, simulation_horizon = tst_cl$t[1] + 2 * tau1,
    immigrant_reproduction_number = R,
    immigrant_gi_exp_decay_rate = eta,
    immigrant_dispersion_parameter = 1,
    offspring_reproduction_number = R * 0.99,
    offspring_gi_exp_decay_rate = 0.8 * eta,
    offspring_dispersion_parameter = 1.5
  )
  expect_true(all(tst$nu[-1] != R * 0.99))
  expect_true(all(diff(tst$t) > 0))

  day <- 24
  f <- 1 / c(day / 2, day)
  w <- 2 * pi * f
  beta <- c(-0.17, -0.59, -0.25, 0.34)

  tst <- simulate_gpm_cluster_process(
    t_seed = tst_cl$t, branching_structure_seed = tst_cl$parent_id,
    observation_horizon = tst_cl$t[1] + tau1, simulation_horizon = tst_cl$t[1] + 2 * tau1,
    immigrant_reproduction_number = R,
    immigrant_gi_exp_decay_rate = eta,
    immigrant_dispersion_parameter = 1,
    offspring_reproduction_number = R * 0.99,
    offspring_gi_exp_decay_rate = 0.8 * eta,
    offspring_dispersion_parameter = 1.5,
    sinusoid_coefficients = beta,
    sinusoid_frequencies = w
  )
  expect_true(all(tst$nu[-1] != R * 0.99))
  expect_true(all(diff(tst$t) > 0))

  tst <- simulate_gpm_cluster_process(
    t_seed = tst_cl$t, branching_structure_seed = tst_cl$parent_id,
    observation_horizon = tst_cl$t[1] + tau1, simulation_horizon = tst_cl$t[1] + 2 * tau1,
    immigrant_reproduction_number = R,
    immigrant_gi_exp_decay_rate = eta,
    immigrant_dispersion_parameter = 1,
    offspring_reproduction_number = R * 0.99,
    offspring_gi_exp_decay_rate = 0.8 * eta,
    offspring_dispersion_parameter = 1.5,
    sinusoid_coefficients = beta,
    sinusoid_frequencies = w,
    cluster_id = "test"
  )
  expect_true(all(tst$nu[-1] != R * 0.99))
  expect_true(all(diff(tst$t) > 0))
})
