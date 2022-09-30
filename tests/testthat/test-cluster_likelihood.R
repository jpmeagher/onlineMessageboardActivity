day <- 24
f <- 1 / c(day, day / 2)
omega <- 2 * pi * f
alpha <- c(-0.17, -0.59, -0.25, 0.34)

set.seed(1010)
sim_cluster <- simulate_gpm_cluster_process(
  t_seed = 8, branching_structure_seed = 0,
  observation_horizon = 8 + 0, simulation_horizon = 8 + 48,
  immigrant_reproduction_number = 4,
  immigrant_gi_exp_decay_rate = 0.2,
  immigrant_dispersion_parameter = Inf,
  offspring_reproduction_number = 0.66,
  offspring_gi_exp_decay_rate = 0.33,
  offspring_dispersion_parameter = 0.1,
  sinusoid_coefficients = alpha,
  sinusoid_frequencies = omega
)

test_that("low level likelihood function returns a value and throws correct errors", {
  checkmate::expect_number(
    branching_point_process_likelihood(
      t = sim_cluster$t,
      beta = sim_cluster$parent_id,
      nu = sim_cluster$nu,
      xi = c(0.2, rep(0.33, nrow(sim_cluster) - 1)),
      a = 8 + 48
    )
  )
  # branching structure
  expect_error(
    branching_point_process_likelihood(
      t = sim_cluster$t,
      beta = c(0, sample(100:120, nrow(sim_cluster) - 1)),
      nu = sim_cluster$nu,
      xi = c(0.2, rep(0.33, nrow(sim_cluster) - 1)),
      a = 8 + 48
    )
  )
  # reproduction numbers
  expect_error(
    branching_point_process_likelihood(
      t = sim_cluster$t,
      beta = sim_cluster$parent_id,
      nu = - sim_cluster$nu,
      xi = c(0.2, rep(0.33, nrow(sim_cluster) - 1)),
      a = 8 + 48
    )
  )
  #memory decay
  expect_error(
    branching_point_process_likelihood(
      t = sim_cluster$t,
      beta = sim_cluster$parent_id,
      nu = sim_cluster$nu,
      xi = c(0.2, rep(0.33, nrow(sim_cluster))),
      a = 8 + 48
    )
  )
  expect_error(
    branching_point_process_likelihood(
      t = sim_cluster$t,
      beta = sim_cluster$parent_id,
      nu = sim_cluster$nu,
      xi = c(0.2, rep(-0.33, nrow(sim_cluster) - 1)),
      a = 8 + 48
    )
  )
  # observation interval
  expect_error(
    branching_point_process_likelihood(
      t = sim_cluster$t,
      beta = sim_cluster$parent_id,
      nu = sim_cluster$nu,
      xi = c(0.2, rep(0.33, nrow(sim_cluster) - 1)),
      a = 30
    )
  )

})

test_that("homogeneous point process likelihood is correct", {
  m <- 0.66
  r <- 0.3

  compensator <- m * sum(pexp(48 + 8 - sim_cluster$t, rate = r))
  log_intensity <- sum(log(m) + dexp(sim_cluster$t[-1] - sim_cluster$t[sim_cluster$parent_id[-1]], rate = r, log = T))

  expect_equal(
    log_intensity - compensator,
    homogeneous_bpp_likelihood(
      t = sim_cluster$t, beta = sim_cluster$parent_id, a = 8 + 48,
      mu = m, eta = r
    )
  )
})
