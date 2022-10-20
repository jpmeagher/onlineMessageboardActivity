day <- 24
f <- 1 / c(day, day / 2)
omega <- 2 * pi * f
alpha <- c(-0.17, -0.59, -0.25, 0.34)

a <- 48

set.seed(1010)
sim_cluster <- simulate_gpm_cluster_process(
  t_seed = 8, branching_structure_seed = 0,
  observation_horizon = 8 + 0, simulation_horizon = 8 + a,
  immigrant_reproduction_number = 4,
  immigrant_gi_exp_decay_rate = 0.2,
  immigrant_dispersion_parameter = Inf,
  offspring_reproduction_number = 0.66,
  offspring_gi_exp_decay_rate = 0.33,
  offspring_dispersion_parameter = 0.1,
  sinusoid_coefficients = alpha,
  sinusoid_frequencies = omega
)
N <- nrow(sim_cluster)

test_that("Likelihood function works", {
  m <- 0.66
  p <- 7
  r <- 0.3

  expect_equal(
    heterogeneous_branching_point_process_marginal_likelihood(
      t = sim_cluster$t, beta = sim_cluster$parent_id,
      a = sim_cluster$t[1] + a,
      mu_vec = rep(m, N), psi_vec = rep(Inf, N), eta_vec = rep(r, N),
    ),
    branching_point_process_likelihood(
      t = sim_cluster$t, beta = sim_cluster$parent_id,
      a = sim_cluster$t[1] + a,
      nu = rep(m, N), xi = rep(r, N),
    )
  )

  expect_equal(
    heterogeneous_branching_point_process_marginal_likelihood(
      t = sim_cluster$t, beta = sim_cluster$parent_id,
      a = sim_cluster$t[1] + a,
      mu_vec = rep(m, N), psi_vec = rep(Inf, N), eta_vec = rep(r, N),
      alpha = alpha, omega = omega
    ),
    branching_point_process_likelihood(
      t = sim_cluster$t, beta = sim_cluster$parent_id,
      a = sim_cluster$t[1] + a,
      nu = rep(m, N), xi = rep(r, N),
      alpha = alpha, omega = omega
    )
  )



  # Estimate ML by MC estimate
  # Take mean likelihood under the prior.
  # Very unstable, but not obviously wrong.

  # ll <- heterogeneous_branching_point_process_marginal_likelihood(
  #   t = sim_cluster$t, beta = sim_cluster$parent_id,
  #   a = a,
  #   mu_vec = rep(m, N), psi_vec = rep(p, N), eta_vec = rep(r, N),
  #   omega = omega, alpha = alpha,
  #   log = TRUE,
  #   perform_checks = TRUE
  # )
  # s <- 10000
  # ll_vec <- replicate(s, {
  #   nu <- rgamma(N, shape = p, rate = p/m)
  #   branching_point_process_likelihood(
  #     t = sim_cluster$t, beta = sim_cluster$parent_id,
  #     a = a,
  #     nu = nu, xi = rep(r, N),
  #     omega = omega, alpha = alpha
  #   )
  # })
  #
  # cum_ml <- sapply(
  #   1:s,
  #   function(i) matrixStats::logSumExp(ll_vec[1:i]) - log(i)
  #   )
  #
  # plot(cum_ml, type = "l")
  # abline(h = ll, col = "red")



})
