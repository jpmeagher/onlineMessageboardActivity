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

test_that("mixture point process likelihood is correct", {
  m_imm <- 0.66
  m_off <- 0.7
  r_imm <- 0.2
  r_off <- 0.4

  compensator <- m_imm * sum(pexp(48 + 8 - sim_cluster$t[1], rate = r_imm))
  compensator <- compensator + m_off * sum(pexp(48 + 8 - sim_cluster$t[-1], rate = r_off))

  imm_off_idx <- sim_cluster$parent_id[-1] == 1
  log_intensity <- sum(
    log(m_imm) +
      dexp(
        sim_cluster$t[-1][imm_off_idx] - sim_cluster$t[sim_cluster$parent_id[-1][imm_off_idx]], rate = r_imm, log = T)
    )
  log_intensity <- log_intensity +
    sum(
      log(m_off) +
        dexp(
          sim_cluster$t[-1][!imm_off_idx] - sim_cluster$t[sim_cluster$parent_id[-1][!imm_off_idx]], rate = r_off, log = T)
    )


  expect_equal(
    log_intensity - compensator,
    mixture_bpp_likelihood(
      t = sim_cluster$t, beta = sim_cluster$parent_id, a = 8 + 48,
      mu_immigrant = m_imm, eta_immigrant = r_imm,
      mu_offspring = m_off, eta_offspring = r_off
    )
  )
})

test_that("time-dependant mixture point process likelihood is correct", {
  m_imm <- 0.66
  m_off <- 0.7
  r_imm <- 0.2
  r_off <- 0.4

  compensator <- m_imm * sum(pexp(48 + 8 - sim_cluster$t[1], rate = r_imm))
  compensator <- compensator + m_imm * sinusoidal_exponential_decay_integral(
    lower_limit = sim_cluster$t[1],
    upper_limit = sim_cluster$t[1] + 48,
    sinusoid_coefficients = alpha,
    sinusoid_frequencies = omega,
    exponential_rate = r_imm
  )
  compensator <- compensator + m_off * sum(pexp(48 + 8 - sim_cluster$t[-1], rate = r_off))
  compensator <- compensator +
    sum(m_off * sinusoidal_exponential_decay_integral(
      lower_limit = sim_cluster$t[-1],
      upper_limit = sim_cluster$t[1] + 48,
      sinusoid_coefficients = alpha,
      sinusoid_frequencies = omega,
      exponential_rate = r_off
    ))

  log_alpha_fun <- log(1 + sinusoidal_function(sim_cluster$t[-1], alpha, omega))

  imm_off_idx <- sim_cluster$parent_id[-1] == 1
  log_rho <- sum(
    log(m_imm) +
      dexp(
        sim_cluster$t[-1][imm_off_idx] - sim_cluster$t[sim_cluster$parent_id[-1][imm_off_idx]], rate = r_imm, log = T)
  )
  log_rho <- log_rho +
    sum(
      log(m_off) +
        dexp(
          sim_cluster$t[-1][!imm_off_idx] - sim_cluster$t[sim_cluster$parent_id[-1][!imm_off_idx]], rate = r_off, log = T)
    )


  expect_equal(
    sum(log_alpha_fun) + log_rho - compensator,
    tide_mixture_bpp_likelihood(
      t = sim_cluster$t, beta = sim_cluster$parent_id, a = 8 + 48,
      mu_immigrant = m_imm, eta_immigrant = r_imm,
      mu_offspring = m_off, eta_offspring = r_off,
      alpha = alpha, omega = omega
    )
  )
})

test_that("Recreate error in model assessment script", {

  tmp_cluster <- sim_cluster[1, ]

  m_imm <- 0.66
  m_off <- 0.7
  r_imm <- 0.2
  r_off <- 0.4

  S <- 10

  par_df <- data.frame(
    `mu[1]` = rnorm(S, mean = m_imm, sd = 0.01), `mu[2]` = rnorm(S, mean = m_off, sd = 0.01),
    `eta[1]` = rnorm(S, mean = r_imm, sd = 0.01), `eta[2]` = rnorm(S, mean = r_off, sd = 0.01),
    `alpha[1]` = rnorm(S, mean = alpha[1], sd = 0.01),
    `alpha[2]` = rnorm(S, mean = alpha[2], sd = 0.01),
    `alpha[3]` = rnorm(S, mean = alpha[3], sd = 0.01),
    `alpha[4]` = rnorm(S, mean = alpha[4], sd = 0.01)
  )
  colnames(par_df) <- c("mu[1]", "mu[2]", "eta[1]", "eta[2]",
                        "alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]")

  alpha_ind <- grepl("alpha", colnames(par_df))
  ll <- mean(sapply(
    1:S,
    function(j){
      tide_mixture_bpp_likelihood(
        t = tmp_cluster$t, beta = tmp_cluster$parent_id, a = 8 + 48,
        mu_immigrant = par_df$`mu[1]`[j], eta_immigrant = par_df$`eta[1]`[j],
        mu_offspring = par_df$`mu[2]`[j], eta_offspring = par_df$`eta[2]`[j],
        alpha = as.numeric(par_df[j, alpha_ind]), omega = omega,
        log = TRUE, perform_checks = T
      )
    }))
  checkmate::expect_number(ll)
})

test_that("time-dependant mixture point process likelihood is correct", {
  m_imm <- 0.66
  m_off <- 0.7
  r_imm <- 0.2
  r_off <- 0.4
  p_imm <- 1
  p_off <- 0.25

  n <- nrow(sim_cluster)
  s <- 1000
  beta_obs <- sim_cluster$parent_id[sim_cluster$t < sim_cluster$t[1] + 3]
  z_obs <- sapply(seq_along(sim_cluster$t), function(i) sum(beta_obs == i))

  post_nu_imm <- gamma_posterior_individual_reproduction_number(
    t = sim_cluster$t[1], observed_offspring = z_obs[1],
    observation_horizon = sim_cluster$t[1] + 3,
    reproduction_number = m_imm, dispersion_parameter = p_imm,
    gi_exp_decay_rate = r_imm,
    sinusoid_coefficients = alpha, sinusoid_frequencies = omega
  )
  post_nu_off <- gamma_posterior_individual_reproduction_number(
    t = sim_cluster$t[-1], observed_offspring = z_obs[-1],
    observation_horizon = sim_cluster$t[1] + 3,
    reproduction_number = m_off, dispersion_parameter = p_off,
    gi_exp_decay_rate = r_off,
    sinusoid_coefficients = alpha, sinusoid_frequencies = omega
  )

  set.seed(101)
  nu <- matrix(NA, nrow = n, ncol = s)
  nu[1, ] <-  rgamma(s, shape = post_nu_imm$shape, rate = post_nu_imm$rate)
  nu[-1, ] <- replicate(s, {
    rgamma(n-1, shape = post_nu_off$shape, rate = post_nu_off$rate)
  })

 tst_lpd <- sapply(
   1:s,
   function(i){
     compensator <- nu[1, i] * sum(pexp(48 + 8 - sim_cluster$t[1], rate = r_imm))
     compensator <- compensator + nu[1, i] * sinusoidal_exponential_decay_integral(
       lower_limit = sim_cluster$t[1],
       upper_limit = sim_cluster$t[1] + 48,
       sinusoid_coefficients = alpha,
       sinusoid_frequencies = omega,
       exponential_rate = r_imm
     )

     compensator <- compensator + sum(nu[-1, i] * pexp(48 + 8 - sim_cluster$t[-1], rate = r_off))
     compensator <- compensator +
       sum(nu[-1, i] * sinusoidal_exponential_decay_integral(
         lower_limit = sim_cluster$t[-1],
         upper_limit = sim_cluster$t[1] + 48,
         sinusoid_coefficients = alpha,
         sinusoid_frequencies = omega,
         exponential_rate = r_off
       ))

     log_alpha_fun <- log(1 + sinusoidal_function(sim_cluster$t[-1], alpha, omega))


     imm_off_idx <- sim_cluster$parent_id[-1] == 1
     log_rho <- sum(
       log(nu[1, i]) +
         dexp(
           sim_cluster$t[-1][imm_off_idx] - sim_cluster$t[sim_cluster$parent_id[-1][imm_off_idx]], rate = r_imm, log = T)
     )
     log_rho <- log_rho +
       sum(
         log(nu[sim_cluster$parent_id[-1][!imm_off_idx], i]) +
           dexp(
             sim_cluster$t[-1][!imm_off_idx] - sim_cluster$t[sim_cluster$parent_id[-1][!imm_off_idx]], rate = r_off, log = T)
       )

     sum(log_alpha_fun) + log_rho - compensator
   })

 expect_equal(
   mean(tst_lpd),
   heterogeneous_tide_mixture_bpp_likelihood(
     t = sim_cluster$t, beta = sim_cluster$parent_id,
     a = sim_cluster$t[1] + 3, b = sim_cluster$t[1] + 48,
     mu_immigrant = m_imm, eta_immigrant = r_imm, psi_immigrant = p_imm,
     mu_offspring = m_off, eta_offspring = r_off, psi_offspring = p_off,
     omega = omega, alpha = alpha,
     log = TRUE, n_samples = 1000, seed = 101,
     perform_checks = TRUE
   )
 )
})
