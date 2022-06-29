#' Simulate a Poisson offspring process
#'
#' Simulate an offspring process for the event at `t_seed`.
#'
#' @param t_seed A real-valued scalar. The time of the seed event.
#' @param sampling_interval An ordered vector of length 2 with real values. The
#'   interval over which we simulate the cluster.
#' @param reproduction_number A positive real valued scalar. The expected number
#'   of offspring from the immigrant event.
#' @param gi_exp_decay_rate A positive real valued scalar. The exponential rate
#'   of decay in the density function of the generation interval distribution
#'   for the immigrant offspring process.
#' @param sinusoid_coefficients A real-valued vector of length \eqn{2 * K}.
#'   \eqn{NULL} default value indicates no absolute-time dependance in the
#'   offspring process.
#' @param sinusoid_frequencies A positive real-valued vector of length \eqn{K}.
#'   ¿\eqn{NULL} default value indicates no absolute-time dependance in the
#'   offspring process. Angular frequencies of the sinusoidal basis.
#' @param dominating_scalar A positive real valued scalar greater than the
#'   absolute-time modulating function for all \eqn{t}. Allows Lewis thinning
#'   for simulating offspring processes. Suitable values are 2 when there is
#'   absolute-time modulation and 2 otherwise.
#' @inheritParams refactor_branching_structure
#'
#' @return An ordered vector of real-valued event times on the sampling interval
#'   for the offspring of t_seed.
simulate_poisson_offspring_process <- function(
  t_seed,
  sampling_interval,
  reproduction_number,
  gi_exp_decay_rate,
  sinusoid_coefficients = NULL, sinusoid_frequencies = NULL,
  dominating_scalar = ifelse(is.null(sinusoid_frequencies), 1, 2),
  perform_checks = TRUE
){
  K <- length(sinusoid_frequencies)
  if (perform_checks) {
    checkmate::assert_number(t_seed, finite = TRUE)
    checkmate::assert_numeric(sampling_interval, lower = t_seed, sorted = TRUE)
    checkmate::assert_number(reproduction_number, lower = 0, finite = TRUE)
    checkmate::assert_number(gi_exp_decay_rate, lower = 0, finite = TRUE)
    checkmate::assert_numeric(sinusoid_coefficients, any.missing = FALSE, len = 2*K, null.ok = TRUE)
    checkmate::assert_numeric(sinusoid_frequencies, any.missing = FALSE, lower = 0, null.ok = TRUE)
    checkmate::assert_number(dominating_scalar, lower = 1, finite = TRUE)
  }

  tmp_F <- stats::pexp(sampling_interval[2] - t_seed, rate = gi_exp_decay_rate) -
    stats::pexp(sampling_interval[1] - t_seed, rate = gi_exp_decay_rate)
  tmp_z <- stats::rpois(n = 1, lambda = dominating_scalar * reproduction_number * tmp_F)
  tmp_t <- t_seed + truncdist::rtrunc(
    n = tmp_z, spec = "exp",
    a = sampling_interval[1] - t_seed, b = sampling_interval[2] - t_seed,
    rate = gi_exp_decay_rate
  )
  if (length(tmp_t) == 0) return(numeric(0))
  tmp_t <- sort(tmp_t)
  if (K == 0) {
    tmp_p <- 1 / dominating_scalar
  } else {
    tmp_p <- (1 + sinusoidal_function(
      t = tmp_t,
      sinusoid_coefficients = sinusoid_coefficients,
      sinusoid_frequencies = sinusoid_frequencies,
      perform_checks = FALSE
      )
    ) / dominating_scalar
  }
  thin_ind <- as.logical(stats::rbinom(tmp_z, size = 1, prob = tmp_p))
  tmp_t[thin_ind]
}

#' Simulate a Gamma-Poisson Mixture Cluster Process
#'
#' Simulate a cluster from the family of generative modes described in Meagher &
#' Friel.
#'
#' @param t_seed A real-valued vector. The time of the seed events.
#' @param branching_structure_seed A vector of non-negative integer values. The
#'   branching structure of seed events.
#' @param observation_horizon A real valued scalar. The time up to which events
#'   in the cluster are observed.
#' @param simulation_horizon A real valued scalar. The time up to which events
#'   in the cluster are simulated.
#' @param immigrant_reproduction_number A positive real valued scalar. The
#'   expected number of offspring from the immigrant event.
#' @param immigrant_gi_exp_decay_rate A positive real valued scalar. The
#'   exponential rate of decay in the density function of the generation
#'   interval distribution for the immigrant offspring process.
#' @param immigrant_dispersion_parameter A positive real valued scalar. The
#'   dispersion parameter for offspring from the immigrant event. Offspring from
#'   the immigrant event are Poisson distributed when this parameter is
#'   Infinite.
#' @param offspring_reproduction_number A positive real valued scalar. The
#'   expected number of offspring from offspring events.
#' @param offspring_gi_exp_decay_rate A positive real valued scalar. The
#'   exponential rate of decay in the density function of the generation
#'   interval distribution for offspring offspring process.
#' @param offspring_dispersion_parameter A positive real valued scalar. The
#'   dispersion parameter for offspring from offspring events. Offspring from
#'   the offspring events are Poisson distributed when this parameter is
#'   Infinite.
#' @param sinusoid_coefficients A real-valued vector of length \eqn{2 * K}.
#'   \eqn{NULL} default value indicates no absolute-time dependance in the
#'   offspring process.
#' @param sinusoid_frequencies A positive real-valued vector of length \eqn{K}.
#'   ¿\eqn{NULL} default value indicates no absolute-time dependance in the
#'   offspring process. Angular frequencies of the sinusoidal basis.
#' @param dominating_scalar A positive real valued scalar greater than the
#'   absolute-time modulating function for all \eqn{t}. Allows Lewis thinning
#'   for simulating offspring processes.
#' @param cluster_id A scalar.An arbitrary cluster label.
#' @inheritParams refactor_branching_structure
#'
#' @return A data frame summarising the simulated Gamma-Poisson mixture cluster
#'   process.
#' @export
simulate_gpm_cluster_process <- function(
    t_seed, branching_structure_seed,
    observation_horizon, simulation_horizon,
    immigrant_reproduction_number,
    immigrant_gi_exp_decay_rate,
    immigrant_dispersion_parameter = Inf,
    offspring_reproduction_number = immigrant_reproduction_number,
    offspring_gi_exp_decay_rate = immigrant_gi_exp_decay_rate,
    offspring_dispersion_parameter = immigrant_dispersion_parameter,
    sinusoid_coefficients = NULL, sinusoid_frequencies = NULL,
    dominating_scalar = ifelse(is.null(sinusoid_frequencies), 1, 2),
    cluster_id = NULL,
    perform_checks = TRUE
){
  N <- length(t_seed)
  K <- length(sinusoid_frequencies)
  if (perform_checks) {
    checkmate::assert_numeric(t_seed, sorted = TRUE, any.missing = FALSE, finite = TRUE)
    checkmate::assert_integerish(branching_structure_seed, lower = 0, any.missing = FALSE)
    checkmate::assert_true(branching_structure_seed[1] == 0)
    checkmate::assert_number(observation_horizon, lower = max(t_seed), finite = TRUE)
    checkmate::assert_number(simulation_horizon, lower = observation_horizon, finite = TRUE)
    checkmate::assert_number(immigrant_reproduction_number, lower = 0, finite = TRUE)
    checkmate::assert_number(immigrant_gi_exp_decay_rate, lower = 0, finite = TRUE)
    checkmate::assert_number(immigrant_dispersion_parameter, lower = 0, finite = FALSE)
    checkmate::assert_number(offspring_reproduction_number, lower = 0, finite = TRUE)
    checkmate::assert_number(offspring_gi_exp_decay_rate, lower = 0, finite = TRUE)
    checkmate::assert_number(offspring_dispersion_parameter, lower = 0, finite = FALSE)
    checkmate::assert_numeric(sinusoid_coefficients, any.missing = FALSE, len = 2*K, null.ok = TRUE)
    checkmate::assert_numeric(sinusoid_frequencies, any.missing = FALSE, lower = 0, null.ok = TRUE)
    checkmate::assert_number(dominating_scalar, lower = 1, finite = TRUE)
  }
  ## Data preparation and transformation
  cluster_df <- data.frame(
    id = 1:N, parent_id = branching_structure_seed, t = t_seed
  )
  cluster_df$z_obs <- sapply(cluster_df$id, function(x) sum(cluster_df$parent_id == x))
  cluster_df$z_sim <- NA
  cluster_df$gen <- 0
  i <- 2
  while (i <= N) {
    cluster_df$gen[i] <- cluster_df$gen[branching_structure_seed[i]] + 1
    i <- i+1
  }
  ## Sample Individual Reproduction numbers
  cluster_df$nu <- NA
  ## Immigrant
  if (is.infinite(immigrant_dispersion_parameter)) {
    cluster_df$nu[1] <- immigrant_reproduction_number
  } else {
    immigrant_post_seed <- gamma_posterior_individual_reproduction_number(
      t = cluster_df$t[1], observed_offspring = cluster_df$z_obs[1],
      observation_horizon = observation_horizon,
      reproduction_number = immigrant_reproduction_number,
      dispersion_parameter = immigrant_dispersion_parameter,
      gi_exp_decay_rate = immigrant_gi_exp_decay_rate,
      sinusoid_coefficients = sinusoid_coefficients, sinusoid_frequencies = sinusoid_frequencies,
      perform_checks = FALSE
    )
    cluster_df$nu[1] <- stats::rgamma(1, immigrant_post_seed$shape, immigrant_post_seed$rate)
  }
  ## Offspring
  if (N > 1) {
    if (is.infinite(offspring_dispersion_parameter)) {
      cluster_df$nu[-1] <- offspring_reproduction_number
    } else {
      offspring_post_seed <- gamma_posterior_individual_reproduction_number(
        t = cluster_df$t[-1], observed_offspring = cluster_df$z_obs[-1],
        observation_horizon = observation_horizon,
        reproduction_number = offspring_reproduction_number,
        dispersion_parameter = offspring_dispersion_parameter,
        gi_exp_decay_rate = offspring_gi_exp_decay_rate,
        sinusoid_coefficients = sinusoid_coefficients, sinusoid_frequencies = sinusoid_frequencies,
        perform_checks = FALSE
      )
      cluster_df$nu[-1] <- stats::rgamma(N-1, offspring_post_seed$shape, offspring_post_seed$rate)
    }
  }
  ## Simulate offspring processes
  ## Immigrant
  t_sim <- simulate_poisson_offspring_process(
    t_seed = cluster_df$t[1],
    sampling_interval = c(observation_horizon, simulation_horizon),
    reproduction_number = cluster_df$nu[1],
    gi_exp_decay_rate = immigrant_gi_exp_decay_rate,
    sinusoid_coefficients = sinusoid_coefficients, sinusoid_frequencies = sinusoid_frequencies,
    dominating_scalar = dominating_scalar,
    perform_checks = FALSE
  )
  M <- length(t_sim)
  cluster_df$z_sim[1] <- M
  nu_sim <- ifelse(
    is.infinite(offspring_dispersion_parameter),
    rep(offspring_reproduction_number, M),
    stats::rgamma(M, shape = offspring_dispersion_parameter, rate = offspring_dispersion_parameter / offspring_reproduction_number)
  )
  if (M > 0) {
    cluster_df <- rbind(
      cluster_df,
      data.frame(
        id = paste(cluster_df$id[1], 1:M, sep = "."),
        parent_id = cluster_df$id[1],
        t = t_sim,
        z_obs = NA, z_sim = NA, gen = cluster_df$gen[i] + 1,
        nu = nu_sim
      )
    )
  }
  ## Offspring
  i <- 2
  while (i <= nrow(cluster_df)) {
    t_sim <- simulate_poisson_offspring_process(
      t_seed = cluster_df$t[i],
      sampling_interval = c(observation_horizon, simulation_horizon),
      reproduction_number = cluster_df$nu[i],
      gi_exp_decay_rate = immigrant_gi_exp_decay_rate,
      sinusoid_coefficients = sinusoid_coefficients, sinusoid_frequencies = sinusoid_frequencies,
      dominating_scalar = dominating_scalar,
      perform_checks = FALSE
    )
    M <- length(t_sim)
    nu_sim <- ifelse(
      is.infinite(offspring_dispersion_parameter),
      rep(offspring_reproduction_number, M),
      stats::rgamma(M, shape = offspring_dispersion_parameter, rate = offspring_dispersion_parameter / offspring_reproduction_number)
      )
    cluster_df$z_sim[i] <- M
    if (M > 0) {
      cluster_df <- rbind(
        cluster_df,
        data.frame(
          id = paste(cluster_df$id[i], 1:M, sep = "."),
          parent_id = cluster_df$id[i],
          t = t_sim,
          z_obs = NA, z_sim = NA, gen = cluster_df$gen[i] + 1,
          nu = nu_sim
          )
      )
    }
    i <- i+1
  }
  ## Bookkeeping to tidy up output
  cluster_df <- dplyr::arrange(cluster_df, t)
  cluster_df$parent_id <- refactor_branching_structure(
    id = cluster_df$id, parent_id = cluster_df$parent_id,
    is_immigrant = cluster_df$parent_id == 0,
    S = nrow(cluster_df), perform_checks = FALSE
  )
  cluster_df$id <- 1:nrow(cluster_df)
  if (!is.null(cluster_id)) {
    cluster_df$id <- paste(cluster_id, cluster_df$id, sep = ".")
    cluster_df$parent_id <- paste(cluster_id, cluster_df$parent_id, sep = ".")
  }
  cluster_df
}
