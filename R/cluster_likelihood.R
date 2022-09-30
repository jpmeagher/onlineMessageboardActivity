#' Branching Point Process Likelihood
#'
#' Computes the branching point process likelihood for the class of branching
#' processes described in Meagher & Friel.
#'
#' @param t A real-valued vector of length N. The times of the point process.
#' @param beta An integer vector of length N. The branching structure of the
#'   point process.
#' @param nu A positive real-valued vector of length N. The fitness associated
#'   with each point of the point process.
#' @param xi A positive real-valued vector of length N. The memory decay
#'   associated with each point of the point process.
#' @param a A real valued scalar. The upper limit of the observation interval.
#' @param alpha A real-valued vector of length 2K. Coefficients of the
#'   sinusoidal basis function for exogenous excitation.
#' @param omega A positive, increasing real-valued vector of length K.
#'   Frequencies defining the sinusoidal basis functions for exogeneous
#'   excitation.
#' @param log A logical scalar. If `TRUE` then the log-likelihood is returned.
#' @param perform_checks A logical scalar. Should inputs be checked before
#'   executing code.
#'
#' @return The (log)-likelihood for the branching point process.
branching_point_process_likelihood <- function(
  t, beta, a,
  nu, xi,
  omega = NULL, alpha = NULL,
  log = TRUE,
  perform_checks = TRUE
  ){

  if (perform_checks) {
    N <- length(t)
    checkmate::assert_numeric(t, any.missing = FALSE)
    checkmate::assert_true(beta[1] == 0)
    checkmate::assert_integerish(beta[-1], len = N - 1, lower = 1, upper = N-1)
    checkmate::assert_numeric(nu, len =  N, lower = 0)
    checkmate::assert_numeric(xi, len = N, lower = 0)
    checkmate::assert_number(a, lower = max(t))
    checkmate::assert_true(!xor(is.null(alpha), is.null(omega)))
  }
  # Compensator
  Psi <- stats::pexp(a - t, rate = xi)
  if (!is.null(alpha)) {
    Psi <- Psi +
      sinusoidal_exponential_decay_integral(
        t, a, alpha, omega, xi, perform_checks
      )
  }
  Psi <- nu * Psi
  # Exogeneous amplitude
  alpha_fun <- 1
  if (!is.null(alpha)) {
    alpha_fun <- alpha_fun +
      sinusoidal_function(
        t[-1], alpha, omega, perform_checks = FALSE
      )
  }
  l_alpha_fun <- log(alpha_fun)
  # Memory kernel
  l_rho <- log(nu[beta[-1]]) + stats::dexp(t[-1] - t[beta[-1]], rate = xi[beta[-1]], log = TRUE)
  # log-likelihood
  ll <- sum(l_alpha_fun, l_rho) - sum(Psi)
  if (!log) return(exp(ll))
  ll
}


#' Homogeneous Branching Point Process Likelihood
#'
#' Computes the branching point process likelihood for the class of branching
#' processes described in Meagher & Friel.
#'
#' @inheritParams branching_point_process_likelihood
#' @param mu A positive real-valued scalar. The reproduction number for each offspring process.
#' @param eta A positive real-valued scalar. The memory decay rate for each offspring process.
#'
#' @return A scalar value. The likelihood for the branching point process.
#' @export
#'
#' @examples
#' # Set a circadian rhythm
#' day <- 24
#' f <- 1 / c(day, day / 2)
#' omega <- 2 * pi * f
#' alpha <- c(-0.17, -0.59, -0.25, 0.34)
#' # Simulate data
#' set.seed(1010)
#' sim_cluster <- simulate_gpm_cluster_process(
#'   t_seed = 8, branching_structure_seed = 0,
#'   observation_horizon = 8 + 0, simulation_horizon = 8 + 48,
#'   immigrant_reproduction_number = 4,
#'   immigrant_gi_exp_decay_rate = 0.2,
#'   immigrant_dispersion_parameter = Inf,
#'   offspring_reproduction_number = 0.66,
#'   offspring_gi_exp_decay_rate = 0.33,
#'   offspring_dispersion_parameter = 0.1,
#'   sinusoid_coefficients = alpha,
#'   sinusoid_frequencies = omega
#' )
#' # Propose some parameters
#' m <- 0.66
#' r <- 0.3
#' # Evaluate the likelihood
#' homogeneous_bpp_likelihood(
#'   t = sim_cluster$t, beta = sim_cluster$parent_id, a = 8 + 48,
#'   mu = m, eta = r
#' )
homogeneous_bpp_likelihood <- function(
    t, beta, a, mu, eta,
    log = TRUE,
    perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_number(mu, lower = 0)
    checkmate::assert_number(eta, lower = 0)
  }
  N <- length(t)
  nu <- rep(mu, N)
  xi <- rep(eta, N)

  branching_point_process_likelihood(
    t = t, beta = beta, a = a,
    nu = nu, xi = xi,
    omega = NULL, alpha = NULL,
    log = log,
    perform_checks = perform_checks
  )
}
