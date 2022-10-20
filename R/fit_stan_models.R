#' Fit a Gamma-Poisson-mixture cluster process (GPMCP) model
#'
#' Sample from the posterior of the TGPMCP model given event times and the
#' underlying branching structure.
#'
#' @param t An ordered vector of real-values. The time of each distinct event.
#' @param branching_structure A vector of non-negative integers. Immigrant
#'   events identified by 0 while offspring link to the index of their parent.
#' @param type A vector of positive integers. Assigns each event to one of a set
#'   number of specified types.
#' @param a A real-valued scalar. The start of the observation interval.
#' @param b A real-valued scalar. The end of the observation interval.
#' @param is_hetero Logical. Does the model allow for heterogeneous offspring
#'   processes?
#' @param K A non-negative integer. The number of sinusoidal components in the
#'   relative activity function,
#' @param f A real valued vector of length K. The frequency of each sinusoidal
#'   component.
#' @param sigma_R A positive real-valued scalar. The standard deviation
#'   parameter for the log-normal prior on R.
#' @param sigma_phi A positive real-valued scalar. The standard deviation
#'   parameter for the log-normal prior on R.
#' @param sigma_eta A positive real-valued scalar. The standard deviation
#'   parameter for the log-normal prior on R.
#' @inheritParams refactor_branching_structure
#' @param pars The parameters returned by the sampling algorithm.
#' @param ... Additional arguments for`rstan::sampling`.
#'
#' @export
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
fit_gpmcp <- function(
    t, branching_structure,
    type = rep(1, length(t)),
    a = floor(min(t)), b = ceiling(max(t)),
    is_hetero = FALSE,
    K = 0, f = numeric(0),
    sigma_R = 1, sigma_phi = 0.25, sigma_eta = 2,
    perform_checks = TRUE,
    pars = c("R", "phi", "eta", "beta", "A", "phase"),
    ...) {
  N <- length(t)
  if (perform_checks) {
    checkmate::assert_numeric(t, lower = a, upper = b, any.missing = FALSE)
    checkmate::assert_integerish(branching_structure, lower = 0, upper = N-1, len = N, any.missing = FALSE)
    checkmate::assert_integerish(type, lower=1, len = N)
    checkmate::assert_number(a, upper = min(t))
    checkmate::assert_number(b, lower = max(t))
    checkmate::assert_logical(is_hetero)
    checkmate::assert_integerish(K, lower = 0)
    checkmate::assert_numeric(f, sorted = TRUE, len = K, any.missing = FALSE, lower = 0)
    checkmate::assert_number(sigma_R, lower = 0, finite = TRUE)
    checkmate::assert_number(sigma_phi, lower = 0, finite = TRUE)
    checkmate::assert_number(sigma_eta, lower = 0, finite = TRUE)
  }
  standata <- list(
    N = N, N_1 = sum(branching_structure != 0),
    a = a, b = b,
    t = t, bs = branching_structure,
    N_type = length(unique(type)), type = type,
    is_hetero = is_hetero,
    K = K, f = structure(f, dim = K),
    sigma_R = sigma_R, sigma_phi = sigma_phi, sigma_eta = sigma_eta
  )

  out <- rstan::sampling(
    stanmodels$gpmcp, data = standata,
    pars = pars, ...
  )
  out
}

#' Fit a branching point process model
#'
#' Fits the Branching point process models described in Meagher & Friel to data.
#'
#' @param t A vector of N positive real values. The time associated with each point ordered by cluster.
#' @param branching_structure A vector of N non-negative integers. The underlying branching structure of the point process.
#' Note that this branching structure is defined within each cluster rather for the overall point process.
#' @param observation_interval A positive scalar or vector of length S. The interval over which each point process is observed. Defined either in terms of the end of observation
#' interval for each of the S clusters, or as a scalar defining the observation interval for each cluster from the initial immigrant point.
#' @param point_type An integer vector of length N. The pre-specified type of each point.
#' @param K A non-negative integer. The number of sinusoidal components in the
#'   relative activity function,
#' @param omega A real valued vector of length K. The angular frequency of each sinusoidal
#'   component.
#' @param omega
#' @param is_hetero A logical vector of length M. Do points of each type have heterogeneous offspring
#'   processes.
#' @param shape_mu A positive real-valued scalar. The shape of the Gamma prior on the reproduction number.
#' @param rate_mu A positive real-valued scalar. The rate of the Gamma prior on the reproduction number.
#' @param shape_eta A positive real-valued scalar. The shape of the Gamma prior on the memory decay rate.
#' @param rate_eta A positive real-valued scalar. The rate of the Gamma prior on the memory decay rate.
#' @param ... Additional arguments for`rstan::sampling`.
#'
#' @return A `stanfit` object. Samples from the posterior distribution over model parameters.
#' @export
#'
#' @examples
fit_branching_point_process <- function(
    t, branching_structure, observation_interval,
    point_type = rep(1, length(t)),
    K = 0, omega = numeric(0),
    is_hetero = rep(FALSE, length(unique(point_type))),
    shape_mu = 4, rate_mu = 8,
    shape_eta = 1, rate_eta = 1,
    perform_checks = TRUE,
    ...
  ){
  # Specify vector lengths
  N <- length(t)
  S <- sum(branching_structure == 0)
  M <- length(unique(point_type))
  # Format the observation interval into a vector
  if (length(observation_interval) == 1) observation_interval <- t[branching_structure == 0] + observation_interval
  # Obtain cluster sizes from the branching structure
  cluster_size <- tabulate(cumsum(branching_structure == 0))
  # Check inputs
  if (perform_checks) {
    checkmate::assert_numeric(t, lower = 0, upper = max(observation_interval), any.missing = FALSE)
    checkmate::assert_integerish(branching_structure, lower = 0, upper = N-1, len = N, any.missing = FALSE)
    checkmate::assert_integerish(point_type, lower=1, len = N)
    checkmate::assert_numeric(observation_interval, len = S)
    checkmate::assert_true(all(observation_interval > t[branching_structure == 0]))
    checkmate::assert_logical(is_hetero, len = M)
    checkmate::assert_integerish(K, lower = 0)
    checkmate::assert_numeric(omega, sorted = TRUE, len = K, any.missing = FALSE, lower = 0)
    checkmate::assert_number(shape_mu, lower = 0, finite = TRUE)
    checkmate::assert_number(rate_mu, lower = 0, finite = TRUE)
    checkmate::assert_number(shape_eta, lower = 0, finite = TRUE)
    checkmate::assert_number(rate_eta, lower = 0, finite = TRUE)
  }
  # Clarify variable structure for Stan
  is_hetero <- structure(is_hetero, dim = M)
  omega <- structure(omega, dim = K)
  # Initialising function
  if (K > 0) {
    init_fun <- function(){
      list(
        mu = NULL, log_psi = NULL, eta = NULL,
        alpha = runif(K * 2, min = -1 /  (2 * K), max = 1 /  (2 * K))
      )
    }
  } else {
    init_fun <- function() {
      list(
        mu = NULL, log_psi = NULL, eta = NULL,
        alpha = NULL
      )
    }
  }
  # Data For Stan program
  standata <- list(
    N = N, S = S, K = K, M = M,
    a = observation_interval, omega = omega,
    t = t, beta = branching_structure,
    n_i = cluster_size, type = point_type,
    is_hetero = is_hetero,
    alpha_eta = shape_eta, beta_eta = shape_eta,
    alpha_mu = rate_mu, beta_mu = rate_mu
  )
  # Perform sampling
  out <- rstan::sampling(
    stanmodels$branching_point_process, data = standata, init = init_fun, ...
  )
  out
}

