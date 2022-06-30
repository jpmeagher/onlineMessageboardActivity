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
  print(f)
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
