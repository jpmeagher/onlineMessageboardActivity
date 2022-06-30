#' Fit a Time-Dependant Doubly-Stochastic Hawkes (TiDe-DSH) process model
#'
#' Sample from the posterior of the TiDe-DSH model given event times and the
#' underlying branching structure.
#'
#' @param t An ordered vector of real-values. The time of each distinct event.
#' @param branching_structure A vector of non-negative integers. Immigrant events identified by 0 while
#' offspring link to the index of their parent.
#' @param a A real-valued scalar. The start of the observation interval.
#' @param b A real-valued scalar. The end of the observation interval.
#' @param is_ds Logical. Is the model doubly stochastic?
#' @param K A non-negative integer. The number of sinusoidal components in the relative activity function,
#' @param f A real valued vector of length K. The frequency of each sinusoidal component.
#' @inheritParams refactor_branching_structure
#' @param pars The parameters return by the sampling algorithm.
#' @param ... Additional arguments for`rstan::sampling`.
#'
#' @export
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
fit_TiDeDSH <- function(
  t, branching_structure,
  a = floor(min(t)), b = ceiling(max(t)),
  is_ds = FALSE,
  K = 0, f = numeric(0),
  perform_checks = TRUE,
  pars = c("mu", "R", "alpha", "eta", "beta", "A", "phase"),
  ...) {
  N <- length(t)
  if (perform_checks) {
    checkmate::assert_numeric(t, lower = a, upper = b, any.missing = FALSE)
    checkmate::assert_integerish(branching_structure, lower = 0, upper = N-1, len = N, any.missing = FALSE)
    checkmate::assert_number(a, upper = min(t))
    checkmate::assert_number(b, lower = max(t))
    checkmate::assert_logical(is_ds)
    checkmate::assert_integerish(K, lower = 0)
    checkmate::assert_numeric(f, sorted = TRUE, len = K, any.missing = FALSE, lower = 0)

  }

  standata <- list(
    N = N,
    N_immigrant = sum(branching_structure == 0),
    N_offspring = sum(branching_structure != 0),
    a = a, b = b,
    t = t,
    v = branching_structure,
    is_ds = is_ds,
    K = K, f = structure(f, dim = K)
  )
  out <- rstan::sampling(
    stanmodels$tidedsh, data = standata,
    pars = pars, ...
  )
  out
}

