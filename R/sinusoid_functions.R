#' Linear Combination of Sinusoidal Basis Functions
#'
#' Compute the value of the function \deqn{\rho (t) = \sum_{k = 1}^K \beta_{2k -
#' 1} sin(\omega_k t ) + \beta_{2k} cos(\omega_k t )}.
#'
#' @param t A real-valued vector or scalar. The time points at which the
#'   function is computed.
#' @param sinusoid_coefficients A real-valued vector of length \eqn{2 * K}.
#'   Denoted \eqn{\beta} above.
#' @param sinusoid_frequencies A positive real-valued vector of length \eqn{K}.
#'   Angular frequencies of the sinusoidal basis. Denoted \eqn{\omega} above.
#' @inheritParams refactor_branching_structure
#'
#' @return A vector of real-valued function evaluations, one for each element of
#'   `t`.
#' @export
sinusoidal_function <- function(
  t,
  sinusoid_coefficients, sinusoid_frequencies,
  perform_checks = TRUE
){
  N <- length(t)
  K <- length(sinusoid_frequencies)
  if (perform_checks) {
    checkmate::assert_numeric(t, any.missing = FALSE, finite = TRUE)
    checkmate::assert_numeric(sinusoid_coefficients, any.missing = FALSE, len = 2*K)
    checkmate::assert_numeric(sinusoid_frequencies, any.missing = FALSE, lower = 0)
  }
  sinusoidal_basis <- matrix(
    sapply(
      sinusoid_frequencies, function(x){
        c(sin(t * x), cos(t * x))
      }), nrow = N, ncol = 2*K)

  c(sinusoidal_basis %*% sinusoid_coefficients)
}


#' Integrate a Linear Combination of Sinusoidal Basis Functions
#'
#' Compute the integral of the function \deqn{\int_{a}^{b} \rho (t) dt =
#' \int_{a}^{b} \sum_{k = 1}^K \beta_{2k - 1} sin(\omega_k t ) + \beta_{2k}
#' cos(\omega_k t ) dt}.
#'
#' @param lower_limit A real-valued vector or scalar. The lower limit of the
#'   integral. The quantity a above.
#' @param upper_limit A real-valued vector or scalar. The upper limit of the
#'   integral. The quantity b above.
#' @inheritParams sinusoidal_function
#'
#' @return A vector of real-valued function evaluations, one for each element of
#'   the longer of `lower_limit` or `upper_limit`.
#' @export
sinusoidal_integral <- function(
  lower_limit, upper_limit,
  sinusoid_coefficients, sinusoid_frequencies,
  perform_checks = TRUE
){
  N <- max(length(lower_limit), length(upper_limit))
  K <- length(sinusoid_frequencies)
  if (perform_checks) {
    checkmate::assert_numeric(lower_limit, any.missing = FALSE, finite = TRUE)
    checkmate::assert_numeric(upper_limit, any.missing = FALSE, finite = TRUE)
    checkmate::assert_true(all(lower_limit <= upper_limit))
    checkmate::assert_numeric(sinusoid_coefficients, any.missing = FALSE, len = 2*K)
    checkmate::assert_numeric(sinusoid_frequencies, any.missing = FALSE, lower = 0)
  }
  integral_basis <- matrix(
    sapply(
      sinusoid_frequencies, function(x){
        c((cos(x * lower_limit) - cos(x * upper_limit)) / x, (sin(x * upper_limit) - sin(x * lower_limit)) / x)
      }), nrow = N, ncol = 2*K)

  c(integral_basis %*% sinusoid_coefficients)
}


#' Integrate a Linear Combination of Exponentially decaying Sinusoidal Basis
#' Functions
#'
#' Compute the integral of the function \deqn{\int_{a}^{b} \rho (t) \exp(- \psi
#' (t - a)) dt = \int_{a}^{b} \exp(- \psi (t - a)) \left( \sum_{k = 1}^K
#' \beta_{2k - 1} sin(\omega_k t ) + \beta_{2k} cos(\omega_k t ) \right) dt}.
#'
#' @inheritParams sinusoidal_integral
#' @param exponential_rate A positive real-valued scalar. The exponential rate
#'   of decay. \eqn{\psi} in the expression above.
#'
#' @return A vector of real-valued function evaluations, one for each element of
#'   the longer of `lower_limit` or `upper_limit`.
#' @export
sinusoidal_exponential_decay_integral <- function(
  lower_limit, upper_limit,
  sinusoid_coefficients, sinusoid_frequencies, exponential_rate,
  perform_checks = TRUE
){
  N <- max(length(lower_limit), length(upper_limit))
  K <- length(sinusoid_frequencies)
  if (perform_checks) {
    checkmate::assert_numeric(lower_limit, any.missing = FALSE, finite = TRUE)
    checkmate::assert_numeric(upper_limit, any.missing = FALSE, finite = TRUE)
    checkmate::assert_true(all(lower_limit <= upper_limit))
    checkmate::assert_numeric(sinusoid_coefficients, any.missing = FALSE, len = 2*K)
    checkmate::assert_numeric(sinusoid_frequencies, any.missing = FALSE, lower = 0)
    checkmate::assert_number(exponential_rate, lower = 0)
  }
  integral_basis <- matrix(
    sapply(
      sinusoid_frequencies, function(x){
        c(
          (exponential_rate / (exponential_rate^2 + x^2)) * (
            exponential_rate * (
                sin(x * lower_limit) - (sin(x * upper_limit) * exp(-exponential_rate * (upper_limit - lower_limit)))
                ) +
              x * (
                cos(x * lower_limit) - (cos(x * upper_limit) * exp(-exponential_rate * (upper_limit - lower_limit)))
              )
          ),
          (exponential_rate / (exponential_rate^2 + x^2)) * (
            x * (
              (sin(x * upper_limit) * exp(-exponential_rate * (upper_limit - lower_limit))) - sin(x * lower_limit)
              ) -
            exponential_rate * (
              (cos(x * upper_limit) * exp(-exponential_rate * (upper_limit - lower_limit))) - cos(x * lower_limit)
            )
          )
        )
      }), nrow = N, ncol = 2*K)

  c(integral_basis %*% sinusoid_coefficients)
}
