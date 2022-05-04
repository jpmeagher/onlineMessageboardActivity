#' Individual reproduction number posterior
#'
#' Parameterise the Gamma distributed posterior on individual reproduction
#' numbers.
#'
#' @param t A real-valued vector or scalar. The time points at which the event
#'   arrived.
#' @param observed_offspring A non-negative integer vector or scalar. The number
#'   of offsring from the observed event.
#' @param observation_horizon  A real valued vector or scalar greater than
#'   \eqn{t}. The time up to which offspring are observed,
#' @param reproduction_number A positive real valued scalar. The expected
#'   number of offspring from each event.
#' @param dispersion_parameter A positive real valued scalar. The
#'   dispersion parameter for offspring from each event.
#' @param gi_exp_decay_rate A positive real valued scalar. The exponential rate
#'   of decay in the density function of the generation interval distribution,
#' @inheritParams sinusoidal_exponential_decay_integral
#'
#' @return A named list of length 2 containing shape and rate parameters for the
#'   Gamma posterior distribution on individual reproduction numbers.
#' @export
gamma_posterior_individual_reproduction_number <- function(
  t, observed_offspring,
  observation_horizon,
  reproduction_number, dispersion_parameter,
  gi_exp_decay_rate,
  sinusoid_coefficients = NULL, sinusoid_frequencies = NULL,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_numeric(t, any.missing = FALSE)
    checkmate::assert_integerish(observed_offspring, lower = 0)
    checkmate::assert_numeric(observation_horizon, any.missing = FALSE)
    checkmate::assert_true(all(t < observation_horizon))
    checkmate::assert_number(reproduction_number, lower = 0)
    checkmate::assert_number(dispersion_parameter, lower = 0)
    checkmate::assert_number(gi_exp_decay_rate, lower = 0)
  }
  shape <-  observed_offspring + dispersion_parameter
  rate <- (dispersion_parameter / reproduction_number) +
    stats::pexp(observation_horizon - t, rate = gi_exp_decay_rate)
  if (!is.null(sinusoid_coefficients)){
    rate <- rate + sinusoidal_exponential_decay_integral(
      lower_limit = t, upper_limit = observation_horizon,
      sinusoid_coefficients = sinusoid_coefficients, sinusoid_frequencies = sinusoid_frequencies,
      exponential_rate = gi_exp_decay_rate,
      perform_checks = perform_checks
    )
  }

  list(shape = shape, rate = rate)
}
