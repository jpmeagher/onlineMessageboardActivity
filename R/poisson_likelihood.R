

#' Liklihood function for a Poisson process
#'
#' Computes the likelihood for a Poisson process with a periodic intensity function
#'
#' @param t A vector of N positive real values. The time associated with each point.
#' @param a A real valued scalar. Defines the interval over which points arrive.
#' @param lambda_0 The average intensity of arrivals.
#' @param omega The frequencies of the sinusoidal basis functions for the periodic function.
#' @param alpha The coefficients on each sinusoidal basis function
#' @param log Logical. Report the log likelihood?
#' @inheritParams refactor_branching_structure
#'
#' @return A real value. The liklihood for the specified process.
#' @export
periodic_poisson_point_process_likelihood <-
  function(
    t, a,
    lambda_0,
    omega = NULL, alpha = NULL,
    log = TRUE,
    perform_checks = TRUE
  ){
    N <- length(t)
    if (perform_checks) {
      checkmate::assert_numeric(t, lower = 0, any.missing = FALSE)
      checkmate::assert_number(a, lower = max(t))
      checkmate::assert_true(!xor(is.null(alpha), is.null(omega)))
    }
    # compensator
    if (is.null(omega)) {
      Gamma_0 <- lambda_0 * a
    } else {
      Gamma_0 <- lambda_0 * (
        a + sinusoidal_integral(
          lower_limit = 0,
          upper_limit = a,
          sinusoid_coefficients = alpha,
          sinusoid_frequencies = omega
        )
      )
    }
    # intensity function
    l_lambda_fun <- rep(log(lambda_0), N)
    alpha_fun <- 1
    if (!is.null(alpha)) {
      alpha_fun <- alpha_fun +
        sinusoidal_function(
          t, alpha, omega, perform_checks = FALSE
        )
    }
    l_alpha_fun <- log(alpha_fun)
    # log-likelihood
    ll <- sum(l_alpha_fun + l_lambda_fun) - Gamma_0
    if (!log) return(exp(ll))
    ll
  }


