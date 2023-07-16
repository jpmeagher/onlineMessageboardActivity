#' Heterogeneous Branching Point Process Marginal Likelihood
#'
#' Computes the heterogeneous branching point process marginal likelihood for
#' the class of branching processes described in Meagher & Friel.
#'
#' @inheritParams branching_point_process_likelihood
#' @param z A non-negative integer-valued vector of length N. The number of
#'   offspring from each point.
#' @param mu_vec A positive real-valued vector of length N. The reproduction
#'   number for each point in the cluster.
#' @param psi_vec A positive real-valued vector of length N. The dispersion
#'   parameter for each point in the cluster.
#' @param eta_vec A positive real-valued vector of length N. The memory decay
#'   rate for each point in the cluster.
#'
#' @return The (log)-likelihood for the branching point process.
#' @export
heterogeneous_branching_point_process_marginal_likelihood <- function(
    t, beta, a,
    mu_vec, psi_vec, eta_vec,
    omega = NULL, alpha = NULL,
    z = NULL,
    log = TRUE,
    perform_checks = TRUE
){
  N <- length(t)
  if (perform_checks) {
    checkmate::assert_numeric(t, any.missing = FALSE)
    checkmate::assert_true(beta[1] == 0)
    checkmate::assert_integerish(beta[-1], len = N - 1, lower = 1, upper = N-1)
    checkmate::assert_numeric(mu_vec, len =  N, lower = 0)
    checkmate::assert_numeric(psi_vec, len =  N, lower = 0)
    checkmate::assert_numeric(eta_vec, len =  N, lower = 0)
    checkmate::assert_number(a, lower = max(t))
    checkmate::assert_integerish(z, lower = 0, len = N, null.ok = T)
    checkmate::assert_true(!xor(is.null(alpha), is.null(omega)))
  }
  inf_idx <- is.infinite(psi_vec)
  # Offspring Counts
  if (is.null(z)) {
    z <- sapply(1:N, function(i) sum(beta == i))
  }
  # Offspring point intensity
  alpha_fun <- 1
  if (!is.null(alpha) & N > 1) {
    alpha_fun <- alpha_fun +
      sinusoidal_function(
        t[-1], alpha, omega, perform_checks = FALSE
      )
  }
  l_alpha_fun <- log(alpha_fun)
  l_rho <- stats::dexp(t[-1] - t[beta[-1]], rate = eta_vec[beta[-1]], log = TRUE)
  l_intensity <- l_alpha_fun + l_rho
  # Rate scaling constants
  K_rate <- stats::pexp(a - t, rate = eta_vec)
  if (!is.null(alpha)) {
    K_rate <- K_rate +
      sinusoidal_exponential_decay_integral(
        t, a, alpha, omega, eta_vec, perform_checks
      )
  }
  # Neg-Bin normalising constants
  l_norm <- lgamma(z + psi_vec) - lgamma(psi_vec) - z * log(mu_vec * K_rate + psi_vec)
  l_norm[inf_idx] <- 0
  # offspring density component
  l_off <- z * log(mu_vec)
  # dispersion density component
  l_disp <- - psi_vec * (log1p(mu_vec * K_rate / psi_vec))
  l_disp[inf_idx] <- - (mu_vec * K_rate)[inf_idx]
  # log-likelihood
  ll <- sum(l_alpha_fun + l_rho) + sum(l_norm + l_disp + l_off)
  if (!log) return(exp(ll))
  ll
}
