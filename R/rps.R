#' Compute the ranked probability score
#'
#' The ranked probability score (RPS) is a proper scoring rule for evaluating
#' the quality of probabilistic forecasts for ordinal discrete variables.
#'
#' @param obs A positive integer scalar. The observed ordinal variable.
#' @param pred_prob A vector of probabilities. The forecast probability for each
#'   possible outcome.
#' @param n_outcomes A positive integer scalar. The total number of possible
#'   outcomes. The RPS is scaled according to `n_outcomes` and so specifying
#'   this quantity explicitly allows model comparison without increasing the
#'   computational overhead
#' @inheritParams refactor_branching_structure
#'
#' @return A scalar value on the unit interval. Smaller values correspond to
#'   better probabilistic forecasts.
#' @export
compute_rps <- function(
  obs, pred_prob,
  n_outcomes = max(length(pred_prob), obs),
  perform_checks = TRUE
  ){
  n_forecast <- max(length(pred_prob), obs)
  pred_prob <- pred_prob / sum(pred_prob)
  if (perform_checks) {
    checkmate::assert_integerish(obs, len = 1, lower = 1, upper = n_forecast, any.missing = FALSE)
    checkmate::assert_numeric(pred_prob, lower = 0, upper = 1, any.missing = FALSE, max.len = n_outcomes)
    checkmate::assert_integerish(n_outcomes, lower = 1, len = 1)
  }
  ## Temporary computation vectors
  tmp_obs <- tmp_pred_prob <- rep(0, n_forecast)
  tmp_obs[obs] <- 1
  tmp_pred_prob[seq_along(pred_prob)] <- pred_prob
  ## Cumulative computation vectors
  cum_obs <- cumsum(tmp_obs)
  cum_pred_prob <- cumsum(tmp_pred_prob)
  ## RPS
  sum((cum_pred_prob - cum_obs)^2) / (n_outcomes - 1)
}
