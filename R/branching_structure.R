#' Refactor Branching Structure
#'
#' Refactors a branching structure such that \eqn{e_j = i} if the event at index \eqn{j}
#' is the offspring of the event at index \eqn{i} and \eqn{e_j = 0} otherwise (if the event is an immigrant.)
#' The code assumes that \eqn{i < j}.
#'
#' @param id An N-dimensional vector of unique event ids.
#' @param parent_id An N-dimensional vector noting the id of the parent for each event.
#' @param is_immigrant An N-dimensional logical vector identifying immigrants.
#' @param S An integer scalar such that \eqn{0 < S \leq N}. Refactoring can be more efficient when we consider the most recent
#' S events as possible parents before executing a full sweep.
#' @param perform_checks A logical scalar. Should inputs be checked before executing code.
#'
#' @return The vector of non-negative integers that define the branching structure.
#' @export
refactor_branching_structure <- function(
  id, parent_id, is_immigrant, S,
  perform_checks = TRUE
){
  N <- length(id)
  if (perform_checks) {
    checkmate::assert_atomic_vector(id, any.missing = FALSE, len = N)
    checkmate::assert_true(length(unique(id)) == N)
    checkmate::assert_atomic_vector(parent_id, any.missing = TRUE, len = N)
    checkmate::assert_logical(is_immigrant, any.missing = FALSE, len = N)
    checkmate::assert_integerish(S, lower = 1, upper = N)
  }
  B <- future.apply::future_sapply(
    1:N,
    function(i){
      if(is_immigrant[i]) {
        b_i <- 0
      } else {
        b_i_tmp1 <- which(parent_id[i] == id[max(1, i - S):(i-1)])
        if (length(b_i_tmp1) == 0) {
          b_i_tmp2 <- which(parent_id[i] == id[1:(i-1)])
          if (length(b_i_tmp2) == 0) {
            b_i <- NA
          } else {
            b_i <- b_i_tmp2
          }
        } else {
          b_i <- max(0, i - (S + 1)) + b_i_tmp1
        }
      }
      b_i
    }
  )
  B
}
