#' Timestamps and branching structure for online message-board events
#'
#' A dataset containing timestamps and branching structure for activity on the
#' \eqn{r/ireland} subreddit for the 42 days (6 weeks) starting on April 1,
#' 2019
#'
#' @format A dataframe with 117,787 observations of 4 variables:
#' \describe{
#'   \item{id}{A unique, positive integer id for each event ordered by timestamp.}
#'   \item{parent_id}{The id for the parent of each event. Set to 0 when the event is a post.}
#'   \item{tree}{The id for each distinct discussion tree seeded by a post.}
#'   \item{time}{The timestamp for each event.}
#'}
"messageboard_df"

#' Training data for modelling online message-board activity
#'
#' A dataset containing discussion trees on the \eqn{r/ireland} subreddit
#' seeded between April 1 and April 28, 2019 and observed up until April 28,
#' 2019. Each seeded tree is included with probability \eqn{p = 0.1}.
#'
#' @format A dataframe with 8,108 observations of 4 variables:
#' \describe{
#'   \item{id}{A unique, positive integer id for each event ordered by timestamp.}
#'   \item{parent_id}{The id for the parent of each event. Set to 0 when the event is a post.}
#'   \item{tree}{The id for each distinct discussion tree seeded by a post.}
#'   \item{time}{The timestamp for each event.}
#'}
"train_df"

#' Testing data for modelling online message-board activity
#'
#' A dataset containing discussion trees on the \eqn{r/ireland} subreddit
#' seeded between April 29 and May 12, 2019 and observed up until May 12, 2019.
#'
#' @format A dataframe with 37,792 observations of 4 variables:
#' \describe{
#'   \item{id}{A unique, positive integer id for each event ordered by timestamp.}
#'   \item{parent_id}{The id for the parent of each event. Set to 0 when the event is a post.}
#'   \item{tree}{The id for each distinct discussion tree seeded by a post.}
#'   \item{time}{The timestamp for each event.}
#'}
"test_df"
