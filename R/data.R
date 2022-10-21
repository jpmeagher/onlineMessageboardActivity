#' Timestamps and branching structure for online message-board events
#'
#' A dataset containing timestamps and branching structure for activity on the
#' \eqn{r/ireland} subreddit for the 42 days (6 weeks) starting on April 1,
#' 2019. Consists of 117,787 events within 38,616 converstaions.
#'
#' @format A dataframe with 117,787 observations of 4 variables:
#' \describe{
#'   \item{id}{A unique, positive integer id for each event ordered by timestamp.}
#'   \item{parent_id}{The id for the parent of each event. Set to 0 when the event is a post.}
#'   \item{discussion}{The id for each distinct discussion seeded by a post.}
#'   \item{time}{The timestamp for each event.}
#'}
"messageboard_df"

#' Training data for models of online message-board activity
#'
#' A dataset containing discussions on the \eqn{r/ireland} subreddit seeded
#' by posts between April 1 and April 28, 2019 and observed for 48 hours after the
#' initial post. Each discussion seeded during the period is included with
#' probability \eqn{p = 0.1} for a total of 704 discussions consisting of 2262
#' events.
#'
#' @format A dataframe with 8,108 observations of 4 variables: \describe{
#'  \item{id}{A unique, positive integer id for each event ordered by
#'  timestamp.}
#'  \item{parent_id}{The id for the parent of each event. Set to 0
#'  when the event is a post.}
#'  \item{discussion}{The id for each distinct discussion
#'  seeded by a post.}
#'  \item{time}{The timestamp for each event.} }
"train_df"

#' Testing data models od online message-board activity
#'
#' A dataset containing discussion on the \eqn{r/ireland} subreddit
#' seeded between April 29 and May 10, 2019 and observed for 48 hours after the
#' initial post. Initially, each discussion seeded during the period is included with
#' probability \eqn{p = 0.1}, but those that had already appeared in the training data are
#' then excluded. This provides a total of 3,601 discussions consisting of 10,864
#' events.
#'
#' @format A dataframe with 37,2022 observations of 4 variables:
#' \describe{
#'   \item{id}{A unique, positive integer id for each event ordered by timestamp.}
#'   \item{parent_id}{The id for the parent of each event. Set to 0 when the event is a post.}
#'   \item{discussion}{The id for each distinct discussion seeded by a post.}
#'   \item{time}{The timestamp for each event.}
#'}
"test_df"
