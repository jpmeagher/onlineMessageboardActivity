#' Timestamps and branching structure for online message-board events
#'
#' A data set containing timestamps and branching structure for posts on the
#' \eqn{r/ireland} subreddit seeded during the 42 days (6 weeks) starting on April 1,
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
#' A data set containing 10% of the discussions on the \eqn{r/ireland} subreddit seeded
#' by posts between April 1 and April 21, 2019 inclusive and observed for 48 hours after the
#' initial post. That is, we have a 2,017 discussions sampled at random.
#'
#' @format A dataframe with 5,891 observations of 4 variables:
#' \describe{
#'  \item{id}{A unique, positive integer id for each event ordered by
#'  time}
#'  \item{parent_id}{The id for the parent of each event. Set to 0
#'  when the event is a post and so does not have a parent.}
#'  \item{discussion}{The id for each distinct discussion
#'  seeded by a post.}
#'  \item{t}{The number of hours since 00:00:00 on April 1, 2019.} }
"train_df"

#' Testing data models for online message-board activity
#'
#' A data set containing 10% of the discussions on the \eqn{r/ireland} subreddit
#' seeded  by posts between April 1 and May 12, 2019 inclusive and observed for 48 hours after the
#' initial post. That is, we have a 3,716 discussions sampled at random after excluding discussions
#' already included in `train_df`.
#'
#' @format A dataframe with 11,321 observations of 4 variables:
#' \describe{
#'   #'  \item{id}{A unique, positive integer id for each event ordered by
#'  timestamp.}
#'  \item{parent_id}{The id for the parent of each event. Set to 0
#'  when the event is a post and so does not have a parent.}
#'  \item{discussion}{The id for each distinct discussion
#'  seeded by a post.}
#'  \item{t}{The number of hours since 00:00:00 on April 1, 2019.} }
"test_df"
