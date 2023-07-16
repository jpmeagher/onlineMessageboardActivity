## code to prepare data sets goes here

## Packages used
library(magrittr)
library(future.apply)
library(dplyr)
library(lubridate)
library(onlineMessageboardActivity)

## Specify that futures will be resolved in parallel
plan("multisession", workers = 6)

## Import the raw_df data
raw_df <- read.csv("~/R/Projects/reddit_data/data/ireland_SR_1549670400_1581206400_comments_v1.csv") %>%
  mutate(Time = as.POSIXct(Time, origin="1970-01-01", tz = "Europe/London")) %>%
  mutate(type = future_sapply(strsplit(Parent, "_"), function(x) x[1])) %>%
  mutate(parent = future_sapply(strsplit(Parent, "_"), function(x) x[2])) %>%
  select(time = Time, type, id = ID, parent)

## Specify data for analysis
## We use 6 weeks of discussions
## We filter for nodes occurring in the 44 days following April 1, 2019
## This allows us to observe the first 48 hours of all relevant discussions
## We filter out nodes belonging to discussions seeded outside the period

# all nodes in the 44 days
messageboard_df <- raw_df %>%
  filter(dmy(01042019) < time & dmy(01042019) + days(44) > time) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent,
      is_immigrant = (type == "t3"),
      S = 2000
    )) %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(discussion = factor(identify_clusters(parent_id))) %>%
  select(id, parent_id, discussion, time)

# filter out nodes from discussions seeded before April 1
messageboard_df <- messageboard_df %>%
  filter(!is.na(discussion)) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent_id,
      is_immigrant = (parent_id == 0),
      S = 2000
    )) %>%
  mutate(id = 1:nrow(.))

## identify discussions seeded after the 6 week period
late_discussions <- messageboard_df %>%
  filter(parent_id == 0) %>%
  filter(time > dmy(01042019) + days(42)) %>%
  use_series(discussion)

# filter out late discussions
messageboard_df <- messageboard_df %>%
  filter(!(discussion %in% late_discussions)) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent_id,
      is_immigrant = (parent_id == 0),
      S = 2000
    )) %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(discussion = factor(discussion))

usethis::use_data(messageboard_df, overwrite = TRUE)

## Training data
max_tau <- 48 # 48 hours
set.seed(123456)

# Train on 10% of discussions seeded in the first 21 days
S_training_total <- messageboard_df %>%
  filter(time < dmy(01042019) + days(21)) %>%
  use_series(discussion) %>%
  unique() %>%
  as.numeric() %>%
  length()
S_training <- round(0.1 * S_training_total)
training_discussions <- sort(sample.int(S_training_total, S_training))

train_df <- messageboard_df %>%
  filter(discussion %in% training_discussions) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id, parent_id = parent_id,
      is_immigrant = parent_id == 0,
      S = 250
    )
  ) %>%
  mutate(id = seq_along(id)) %>%
  mutate(t = as.numeric(difftime(time, dmy_hms(010419000000, tz = "Europe/London"), units = "hours"))) %>%
  group_by(discussion) %>%
  mutate(tau = t - min(t)) %>%
  ungroup() %>%
  filter(tau < max_tau) %>%
  select(id, parent_id, discussion, t) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id, parent_id = parent_id,
      is_immigrant = parent_id == 0,
      S = 250
    )
  ) %>%
  mutate(id = seq_along(id))

usethis::use_data(train_df, overwrite = TRUE)

## Testing Data

# Test on 10% of discussions over the full 6 weeks
S_testing_total <- messageboard_df %>%
  use_series(discussion) %>%
  unique() %>%
  as.numeric() %>%
  length()

S_testing <- round(0.1 * S_testing_total)
testing_discussions <- sort(sample(seq.int(S_testing_total)[-training_discussions], S_testing))

test_df <- messageboard_df %>%
  filter(discussion %in% testing_discussions) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id, parent_id = parent_id,
      is_immigrant = parent_id == 0,
      S = 250
    )
  ) %>%
  mutate(id = seq_along(id)) %>%
  mutate(t = as.numeric(difftime(time, dmy_hms(010419000000, tz = "Europe/London"), units = "hours"))) %>%
  select(id, parent_id, discussion, t) %>%
  group_by(discussion) %>%
  mutate(tau = t - min(t)) %>%
  ungroup() %>%
  filter(tau < max_tau) %>%
  select(id, parent_id, discussion, t) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id, parent_id = parent_id,
      is_immigrant = parent_id == 0,
      S = 250
    )
  ) %>%
  mutate(id = seq_along(id))

usethis::use_data(test_df, overwrite = TRUE)
