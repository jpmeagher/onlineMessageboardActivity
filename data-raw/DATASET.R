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

# ## Explicitly identify posts and comments to improve clarity
# new_type <- rep("comment", nrow(raw_df))
# new_type[raw_df$type == "t3"] <- "post"
# raw_df <- raw_df %>%
#   mutate(type = factor(new_type, levels = c("post", "comment")))
# rm(new_type)


## Specify data for analysis
messageboard_df <- raw_df %>%
  filter(dmy(01042019) < time & dmy(01042019) + days(42) > time) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent,
      is_immigrant = (type == "t3"),
      S = 2000
    )) %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(discussion = factor(identify_clusters(parent_id))) %>%
  select(id, parent_id, discussion, time) %>%
  filter(!is.na(discussion)) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent_id,
      is_immigrant = (parent_id == 0),
      S = 2000
    )) %>%
  mutate(id = 1:nrow(.))

usethis::use_data(messageboard_df, overwrite = TRUE)

## Training data
max_tau <- 48
set.seed(123456)

S_training_total <- messageboard_df %>%
  filter(time < dmy(08042019)) %>%
  use_series(discussion) %>%
  unique() %>%
  as.numeric() %>%
  max()
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
  select(id, parent_id, discussion, t)

usethis::use_data(train_df, overwrite = TRUE)

## Testing Data

S_testing_total <- messageboard_df %>%
  filter(time < dmy(11052019)) %>%
  use_series(discussion) %>%
  unique() %>%
  as.numeric() %>%
  max()

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
  select(id, parent_id, discussion, t)

usethis::use_data(test_df, overwrite = TRUE)

# messageboard_df %>%
#   mutate(t = as.numeric(difftime(time, dmy_hms(010419000000, tz = "Europe/London"), units = "hours"))) %>%
#   group_by(discussion) %>%
#   mutate(tau = t - min(t)) %>%
#   filter(tau > max_tau) %>%
#   nrow()
#
# messageboard_df %>%
#   mutate(t = as.numeric(difftime(time, dmy_hms(010419000000, tz = "Europe/London"), units = "hours"))) %>%
#   group_by(discussion) %>%
#   mutate(tau = t - min(t)) %>%
#   filter(tau > max_tau) %>%
#   use_series(discussion) %>%
#   unique() %>%
#   length()
#
#
#
# 1191 / sum(messageboard_df$parent_id != 0)
# 512 / sum(messageboard_df$parent_id == 0)
# 512 / S_testing_total
