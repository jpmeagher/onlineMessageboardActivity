## code to prepare `DATASET` dataset goes here

## Packages used
library(magrittr)
library(future.apply)
library(dplyr)
library(lubridate)

## Specify that futures will be resolved in parallel
plan("multisession", workers = 6)

## Import the raw_df data
raw_df <- read.csv("~/R/Projects/reddit_data/data/ireland_SR_1549670400_1581206400_comments_v1.csv") %>%
  mutate(Time = as.POSIXct(Time, origin="1970-01-01", tz = "Europe/London")) %>%
  mutate(type = future_sapply(strsplit(Parent, "_"), function(x) x[1])) %>%
  mutate(parent = future_sapply(strsplit(Parent, "_"), function(x) x[2])) %>%
  select(time = Time, type, id = ID, parent)

## Explicitly identify posts and comments to improve clarity
new_type <- rep("comment", nrow(raw_df))
new_type[raw_df$type == "t3"] <- "post"
raw_df <- raw_df %>%
  mutate(type = factor(new_type, levels = c("post", "comment")))
rm(new_type)

## Specify data for analysis
messageboard_df <- raw_df %>%
  filter(dmy(01042019) < time & dmy(01042019) + days(42) > time) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent,
      is_immigrant = (type == "post"),
      S = 2000
    )) %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(tree = factor(identify_clusters(parent_id))) %>%
  select(id, parent_id, tree, time) %>%
  filter(!is.na(tree)) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent_id,
      is_immigrant = (parent_id == 0),
      S = 2000
    )) %>%
  mutate(id = 1:nrow(.))

usethis::use_data(messageboard_df, overwrite = TRUE)


## train/test split
p <- 0.1
set.seed(101)
train_trees <- messageboard_df %>%
  filter(time < dmy(29042019)) %>%
  use_series(tree) %>%
  unique() %>%
  length() %>%
  rbinom(size = 1, p = p) %>%
  as.logical() %>%
  which()

train_df <- messageboard_df %>%
  filter(time < dmy(29042019)) %>%
  filter(tree %in% train_trees) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent_id,
      is_immigrant = (parent_id == 0),
      S = 2000
    )) %>%
  mutate(id = 1:nrow(.))

usethis::use_data(train_df, overwrite = TRUE)

test_trees <- messageboard_df %>%
  filter(time >= dmy(29042019)) %>%
  filter(parent_id == 0) %>%
  use_series(tree) %>%
  as.numeric()

test_df <- messageboard_df %>%
  filter(tree %in% test_trees) %>%
  mutate(
    parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent_id,
      is_immigrant = (parent_id == 0),
      S = 2000
    )) %>%
  mutate(id = 1:nrow(.))

usethis::use_data(test_df, overwrite = TRUE)
