
<!-- README.md is generated from README.Rmd. Please edit that file -->

# onlineMessageboardActivity

<!-- badges: start -->
<!-- badges: end -->

The goal of onlineMessageboardActivity is to distribute the data and
algorithms presented in “Modelling online discussion under circadian
rhythms and superspreading dynamics” by Joe Meagher and Nial Friel.

## Installation

You can install the development version of onlineMessageboardActivity
from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jpmeagher/onlineMessageboardActivity")
```

``` r
library(onlineMessageboardActivity)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(rstan)
#> Loading required package: StanHeaders
#> Loading required package: ggplot2
#> rstan (Version 2.21.7, GitRev: 2e1f913d3ca3)
#> For execution on a local, multicore CPU with excess RAM we recommend calling
#> options(mc.cores = parallel::detectCores()).
#> To avoid recompilation of unchanged Stan programs, we recommend calling
#> rstan_options(auto_write = TRUE)
## Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

## Data Examples

Data used in the manuscript is called as follows:

``` r
head(messageboard_df)
#>   id parent_id discussion                time
#> 1  1         0          1 2019-04-01 00:01:52
#> 2  2         0          2 2019-04-01 00:04:13
#> 3  3         0          3 2019-04-01 00:05:03
#> 4  4         0          4 2019-04-01 00:05:11
#> 5  5         0          5 2019-04-01 00:05:48
#> 6  6         0          6 2019-04-01 00:06:38
head(train_df)
#> # A tibble: 6 × 4
#>      id parent_id discussion     t
#>   <int>     <dbl> <fct>      <dbl>
#> 1     1         0 6          0.111
#> 2     2         0 9          0.141
#> 3     3         0 12         0.172
#> 4     4         0 19         0.388
#> 5     5         1 6          0.43 
#> 6     6         5 6          0.474
head(test_df)
#> # A tibble: 6 × 4
#>      id parent_id discussion     t
#>   <int>     <dbl> <fct>      <dbl>
#> 1     1         0 32         0.712
#> 2     2         0 48         1.16 
#> 3     3         0 50         1.20 
#> 4     4         0 55         1.38 
#> 5     5         0 90         3.58 
#> 6     6         0 93         3.75
```

## Algorithm examples

Fit a model as follows

``` r
## Prepare data
train_df <- train_df %>% 
  arrange(discussion) %>% 
  group_by(discussion) %>% 
  mutate(
    within_discussion_parent_id = refactor_branching_structure(
      id = id, parent_id = parent_id, is_immigrant = parent_id == 0, S = length(id)
    )
  ) %>% 
  mutate(within_discussion_id = seq_along(id)) %>% 
  ungroup()
a <- 48
# Specify Hyper-parameters
period <- 24 # 24 hour period
K <- 2
omega <- structure(2 * pi * (1:K / period), dim = K)
## Fit a model
fit <- fit_branching_point_process(
  t = train_df$t, branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a,
  seed = 101
)
# see ?fit_branching_point_process for options
print(fit)
#> Inference for Stan model: branching_point_process.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>            mean se_mean   sd     2.5%      25%      50%      75%    97.5% n_eff
#> mu[1]      0.65    0.00 0.02     0.62     0.64     0.65     0.66     0.69  3119
#> eta[1]     0.39    0.00 0.01     0.37     0.38     0.39     0.40     0.41  2802
#> psi[1]      Inf     NaN  NaN      Inf      Inf      Inf      Inf      Inf   NaN
#> lp__   -4529.34    0.02 0.97 -4532.01 -4529.72 -4529.03 -4528.64 -4528.39  1646
#>        Rhat
#> mu[1]     1
#> eta[1]    1
#> psi[1]  NaN
#> lp__      1
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Dec 20 20:56:42 2022.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

See the supplementary information accompanying the manuscript for more
detail
