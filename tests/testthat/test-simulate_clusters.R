testthat::test_that("Simulated Hawkes process satisfies sanity checks", {
  set.seed(1011)
  tst <- simulate_hawkes_cluster(
    t_i = 0, reproduction_number = 0.99, gi_exp_decay_rate = 1, horizon = 2
  )
  N <- nrow(tst)
  testthat::expect_true(
    all(tst$t[1:(N-1)] < tst$t[2:N])
  )
  testthat::expect_true(
    all(tst$id > tst$parent_id)
  )
  testthat::expect_error(
    simulate_hawkes_cluster(
      t_i = 4, reproduction_number = 0.5, gi_exp_decay_rate = 1, horizon = 2
    )
  )
  testthat::expect_warning(
    simulate_hawkes_cluster(
      t_i = -1, reproduction_number = 1, gi_exp_decay_rate = 1, horizon = 2
    )
  )
  testthat::expect_error(
    simulate_hawkes_cluster(
      t_i = 0, reproduction_number = 0.5, gi_exp_decay_rate = -1, horizon = 2
    )
  )
  testthat::expect_error(
    simulate_hawkes_cluster(
      t_i = 0, reproduction_number = -1, gi_exp_decay_rate = 10, horizon = 2
    )
  )

  simulate_hawkes_cluster(
    t_i = 0,
    reproduction_number = 0.67, gi_exp_decay_rate = 0.27,
    horizon = 3
  )
})
