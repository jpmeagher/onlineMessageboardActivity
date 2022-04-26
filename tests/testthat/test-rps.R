test_that("Check bog standard rps", {
  y <- 3
  prob <- dnbinom(1:20, size = 0.1, mu = 10)
  prob_sc <- prob / sum(prob)
  ## verification package does not scale probabilities automatically
  expect_equal(
    compute_rps(obs = y, pred_prob = prob),
    verification::rps(obs = y, pred = t(prob_sc))$rps
  )

})

test_that("Check rps when outcome outside forecast range", {
  y <- 40
  prob <- dnbinom(1:20, size = 0.1, mu = 10)
  prob_sc <- rep(0, y)
  prob_sc[seq_along(prob)] <- prob / sum(prob)

  expect_equal(
    compute_rps(obs = y, pred_prob = prob),
    verification::rps(obs = y, pred = t(prob_sc))$rps
  )
})

test_that("Check rps when we expand the forecast range", {
  y <- 40
  prob <- dnbinom(1:20, size = 0.1, mu = 10)
  prob_sc <- rep(0, 1000)
  prob_sc[seq_along(prob)] <- prob / sum(prob)

  expect_equal(
    compute_rps(obs = y, pred_prob = prob, n_outcomes = 1000),
    verification::rps(obs = y, pred = t(prob_sc))$rps
  )
})
