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

test_that("PWM CRPS seems to work well", {
  y <- 3
  z <- rnbinom(100, size = 0.1, mu = 3)

  expect_equal(
    mean(abs(y - z)) + mean(z) - 2 * sum((0:99) * sort(z)) / (100 * (99)),
    compute_pwm_crps(
      obs = y, pred = z
    )
  )
})
