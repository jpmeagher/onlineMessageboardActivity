test_that("individual reproduction number posterior", {
  a <- 10
  arr <- sort(runif(10, -a, a))
  z <- rnbinom(10, size = 0.1, mu = 2)

  R <- 0.66
  k <- 0.1

  psi <- 0.25

  day <- 24
  f <- 1 / c(7 * day, day, day / 2)
  w <- 2 * pi * f
  beta <- runif(2 * length(w))

  post <- gamma_posterior_individual_reproduction_number(
    t = arr, observed_offspring = z,
    observation_horizon = a,
    reproduction_number = R, dispersion_parameter = k,
    gi_exp_decay_rate = psi,
    sinusoid_coefficients = beta, sinusoid_frequencies = w
  )

  expect_equal(
    k + z, post$shape
  )

  expect_equal(
    (k / R) + pexp(a - arr, rate = psi) +
      sinusoidal_exponential_decay_integral(
        lower_limit = arr, upper_limit = a,
        sinusoid_coefficients = beta, sinusoid_frequencies = w,
        exponential_rate = psi
      ),
    post$rate
  )

  post <- gamma_posterior_individual_reproduction_number(
    t = arr, observed_offspring = z,
    observation_horizon = a,
    reproduction_number = R, dispersion_parameter = k,
    gi_exp_decay_rate = psi
  )

  expect_equal(
    k + z, post$shape
  )

  expect_equal(
    (k / R) + pexp(a - arr, rate = psi) ,
    post$rate
  )
})
