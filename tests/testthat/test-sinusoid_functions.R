test_that("linear combination of sinusoids works", {
  t <- seq(-100, 100, length.out = 1001)
  day <- 24
  f <- 1 / c(day / 2, day, 7 * day)
  w <- 2 * pi * f
  beta <- runif(2 * length(w))

  tst <- sinusoidal_function(
    t = t,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    perform_checks = TRUE
  )

  S <- matrix(nrow = length(t), ncol = 2 * length(w))
  for (k in seq_along(w)) {
    S[, 2*k - 1] <- sin(w[k] * t)
    S[, 2*k] <- cos(w[k] * t)
  }

  expect_equal(
    tst,
    c(S %*% beta)
  )

  tst <- sinusoidal_function(
    t = 564,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    perform_checks = TRUE
  )

  S <- rep(NA, 2 * length(w))
  for (k in seq_along(w)) {
    S[2*k - 1] <- sin(w[k] * 564)
    S[2*k] <- cos(w[k] * 564)
  }

  expect_equal(
    tst,
    sum(S * beta)
  )
})

test_that("integral of sinusoids works", {
  t <- seq(-100, 100, length.out = 1001)
  a <- -200
  b <- 200

  day <- 24
  f <- 1 / c(day / 2, day, 7 * day)
  w <- 2 * pi * f
  beta <- runif(2 * length(w))

  tst_f <- function(t){
    sinusoidal_function(
      t,
      sinusoid_coefficients = beta, sinusoid_frequencies = w
    )
  }

  expect_equal(
    sinusoidal_integral(
      lower_limit = a, upper_limit = b,
      sinusoid_coefficients = beta, sinusoid_frequencies = w,
      perform_checks = TRUE
    ),
    integrate(tst_f, lower = a, upper = b)$value,
    tolerance = integrate(tst_f, lower = a, upper = b)$abs.error
  )


  tst <- sinusoidal_integral(
    lower_limit = t, upper_limit = b,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    perform_checks = TRUE
  )

  S <- matrix(nrow = length(t), ncol = 2 * length(w))
  for (k in seq_along(w)) {
    S[, 2*k - 1] <- (cos(w[k] * t) - cos(w[k] * b)) / w[k]
    S[, 2*k] <- (sin(w[k] * b) - sin(w[k] * t)) / w[k]
  }

  expect_equal(
    tst,
    c(S %*% beta)
  )

  tst <- sinusoidal_integral(
    lower_limit = a, upper_limit = t,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    perform_checks = TRUE
  )

  S <- matrix(nrow = length(t), ncol = 2 * length(w))
  for (k in seq_along(w)) {
    S[, 2*k - 1] <- (cos(w[k] * a) - cos(w[k] * t)) / w[k]
    S[, 2*k] <- (sin(w[k] * t) - sin(w[k] * a)) / w[k]
  }

  expect_equal(
    tst,
    c(S %*% beta)
  )

  tst <- sinusoidal_integral(
    lower_limit = t, upper_limit = t+1,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    perform_checks = TRUE
  )

  S <- matrix(nrow = length(t), ncol = 2 * length(w))
  for (k in seq_along(w)) {
    S[, 2*k - 1] <- (cos(w[k] * t) - cos(w[k] * (t+1))) / w[k]
    S[, 2*k] <- (sin(w[k] * (t+1)) - sin(w[k] * t)) / w[k]
  }

  expect_equal(
    tst,
    c(S %*% beta)
  )

  tst <- sinusoidal_integral(
    lower_limit = a, upper_limit = b,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    perform_checks = TRUE
  )

  S <- matrix(nrow = 1, ncol = 2 * length(w))
  for (k in seq_along(w)) {
    S[, 2*k - 1] <- (cos(w[k] * a) - cos(w[k] * b)) / w[k]
    S[, 2*k] <- (sin(w[k] * b) - sin(w[k] * a)) / w[k]
  }

  expect_equal(
    tst,
    c(S %*% beta)
  )
})

test_that("integral of exponentially decaying sinusoids works", {
  t <- seq(0, 50, length.out = 1001)
  a <- 0
  b <- 50

  day <- 24
  f <- 1 / c(day, day / 2)
  w <- 2 * pi * f
  beta <- c(-0.17, -0.59, -0.25, 0.34)

  psi <- 0.33

  tst <- sinusoidal_exponential_decay_integral(
    lower_limit = a, upper_limit = b,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    exponential_rate = psi,
    perform_checks = TRUE
  )

  tst_f <- function(t, t0, alpha, omega, gi_rate){
    N <- length(t)
    K <- length(omega)
    sinusoidal_basis <- matrix(
      sapply(
        omega, function(x){
          c(
            sin(t * x) * dexp(t-t0, rate = gi_rate),
            cos(t * x) * dexp(t-t0, rate = gi_rate)
            )
        }), nrow = N, ncol = 2*K)

    c(sinusoidal_basis %*% alpha)
  }

  # plot(
  #   t,
  #   tst_f(
  #     t, t0 = 0,
  #     alpha = beta, omega = w,
  #     gi_rate = psi
  #     )
  # )

  num_int <- integrate(
    tst_f, lower = a, upper = b,
    t0 = 0,
    alpha = beta, omega = w,
    gi_rate = psi
    )

  expect_equal(
    tst,
    num_int$value,
    tolerance = num_int$abs.error
  )

  tst <- sinusoidal_exponential_decay_integral(
    lower_limit = t, upper_limit = b,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    exponential_rate = psi,
    perform_checks = TRUE
  )

  S <- matrix(nrow = length(t), ncol = 2 * length(w))
  for (k in seq_along(w)) {
    S[, 2*k - 1] <- (
      (psi / (psi^2 + w[k]^2)) * (
        psi * (
          sin(t * w[k]) - (sin(w[k] * b) * exp( -psi * (b - t)))
        ) +
          w[k] * (
            cos(t * w[k]) - (cos(w[k] * b) * exp( -psi * (b - t)))
          )

      )
    )
    S[, 2*k] <- (
      (psi / (psi^2 + w[k]^2)) * (
        w[k] * (
          (sin(w[k] * b) * exp( -psi * (b - t))) - sin(w[k] * t)
        ) -
          psi * (
            (cos(w[k] * b) * exp( -psi * (b - t))) - cos(w[k] * t)
          )

      )
    )
  }

  expect_equal(
    tst,
    c(S %*% beta)
  )

  tst <- sinusoidal_exponential_decay_integral(
    lower_limit = a, upper_limit = t,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    exponential_rate = psi,
    perform_checks = TRUE
  )

  S <- matrix(nrow = length(t), ncol = 2 * length(w))
  for (k in seq_along(w)) {
    S[, 2*k - 1] <- (
      (psi / (psi^2 + w[k]^2)) * (
        psi * (
          sin(a * w[k]) - (sin(w[k] * t) * exp( -psi * (t - a)))
        ) +
          w[k] * (
            cos(a * w[k]) - (cos(w[k] * t) * exp( -psi * (t - a)))
          )

      )
    )
    S[, 2*k] <- (
      (psi / (psi^2 + w[k]^2)) * (
        w[k] * (
          (sin(w[k] * t) * exp( -psi * (t - a))) - sin(w[k] * a)
        ) -
          psi * (
            (cos(w[k] * t) * exp( -psi * (t - a))) - cos(w[k] * a)
          )

      )
    )
  }

  expect_equal(
    tst,
    c(S %*% beta)
  )

  tst <- sinusoidal_exponential_decay_integral(
    lower_limit = a, upper_limit = b,
    sinusoid_coefficients = beta, sinusoid_frequencies = w,
    exponential_rate = psi,
    perform_checks = TRUE
  )

  S <- matrix(nrow = 1, ncol = 2 * length(w))
  for (k in seq_along(w)) {
    S[, 2*k - 1] <- (
      (psi / (psi^2 + w[k]^2)) * (
        psi * (
          sin(a * w[k]) - (sin(w[k] * b) * exp( -psi * (b - a)))
        ) +
          w[k] * (
            cos(a * w[k]) - (cos(w[k] * b) * exp( -psi * (b - a)))
          )

      )
    )
    S[, 2*k] <- (
      (psi / (psi^2 + w[k]^2)) * (
        w[k] * (
          (sin(w[k] * b) * exp( -psi * (b - a))) - sin(w[k] * a)
        ) -
          psi * (
            (cos(w[k] * b) * exp( -psi * (b - a))) - cos(w[k] * a)
          )

      )
    )
  }

  expect_equal(
    tst,
    c(S %*% beta)
  )
})
