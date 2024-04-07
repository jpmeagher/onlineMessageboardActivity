test_that("poisson likelihood matches scaled base R function", {
  y <- train_df %>%
    dplyr::filter(parent_id == 0) %>%
    magrittr::use_series(t)

  N <- length(y)
  a <- y %>% max() %>% ceiling()

  expect_equal(
    periodic_poisson_point_process_likelihood(
      t = y,
      a = a,
      lambda_0 = N / a,
      omega = NULL, alpha = NULL,
      log = TRUE,
      perform_checks = TRUE
    ),
    dpois(N, lambda = N, log = TRUE) +
      lfactorial(N) -
      (N * log(a))
  )
})


test_that("poisson likelihood matches calculation", {
  day <- 24
  f <- 1 / c(day, day / 2)
  omega <- 2 * pi * f
  alpha <- c(-0.17, -0.59, -0.25, 0.34)

  y <- train_df %>%
    dplyr::filter(parent_id == 0) %>%
    magrittr::use_series(t)

  N <- length(y)
  a <- y %>% max() %>% ceiling()

  compensator <- (N / a) * (
    a + sinusoidal_integral(
      0, a, alpha, omega
    ))

  intensity <- (N / a) * (
    1 + sinusoidal_function(
      y, alpha, omega
    )
  )



  expect_equal(
    periodic_poisson_point_process_likelihood(
      t = y,
      a = a,
      lambda_0 = N / a,
      omega = omega,
      alpha = alpha,
      log = TRUE,
      perform_checks = TRUE
    ),
    sum(log(intensity)) - compensator
  )
})

