testthat::test_that("branching structure correctly refactored",{
    tst_N <- 200
    tst_id <- as.character(1:tst_N)
    tst_B <- sapply(
      1:tst_N, function(i){
        u <- runif(1)
        if (u < 0.25) {
          b <- 0
        } else {
          b <- sample(1:(i-1), 1)
        }
        b
      })
    tst_parent_id <- rep(NA, tst_N)
    tst_parent_id[tst_B == 0] <- "imm"
    tst_parent_id[tst_B != 0] <- tst_id[tst_B[tst_B != 0]]

    testthat::expect_equal(
      refactor_branching_structure(
        id = tst_id, parent_id = tst_parent_id, is_immigrant = (tst_B == 0), S = 10
      ),
      tst_B
    )

    tst_parent_id[tst_B == 0] <- NA
    testthat::expect_equal(
      refactor_branching_structure(
        id = tst_id, parent_id = tst_parent_id, is_immigrant = (tst_B == 0), S = 10
      ),
      tst_B
    )

    tst_id <- NA
    testthat::expect_error(
      refactor_branching_structure(
        id = tst_id, parent_id = tst_parent_id, is_immigrant = (tst_B == 0), S = 10
      )
    )
  })

testthat::test_that("trees labelled correctly",{
  tst_N <- 200
  tst_B <- sapply(
    1:tst_N, function(i){
      u <- runif(1)
      if (u < 0.25) {
        b <- 0
      } else {
        b <- sample(1:(i-1), 1)
      }
      b
    })
  tst_trees <- rep(NA, tst_N)
  k <- 1
  for (i in 1:tst_N) {
    if(tst_B[i] == 0) {
      tst_trees[i] <- k
      k <- k+1
    } else {
      tst_trees[i] <- tst_trees[tst_B[i]]
    }
  }

  expect_equal(
    tst_trees,
    identify_clusters(tst_B)
  )
})

