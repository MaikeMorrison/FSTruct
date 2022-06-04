
test_that("Q_bootstrap works regardless of K provided", {

  set.seed(1)
  A <- matrix(rbeta(n = 100, shape1 = 1, shape2 = 1), ncol = 5)
  A <- A/rowSums(A)

  B <- matrix(rbeta(n = 100, shape1 = 1, shape2 = 1), ncol = 5)
  B <- B/rowSums(B)

  boot_noK <- Q_bootstrap(matrices = list(A, B), n_replicates = 100, seed = 1)
  boot_K <- Q_bootstrap(matrices = list(A, B), n_replicates = 100, K = 5, seed = 1)
  boot_Kvec <- Q_bootstrap(matrices = list(A, B), n_replicates = 100, K = c(5,5), seed = 1)

  expect_true(all(boot_noK$statistics == boot_K$statistics))
  expect_true(all(boot_noK$statistics == boot_Kvec$statistics))
  expect_true(all(boot_K$statistics == boot_Kvec$statistics))
})

test_that("Q_bootstrap works on a single matrix",{
  A <- Q_simulate(alpha = 1, lambda = rep(.2, 5),
                  rep = 1, popsize = 10, seed = 1)

  # specify K
  boot <- Q_bootstrap(matrices = A, n_replicates = 10, K = 5, seed = 1)
  expect_equal(boot$statistics$ratio[1], 0.51565959)

  # don't specify K
  boot <- Q_bootstrap(matrices = A[, 6:10], n_replicates = 10, seed = 1)
  expect_equal(boot$statistics$ratio[1], 0.51565959)
})

test_that("Q_bootstrap long matrix option works",{
  Q <- Q_simulate(alpha = 1, lambda = rep(.2, 5),
                  rep = 15, popsize = 10, seed = 1)
  boot <-  Q_bootstrap(matrices = Q, n_replicates = 100,
                       K = 5, seed = 1, group = "Pop")
  expect_equal(boot$statistics$ratio[1], 0.51565959)
})


