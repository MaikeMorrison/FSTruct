
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

test_that("Q_bootstrap works regardless of K provided", {

  set.seed(1)
  A <- matrix(rbeta(n = 100, shape1 = 1, shape2 = 1), ncol = 5)
  A <- A/rowSums(A)

  B <- matrix(rbeta(n = 100, shape1 = 1, shape2 = 1), ncol = 5)
  B <- B/rowSums(B)

  boot <- Q_bootstrap(matrices = list(A, B), n_replicates = 100, seed = 1)

  boot$plot_boxplot
})
