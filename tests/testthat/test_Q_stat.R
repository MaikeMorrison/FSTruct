
test_that("Q_stat gives three numbers",{

  set.seed(1)
  A <- matrix(rbeta(n = 100, shape1 = 1, shape2 = 1), ncol = 5)
  A <- A/rowSums(A)

  expect_true(all(sapply(Q_stat(A), is.numeric)))
})

test_that("Q_stat gives correct answer based on provided K",{

  set.seed(1)
  A <- matrix(rbeta(n = 100, shape1 = 1, shape2 = 1), ncol = 5)
  A <- A/rowSums(A)

  expect_true(all(unlist(Q_stat(A)) == unlist(Q_stat(A, K = 5))))
})

test_that("Q_stat ratio = Fst/FstMax",{
  set.seed(1)
  A <- matrix(rbeta(n = 100, shape1 = 1, shape2 = 1), ncol = 5)
  A <- A/rowSums(A)
  stat = Q_stat(A)

  expect_true(stat$ratio == stat$Fst / stat$FstMax)
})

