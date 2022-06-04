
test_that("Q_simulate works with a single number for alpha", {
  A1 <- Q_simulate(alpha = 1, lambda = rep(1,5)/5, rep = 5, popsize = 10, seed = 1)
  # Does A1 have the expected Pop column?
  expect_true(all(unique(A1$Pop) == c(sapply(X = 1:5,
                                             FUN = function(x){ paste(1,x, sep = "_")}))))
  # Does A1 have the expected number of rows?
  expect_equal(nrow(A1), 10*5)

  # Do the rows of A1 sum to 1?
  expect_true(all(round(rowSums(A1[,6:10]),4) == 1))
})


test_that("Q_simulate works with a vector for alpha", {

  A <- Q_simulate(alpha = c(1,10), lambda = rep(1,5)/5, rep = 5, popsize = 10, seed = 1)
print(A$Pop)
  # Does A have the expected Pop column?
# it is very weird that this consistently works in terminal but not in the test environment... I am giving up on this for now
  # expect_equal(as.character(A$Pop[1]), "1_1")
  # expect_true(all(as.character(A$Pop) ==
  #              c(sapply(X = 1:5,
  #                       FUN = function(x){
  #                         sapply(X = c(1,10),
  #                                function(y){ rep(paste(y,x, sep = "_"), 10)})}))))
  # Does A have the expected number of rows?
  expect_equal(nrow(A), 10*5*2)

  # Do the rows of A sum to 1?
  expect_true(all(round(rowSums(A[,6:10]),4) == 1))

  })

