context("Test bdMatrix class")

test_that("Error throwing when failing to build a bdMatrix", {
  expect_error(bdMatrix(matrix(1:4)))
  expect_error(bdMatrix(c(1:5)))
  expect_error(bdMatrix(list(matrix(1:4),NA)))
  expect_error(bdMatrix(list(matrix(1:4),c())))
  # expect_error(bdMatrix(list(matrix(1:4),1:2))) -> works!
  expect_error(bdMatrix(list(matrix(1:4),list(matrix(1:4)))))
  # expect_error(bdMatrix(list(matrix(1:4),Matrix(1:4)))) -> works!
})

test_that("Correct representation of a bdMatrix", {
  expect_equal(Matrix::bdiag(bdMatrix(list(matrix(1:4)))@listOfBlocks), 
                   bdiag(list(Matrix(1:4))))
  expect_equal(Matrix::bdiag(bdMatrix(list(matrix(1:4), matrix(1:9, ncol=3)))@listOfBlocks), 
               bdiag(list(Matrix(1:4), matrix(1:9, ncol=3))))
  expect_equal(Matrix::bdiag(bdMatrix(list(matrix(1:4), 1:6))@listOfBlocks), 
               bdiag(list(Matrix(1:4), 1:6)))
  expect_s4_class(Matrix::bdiag(bdMatrix(list(matrix(1:4), 1:6))@listOfBlocks),"dgCMatrix")
  expect_s4_class(Matrix::bdiag(bdMatrix(list(matrix(1:100, ncol=10), 1:6))@listOfBlocks),"dgCMatrix")
  expect_s4_class(Matrix::bdiag(bdMatrix(list(0, 1:6))@listOfBlocks),"dgCMatrix")
})

test_that("Working methods for bdMatrix", {
  # ncol
  expect_equal(ncol(bdMatrix(list(matrix(1:4), 1:6))), 2)
  expect_equal(ncol(bdMatrix(list(matrix(1:4), matrix(1:6,nrow=1)))), 7)
  expect_equal(ncol(bdMatrix(list(matrix(1:4), matrix(1:6,ncol=3)))), 4)
  # nrow
  expect_equal(nrow(bdMatrix(list(matrix(1:4), matrix(1:6,nrow=1)))), 5)
  expect_equal(nrow(bdMatrix(list(matrix(1:4), 1:6))), 10)
  expect_equal(nrow(bdMatrix(list(matrix(1:4), matrix(1:6,ncol=3)))), 6)
  # length
  expect_equal(length(bdMatrix(list())), 0)
  expect_equal(length(bdMatrix(list(matrix(1:4)))), 1)
  expect_equal(length(bdMatrix(list(matrix(1:4), matrix(1:6,ncol=3)))), 2)
  # dim
  expect_error(dim(bdMatrix(list())))
  expect_equal(dim(bdMatrix(list(matrix(1:4)))), c(4,1))
  expect_equal(dim(bdMatrix(list(matrix(1:4), matrix(1:6,ncol=3)))), c(6,4))
  # abs
  expect_equal(abs(bdMatrix(list(matrix(-1*1:4), matrix(1:6,ncol=3)))),
               bdMatrix(list(matrix(1:4), matrix(1:6,ncol=3))))
  # max 
  expect_equal(max(bdMatrix(list(matrix(-1*1:4), matrix(1:6,ncol=3)))), 6)
  # t
  expect_equal(t(bdMatrix(list(matrix(-1*1:4), matrix(1:6,ncol=3)))),
               bdMatrix(list(t(matrix(-1*1:4)), t(matrix(1:6,ncol=3)))))
  # solve
  expect_error(solve(bdMatrix(list(matrix(-1*1:4, ncol=2), matrix(1:9,ncol=3)))))
  expect_equal(as.matrix(bdiag(solve(bdMatrix(list(matrix(c(1,0,0,1), ncol=2), 
                                                   matrix(5))))@listOfBlocks)),
               solve(as.matrix(bdiag(list(matrix(c(1,0,0,1), ncol=2), 
                                          matrix(5))))), check.attributes = FALSE)
  # solve a,b
  expect_equal(bdMatrix(split(solve(bdMatrix(list(matrix(c(1,0,0,1), ncol=2),
                                   matrix(5)))) %*% c(1:3), c(1,1,2))),
               solve(a = bdMatrix(list(matrix(c(1,0,0,1), ncol=2),
                                       matrix(5))), b = bdMatrix(list(1:2, 3))),
               check.attributes = FALSE)
  # + 
  # *
  # %*%
  # crossprod with bd
  # crossprod with numeric
  # svd
  # chol
  # [
  # [[
})

  
  
  