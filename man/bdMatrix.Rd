% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{bdMatrix}
\alias{bdMatrix}
\title{Constructor for block diagonal matrices}
\usage{
bdMatrix(listOfBlocks)
}
\arguments{
\item{listOfBlocks}{a list of matrices or numeric vectors, which represent the blocks in the given order}
}
\description{
Constructor for block diagonal matrices
}
\examples{
## construct a bdMatrix
bdMatrix(listOfBlocks = list(matrix(1:4, 2), diag(3)))
## see what the whole matrix looks like by using Matrix::bdiag of the list elements
bdiag(bdMatrix(listOfBlocks = list(matrix(1:4, 2), diag(3)))@listOfBlocks)

}

