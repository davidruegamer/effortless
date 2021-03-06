% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_bdMatrix.R
\docType{methods}
\name{solve,bdMatrix,ANY-method}
\alias{solve,bdMatrix,ANY-method}
\title{solve method for class bdMatrix}
\usage{
\S4method{solve}{bdMatrix,ANY}(a)
}
\arguments{
\item{a}{object of class \code{bdMatrix}}
}
\description{
solve method for class bdMatrix
}
\details{
Due to the block diagonal structure of \code{a}, solving can be performed separately on block levels
if \code{a} only consists of quadratic blocks.
}

