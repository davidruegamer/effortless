% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R, R/methods_bdMatrix.R
\docType{class}
\name{bdMatrix-class}
\alias{*,bdMatrix,numeric-method}
\alias{+,bdMatrix,dgCMatrix-method}
\alias{[,bdMatrix,ANY,ANY,ANY-method}
\alias{[[,bdMatrix-method}
\alias{\%*\%,bdMatrix,bdMatrix-method}
\alias{\%*\%,bdMatrix,numeric-method}
\alias{abs,bdMatrix-method}
\alias{bdMatrix-class}
\alias{chol,bdMatrix-method}
\alias{crossprod,bdMatrix,ANY-method}
\alias{dim,bdMatrix-method}
\alias{forceSymmetric,bdMatrix,ANY-method}
\alias{length,bdMatrix-method}
\alias{max,bdMatrix-method}
\alias{ncol,bdMatrix-method}
\alias{nrow,bdMatrix-method}
\title{Matrix class and methods for block diagonal matrices}
\usage{
\S4method{nrow}{bdMatrix}(x)

\S4method{length}{bdMatrix}(x)

\S4method{ncol}{bdMatrix}(x)

\S4method{dim}{bdMatrix}(x)

\S4method{abs}{bdMatrix}(x)

\S4method{[}{bdMatrix,ANY,ANY,ANY}(x, i, j, ..., drop = FALSE)

\S4method{[[}{bdMatrix}(x, i, j, ...)

\S4method{*}{bdMatrix,numeric}(e1, e2)

\S4method{crossprod}{bdMatrix,ANY}(x)

\S4method{\%*\%}{bdMatrix,bdMatrix}(x, y)

\S4method{\%*\%}{bdMatrix,numeric}(x, y)

\S4method{chol}{bdMatrix}(x)

\S4method{forceSymmetric}{bdMatrix,ANY}(x)

\S4method{+}{bdMatrix,dgCMatrix}(e1, e2)

\S4method{max}{bdMatrix}(x)
}
\arguments{
\item{x}{object of class \code{bdMatrix}}

\item{i}{integer for subsetting rows of \code{x}}

\item{j}{integer for subsetting columns of \code{x}}

\item{...}{further arguments passed to ordinary subsetting functions}

\item{drop}{logical. If TRUE the result is coerced to the lowest possible dimension.}

\item{e1}{see \code{x}}

\item{e2}{see \code{x}}

\item{y}{see \code{x}}
}
\description{
Matrix class and methods for block diagonal matrices
}

