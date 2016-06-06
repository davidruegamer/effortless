#' Matrix class and methods for block diagonal matrices
#' 
#' @param x object of class \code{bdMatrix}
#' @param i integer for subsetting rows of \code{x}
#' @param j integer for subsetting columns of \code{x}
#' @param ... further arguments passed to ordinary subsetting functions
#' @param drop logical. If TRUE the result is coerced to the lowest possible dimension.
#' @param e1 see \code{x}
#' @param e2 see \code{x}
#' @param y see \code{x}
#' 
#' @import Matrix
#' @import methods
#' @importFrom methods as new
#' @rdname bdMatrix-class
#' 
#' @export
setClass("bdMatrix",
         slots = list(listOfBlocks = "list"),
         prototype = prototype(listOfBlocks = list(matrix(1:4, 2), diag(3))),
         validity = function(object) {
           
           invalids <- character(0)
           
           if(!is.list(object@listOfBlocks)) 
             invalids <- "listOfBlocks must be a list."
           
           if(!all(sapply(object@listOfBlocks, function(x) is.numeric(x) | is(x, "Matrix")))) 
             invalids <- c(invalids, "listOfBlocks must only contain matrices and vectors.")
           
           if(length(invalids)) invalids else TRUE
         },
         contains = "list"
)

#' Constructor for block diagonal matrices
#' 
#' @param listOfBlocks a list of matrices or numeric vectors, which represent the blocks in the given order
#'
#' @import Matrix
#' 
#' @examples
#' ## construct a bdMatrix
#' bdMatrix(listOfBlocks = list(matrix(1:4, 2), diag(3)))
#' ## see what the whole matrix looks like by using Matrix::bdiag of the list elements
#' bdiag(bdMatrix(listOfBlocks = list(matrix(1:4, 2), diag(3)))@listOfBlocks)
#' 
#' @export
bdMatrix <- function(listOfBlocks)
{
  
  new("bdMatrix", listOfBlocks = listOfBlocks)
  
}

#' Matrix class and methods for block structured kronecker sums
#' 
#' @param x object of class \code{kroneckersumBlockMatrix}
#' @param e1 see \code{x}
#' @param e2 see \code{x}
#' 
#' @rdname kroneckersumBlockMatrix-class
#' @import Matrix
#' @export
setClass("kroneckersumBlockMatrix",
         slots = list(matLeft = "ddiMatrix", matRight = "Matrix"),
         prototype = prototype(matLeft = Diagonal(5), 
                               matRight = Matrix(1:9, ncol=3)),
         validity = function(object) {
           
           invalids <- character(0)
           
           if(!inherits(object@matRight, "Matrix")) 
             invalids <- "matRight must be a class inheriting from 'Matrix'."
           
           if(!(diff(dim(object@matLeft))==0 & diff(dim(object@matRight))==0))
             invalids <- c(invalids, "Matrices must be quadratic.")
           
           # if(all(sapply(object@blockInd, function(x) x %in% (1:nrow(object@matLeft))))) 
           #   invalids <- c(invalids, 
           #                 "blockInd must be numeric values in the range of 1:nrow(matLeft)")
           
           if(length(invalids)) invalids else TRUE
         }
)

#' Constructor for kroneckersumBlockMatrix objects
#' 
#' @param X1 a diagonal matrix of class \code{ddiMatrix}
#' @param X2 a \code{Matrix} object
#'
#' @import Matrix
#' 
#' @examples 
#' matLeft = Diagonal(5)
#' matRight = Matrix(1:9, ncol=3)
#' kroneckersumBlockMatrix(matLeft, matRight)
#' 
#' @export
kroneckersumBlockMatrix <- function(X1, X2)
{
  
  new("kroneckersumBlockMatrix", matLeft=X1, matRight=X2)
  
}

#' Matrix class and methods for block structured row-wise tensor products
#' 
#' @param x object of class \code{rowtensorBlockMatrix}
#' @param i integer for subsetting rows of \code{x}
#' @param j integer for subsetting columns of \code{x}
#' @param ... further arguments passed to ordinary subsetting functions
#' @param drop logical. If TRUE the result is coerced to the lowest possible dimension.
#' @param e1 see \code{x}
#' @param e2 see \code{x}
#' @param y see \code{x}
#' 
#' @rdname rowtensorBlockMatrix-class
#' @import Matrix
#' @export
setClass("rowtensorBlockMatrix",
         slots = list(matLeft = "bdMatrix", matRight = "Matrix"),
         prototype = prototype(matLeft = bdMatrix(list(c(1), matrix(c(1,1,1), ncol=1))), 
                               matRight = Matrix(1:8, ncol=2)),
         validity = function(object) {
           
           invalids <- character(0)
           
           if(!inherits(object@matRight, "Matrix")) 
             invalids <- "matRight must be a class inheriting from 'Matrix'."
           
           if(nrow(object@matLeft)!=NROW(object@matRight)) 
             invalids <- c(invalids, "Matrices must have equal number of rows.")
           
           if(any(sapply(object@matLeft, NCOL) != 1)) 
             stop("An object of class rowtensorBlockMatrix can not have blocks with more than one column.")
           
           if(any(unlist(object@matLeft)==0)) stop("Zero entries in blocks are not allowed.")
           
           # if(all(sapply(object@blockInd, function(x) x %in% (1:nrow(object@matLeft))))) 
           #   invalids <- c(invalids, 
           #                 "blockInd must be numeric values in the range of 1:nrow(matLeft)")
           
           if(length(invalids)) invalids else TRUE
         }
)

#' Constructor for rowtensorBlockMatrix objects
#'
#' @param X1 a \code{bdMatrix} object
#' @param X2 a \code{Matrix} object
#' 
#' @examples 
#' matLeft = bdMatrix(list(c(1), matrix(c(1,1,1), ncol=1)))
#' matRight = Matrix(1:8, ncol=2)
#' rowtensorBlockMatrix(matLeft, matRight)
#'
#' @import Matrix
#' @export
rowtensorBlockMatrix <- function(X1, X2)
{
  
  new("rowtensorBlockMatrix", matLeft=X1, matRight=X2)
  
}
