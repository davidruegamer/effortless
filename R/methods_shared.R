#' crossprod function for crossproduct of a bdMatrix and a kroneckersumBlockMatrix
#' 
#' @param x object of class \code{bdMatrix}
#' @param y object of class \code{kroneckersumBlockMatrix}
#' 
#' @details if \code{x} and \code{y} have the same block sizes, the crossprod function
#' calculates the crossproduct of each block separately
#' 
setMethod("crossprod", signature(x="bdMatrix", y="kroneckersumBlockMatrix"),
          function(x, y) {
            
            stopifnot(ncol(x) == nrow(y))
            
            # stop if block sizes are not the same
            if(any(dim(y@matRight)!=dim(x[[1]])) | !hasEqualBlocksizes(x))
              stop("crossprod for bdMatrix and kroneckersumBlockMatrix only possible if block sizes are equal.")
            
            diagPart <- split(rep(diag(y@matLeft), each=ncol(y@matRight)), 
                              rep(1:ncol(y@matLeft), each=ncol(y@matRight)))
            bdMatrix(lapply(1:length(x), function(i) crossprod(x[[i]], 
                                                               (diag(diagPart[[i]]) + 
                                                                  y@matRight))))
            
          } )

#' Function to sum a bdMatrix and a kroneckersumBlockMatrix
#' 
#' @param e1 object of class \code{bdMatrix}
#' @param e2 object of class \code{kroneckersumBlockMatrix}
#' 
#' @details \code{x} and \code{y} must have the same number of blocks
#' WARNING: Sizes of blocks are not checked. This can lead to unwanted results.
#' 
setMethod("+", signature(e1="bdMatrix", e2="kroneckersumBlockMatrix"),
          function(e1, e2) {
            
            stopifnot(all(dim(e2@matLeft)==length(e1)))
            
            bdMatrix(lapply(1:length(e1), function(k) e1[[k]] + 
                              diag(rep(e2@matLeft[k,k], each=ncol(e2@matRight))) + 
                              e2@matRight
            ))
            
          } )


#' Matrix multiplication function for a rowtensorBlockMatrix and a bdMatrix
#' 
#' @param x object of class \code{bdMatrix}
#' @param y object of class \code{kroneckersumBlockMatrix}
#' 
#' @details if \code{x} and \code{y} have the same block sizes, the function
#' calculates the ordinary matrix product of each block separately
#' 
setMethod("%*%", signature(x="rowtensorBlockMatrix", y="bdMatrix"),
          function(x, y) {
            
            if(length(unique(c(ncol(x@matRight), sapply(y@listOfBlocks, NROW))))!=1)
              stop("Mismatch of blocks.")
            
            if(length(unique(c(1, sapply(y@listOfBlocks, NCOL))))==1)
              return(x %*% as(y, "vector"))
            
            blocks <- sapply(x@matLeft@listOfBlocks, NROW)
            sta <- c(1, cumsum(blocks)[-length(blocks)] + 1)
            end <- cumsum(blocks)
            
            bdMatrix(mclapply(1:length(x@matLeft), function(i)
              x@matLeft[[i]] * x@matRight[sta[i]:end[i], , drop = FALSE] %*% y[[i]])
            )

          } )

