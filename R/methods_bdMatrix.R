setMethod("nrow", c("bdMatrix"),
          function(x) sum( sapply(x@listOfBlocks, NROW)))

setMethod("length", c("bdMatrix"),
          function(x) length(x@listOfBlocks))

setMethod("ncol", c("bdMatrix"),
          function(x) sum( sapply(x@listOfBlocks, NCOL)))

setMethod("dim", c("bdMatrix"),
          function(x) c(nrow(x),ncol(x)))

setMethod("[", c("bdMatrix"),
          function(x, i, j, ..., drop=FALSE)
          {
            
            indsplit <- getBlockIndices(x)
            
            iMiss <- missing(i) 
            jMiss <- missing(j)
            
            x <- lapply(1:length(x), 
                        function(k){
                          ii <- if(!iMiss) which(indsplit$rows[[k]] %in% i) else 1:length(indsplit$rows[[k]])
                          jj <- if(!jMiss) which(indsplit$cols[[k]] %in% j) else 1:length(indsplit$cols[[k]])
                          
                          if(length(ii)==0 | length(jj)==0) return(NULL) else
                            x[[k]][ii, jj, ..., drop=drop]
                        }
            )
            
            bdMatrix(x[!sapply(x, is.null)])
            
          }
)

setMethod("[[", c("bdMatrix"),
          function(x, i, ...)
          {
            
            bdMatrix(x@listOfBlocks[[i]])
            
          }
)

setMethod("*",
          signature(e1="bdMatrix", e2="numeric"),
          function(e1, e2){

            e1@listOfBlocks <- lapply(e1@listOfBlocks, function(l)l*e2)

            return(e1)
          }
)

setMethod("crossprod", signature(x="bdMatrix"),
          function(x) {
            
            bdMatrix(lapply(x@listOfBlocks, function(y) crossprod(y)))
            
          } )

setMethod("%*%", signature(x = "bdMatrix", y = "bdMatrix"),
          function(x, y)
          {
            
            if( length(unique(c(sapply(x, NCOL), sapply(y, NROW)))) != 1 )
              stop("Multiplication only implemented if all blocks have matching sizes.")
            bdMatrix(lapply(1:length(x), function(i) x[[i]] %*% y[[i]]))
            
          })

setMethod("chol", signature(x="bdMatrix"),
          function(x) {
            
            bdMatrix(lapply(x@listOfBlocks, function(y) chol(y)))
            
          } )

setMethod("svd", signature(x="bdMatrix"),
          function(x, nu = min(nrow(x), p = ncol(x)), nv = min(nrow(x), p = ncol(x))) {
            
            res <- mclapply(x@listOfBlocks, function(y) svd(y, nu, nv))
            list(d = unlist(sapply(res,"[[","d")),
                 u = bdiag(lapply(res,"[[","u")),
                 v = bdiag(lapply(res,"[[","v")))
            
          } )

setMethod("forceSymmetric", c("bdMatrix"),
          function(x)
          {
            
            bdMatrix(lapply(x, forceSymmetric))
            
          }
)

setMethod("+", signature(e1="bdMatrix", e2="dgCMatrix"),
          function(e1, e2) {
            
            stopifnot(all.equal(dim(e1),dim(e2))=="TRUE")
            
            minDim <- min(sapply(e1@listOfBlocks, function(l) min(dim(l))))
            if(any(as.numeric(triu(e2, k = minDim))!=0)){ 
              
              warning("dgCMatrix with non-zero values outside the block range of bdMatrix. 
                                  Coercing bdMatrix to dgCMatrix.")
              
              res <- bdiag(e1@listOfBlocks) + e2
              
            }else{
              
              indsplits <- getBlockIndices(e1)
              res <- bdMatrix(lapply(1:length(e1), function(k) e1[[k]] + e2[indsplits$rows[[k]],
                                                                            indsplits$cols[[k]]]))
              
            }
            
            return(res)
            
          } )

#' solve method for class bdMatrix
#' 
#' @param a object of class \code{bdMatrix}
#' @details Due to the block diagonal structure of \code{a}, solving can be performed separately on block levels
#' if \code{a} only consists of quadratic blocks.
#' @import parallel
setMethod("solve", c("bdMatrix"), 
          function(a) {
            
            if(any(sapply(a, function(xx) diff(dim(xx))!=0))) stop("solve only implemented for quadratic blocks.")
            
            bdMatrix(mclapply(a@listOfBlocks, solve))
            
          }
)

#' transpose method for bdMatrix objects
#' 
#' @param x object of class \code{bdMatrix}
#' 
#' @import parallel
setMethod("t", c("bdMatrix"), 
          function(x) {
            
            bdMatrix(mclapply(x@listOfBlocks, t))
            
          }
)