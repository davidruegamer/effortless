setMethod("nrow", c("bdMatrix"),
          function(x) sum( sapply(x@listOfBlocks, NROW)))

setMethod("length", signature("bdMatrix"), 
          function(x) length(x@listOfBlocks))

setMethod("ncol", c("bdMatrix"),
          function(x) sum( sapply(x@listOfBlocks, NCOL)))

setMethod("dim", c("bdMatrix"),
          function(x) c(nrow(x),ncol(x)))

setMethod("abs", c("bdMatrix"),
          function(x) bdMatrix(lapply(x@listOfBlocks, abs)))

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
            
            x@listOfBlocks[[i]]
            
          }
)

setMethod("*",
          signature(e1="bdMatrix", e2="numeric"),
          function(e1, e2){

            e1@listOfBlocks <- mclapply(e1@listOfBlocks, function(l)l*e2)

            return(e1)
          }
)

setMethod("crossprod", signature(x="bdMatrix"),
          function(x) {
            
            bdMatrix(mclapply(x@listOfBlocks, function(y) crossprod(y)))
            
          } )

setMethod("%*%", signature(x = "bdMatrix", y = "bdMatrix"),
          function(x, y)
          {
            
            if( length(unique(c(sapply(x@listOfBlocks, NCOL), sapply(y@listOfBlocks, NROW)))) != 1 )
              stop("Multiplication only implemented if all blocks have matching sizes.")
            bdMatrix(mclapply(1:length(x), function(i) x[[i]] %*% y[[i]]))
            
          })

setMethod("%*%", signature(x = "bdMatrix", y = "numeric"),
          function(x, y)
          {
            
            if(ncol(x)!=length(y))
              stop("Dimension mismatch.")
            
            blocks <- sapply(x@listOfBlocks, NROW)

            sta <- c(1, cumsum(blocks)[-length(blocks)] + 1)
            end <- cumsum(blocks)
            
            unlist(mclapply(1:length(x), function(i) x[[i]] %*% y[sta[i]:end[i]]))
            
          })

setMethod("chol", signature(x="bdMatrix"),
          function(x) {
            
            bdMatrix(mclapply(x@listOfBlocks, function(y) chol(y)))
            
          } )

#' @export
setGeneric("svd")
# setGeneric("svd", def = function(x, nu, nv) svd(x, nu, nv),
#            signature(x = "matrix", nu = "numeric", nv = "numeric"))

#' @export
setMethod("svd", signature(x="bdMatrix"), 
          function(x, nu = min(nrow(x), p = ncol(x)), nv = min(nrow(x), p = ncol(x))) 
            {
  
  res <- mclapply(x@listOfBlocks, function(y) svd(y, nu = nu, nv = nv))
  list(d = unlist(lapply(res,"[[","d")),
       u = if(nu!=0) bdiag(lapply(res,"[[","u")) else NULL,
       v = if(nv!=0) bdiag(lapply(res,"[[","v")) else NULL)
  
          }
)

# setMethod("svd", signature(x="bdMatrix", nu="numeric", nv="numeric"), svd.bdMatrix)
# 
# svd.bdMatrix <- function(x) {
#   
#   res <- mclapply(x@listOfBlocks, function(y) svd(y))
#   list(d = unlist(lapply(res,"[[","d")),
#        u = bdiag(lapply(res,"[[","u")),
#        v = bdiag(lapply(res,"[[","v")))
#   
# }

# setMethod("svd", signature(x="bdMatrix"), svd.bdMatrix)

setMethod("forceSymmetric", c("bdMatrix"),
          function(x)
          {
            
            bdMatrix(mclapply(x, forceSymmetric))
            
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
              res <- bdMatrix(mclapply(1:length(e1), function(k) e1[[k]] + e2[indsplits$rows[[k]],
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
#' @export
setMethod("solve", c("bdMatrix"), 
          function(a) {
            
            if(any(sapply(a, function(xx) diff(dim(xx))!=0))) stop("solve only implemented for quadratic blocks.")
            
            res <- suppressWarnings(try(mclapply(a@listOfBlocks, solve)))
            if(length(cres <- which(sapply(res,class)=="try-error"))>0) stop(res[cres]) else bdMatrix(res)
            
            
          }
)

#' solve method for class bdMatrix
#' 
#' @param a object of class \code{bdMatrix}
#' @param b object of class \code{bdMatrix}
#' @param ... further arguments passed to \code{solve} (see \code{?solve})
#' 
#' @details Due to the block diagonal structure of \code{a}, solving can be performed separately on block levels
#' if \code{a} only consists of quadratic blocks. However, \code{solve} does not exploit the second argument \code{b}
#' as it is done in solve with vectors or matrices. Instead, \code{solve(a,b)} simply calculates \code{solve(a)\%*\%b}
#' blockwise for corresponding blocks.
#' @import parallel
#' @export
setMethod("solve", c("bdMatrix", "bdMatrix"), 
          function(a, b, ...) {
            
            res <- try(bdMatrix(mclapply(1:length(a), function(i) solve(a@listOfBlocks[[i]], 
                                                             as.vector(b@listOfBlocks[[i]]), ...))))
            if(class(res)=="try-error") return(solve(a)%*%b) else return(res)
            
          }
)

#' transpose method for bdMatrix objects
#' 
#' @param x object of class \code{bdMatrix}
#' 
#' @import parallel
#' @export
setMethod("t", c("bdMatrix"), 
          function(x) {
            
            bdMatrix(mclapply(x@listOfBlocks, t))
            
          }
)

setMethod("max", c("bdMatrix"),
          function(x) 
          {
            
            max(sapply(x@listOfBlocks, max))
            
          })