setMethod("[", c("rowtensorBlockMatrix"),
          function(x, i, j, ..., drop=TRUE)
          {
            if(!missing(j)){
              
              j1 <- rep(1:ncol(x@matLeft), each=ncol(x@matRight))[j]
              j2 <- rep(1:ncol(x@matRight), ncol(x@matLeft))[j]
              rowtensorBlockMatrix(x@matLeft[i = i, j=j1, ..., drop = drop], 
                                   x@matRight[i = i, j=j2, ..., drop = drop]) 
              
            }else
              rowtensorBlockMatrix(x@matLeft[i = i, ..., drop = drop], 
                                   x@matRight[i = i, ..., drop = drop])
          }
)

evalMat <- function(x) if(class(x)=="rowtensorBlockMatrix") 
  t(KhatriRao(t(bdiag(x@matLeft)),t(x@matRight))) else stop("Function only defined for rowtensorBlockMatrix class.")

#' crossprod method for rowtensorBlockMatrix objects
#' 
#' @param x object of class \code{rowtensorBlockMatrix}
#' @param y see \code{x}
#' 
#' @import parallel
setMethod("crossprod", signature(x="rowtensorBlockMatrix", y="rowtensorBlockMatrix"),
          function(x, y) {
            
            warning("Method crossprod for signature rowtensorBlockMatrix,rowtensorBlockMatrix does only work if x is a multiple of y.")
            
            bdMatrix(mclapply(1:ncol(x@matLeft), function(i)
              crossprod(x@matLeft[[i]] * x@matRight,
                        y@matLeft[[i]] * y@matRight)
            ))
            
          } )

#' crossprod method for rowtensorBlockMatrix object
#' 
#' @param x object of class \code{rowtensorBlockMatrix}
#' 
#' @import parallel
setMethod("crossprod", signature(x="rowtensorBlockMatrix"),
          function(x) {
            
            bdMatrix(mclapply(1:ncol(x@matLeft), function(i)
              crossprod(x@matLeft[x@matLeft[,i]!=0,i] * x@matRight[x@matLeft[,i]!=0,])
            ))
            
          } 
)



setMethod("*", signature(e1="rowtensorBlockMatrix", e2="numeric"),
          function(e1, e2){
            
            e1@matLeft <- e1@matLeft*e2
            
            return(e1)
          }
)

setMethod("dim", c("rowtensorBlockMatrix"),
          function(x) c(nrow(x@matRight),length(x@matLeft)*ncol(x@matRight)))

setGeneric("rankMatrix")

rankMatrix.rowtensorBlockMatrix <- function(x, tol = NULL,
                                            method = c("tolNorm2", "qr.R", "qrLINPACK", "qr",
                                                       "useGrad", "maybeGrad"),
                                            sval = svd(x, 0, 0)$d, warn.t = TRUE) 
{
  
  ncol(x@matLeft)*rankMatrix(x@matRight, tol = tol, method = method, sval = sval, warn.t = warn.t)
  
} 

setMethod("rankMatrix", signature(x="rowtensorBlockMatrix"), rankMatrix.rowtensorBlockMatrix)
