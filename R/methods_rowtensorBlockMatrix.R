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
  t(KhatriRao(t(bdiag(x@matLeft@listOfBlocks)),
              t(x@matRight))
    ) else stop("Function only defined for rowtensorBlockMatrix class.")

#' crossprod method for rowtensorBlockMatrix objects
#' 
#' @param x object of class \code{rowtensorBlockMatrix}
#' @param y see \code{x}
#' 
#' @import parallel
setMethod("crossprod", signature(x="rowtensorBlockMatrix", y="rowtensorBlockMatrix"),
          function(x, y) {
            
            if(nrow(x)!=nrow(y)) stop("Mismatch of dimensions.")
            if(length(x@matLeft)!=length(y@matLeft))
              stop("crossprod only implemented for equal number of blocks in x and y.")
            
            blocks <- sapply(x@matLeft@listOfBlocks, NROW)
            
            if(length(unique(blocks - 
                             sapply(y@matLeft@listOfBlocks, NROW))) !=1 )
              stop("Block dimensions do not match.")
            
            sta <- c(1, cumsum(blocks)[-length(blocks)] + 1)
            end <- cumsum(blocks)
            
            bdMatrix(mclapply(1:length(x@matLeft), function(i)
              crossprod(x@matLeft[[i]] * x@matRight[sta[i]:end[i], , drop = FALSE],
                        y@matLeft[[i]] * y@matRight[sta[i]:end[i], , drop = FALSE])
            ))
            
          } )

#' crossprod method for rowtensorBlockMatrix object
#' 
#' @param x object of class \code{rowtensorBlockMatrix}
#' 
#' @import parallel
setMethod("crossprod", signature(x="rowtensorBlockMatrix"),
          function(x) {
            
            crossprod(x = x, y = x)
            
          } 
)

setMethod("crossprod", signature(x="rowtensorBlockMatrix","numeric"),
          function(x, y) {
            
            blocks <- sapply(x@matLeft@listOfBlocks, NROW)
            sta <- c(1, cumsum(blocks)[-length(blocks)] + 1)
            end <- cumsum(blocks)
            
            bdMatrix(mclapply(1:length(x@matLeft), function(i)
              crossprod(x@matRight[sta[i]:end[i], , drop = FALSE], x@matLeft[[i]] * y[sta[i]:end[i]])
              ))
            
          } 
)

setMethod("%*%", signature(x="rowtensorBlockMatrix","numeric"),
          function(x, y) {
            
            blocks <- sapply(x@matLeft@listOfBlocks, NROW)
            sta <- c(1, cumsum(blocks)[-length(blocks)] + 1)
            end <- cumsum(blocks)
            
            blocksY <- cumsum(rep(ncol(x@matRight), length(x@matLeft)))
            staY <- c(1,blocksY[-length(blocksY)] + 1)
            endY <- blocksY
            
            res <- bdMatrix(mclapply(1:length(x@matLeft), function(i)
              x@matLeft[[i]] * x@matRight[sta[i]:end[i], , drop = FALSE] %*% y[staY[i]:endY[i]])
            )
            
            if(length(unique(c(1, sapply(res@listOfBlocks, ncol))))==1)
              return(as(res, "vector")) else return(res)
            
          } 
)

setMethod("*", signature(e1="rowtensorBlockMatrix", e2="numeric"),
          function(e1, e2){
            
            blocks <- sapply(e1@matLeft@listOfBlocks, NROW)
            sta <- c(1, cumsum(blocks)[-length(blocks)] + 1)
            end <- cumsum(blocks)
            
            for(i in 1:length(e1@matLeft)) e1@matLeft[[i]] <- e1@matLeft[[i]]*e2[sta[i]:end[i]]
            
            return(e1)
          }
)

setMethod("dim", c("rowtensorBlockMatrix"),
          function(x) c(nrow(x@matRight),length(x@matLeft)*ncol(x@matRight)))

setMethod("abs", signature("rowtensorBlockMatrix"),
          function(x) 
            {
            
            x@matLeft <- abs(x@matLeft)
            x@matRight <- abs(x@matRight)
            return(x)
            
          })


setMethod("max", c("rowtensorBlockMatrix"),
          function(x) 
          {
            
            m1 <- max(x@matLeft)
            m2 <- max(x@matRight)
            return(m1*m2)
            
          })

setGeneric("rankMatrix")

#' rankMatrix function for rowtensorBlockMatrix objects
#' 
#' @export
rankMatrix.rowtensorBlockMatrix <- function(x, tol = NULL,
                                            method = c("tolNorm2", "qr.R", "qrLINPACK", "qr",
                                                       "useGrad", "maybeGrad"),
                                            sval,
                                            warn.t = TRUE) 
{
  
  message("Ignoring sval argument in rankMatrix call.")
  ncol(x@matLeft)*rankMatrix(x = x@matRight, tol = tol, method = method, warn.t = warn.t)
  
} 

setMethod("rankMatrix", signature(x = "rowtensorBlockMatrix", tol = "numeric", 
                                  method = "character", sval = "function", warn.t = "logical"
), 
rankMatrix.rowtensorBlockMatrix)

# rankMatrix.rowtensorBlockMatrix <- function(x,
#                                             method = c("tolNorm2", "qr.R", "qrLINPACK", "qr",
#                                                        "useGrad", "maybeGrad"),
#                                             warn.t = TRUE) 
# {
#   
#   message("Ignoring sval argument in rankMatrix call.")
#   ncol(x@matLeft)*rankMatrix(x = x@matRight, method = method, warn.t = warn.t)
#   
# } 
# 
# setMethod("rankMatrix", signature(x = "rowtensorBlockMatrix"#, tol = "numeric", 
#                                   #method = "character", sval = "function", warn.t = "logical"
#                                   ), 
#           rankMatrix.rowtensorBlockMatrix)

# svd.rowtensorBlockMatrix <- function(x, nu = min(nrow(x), p = ncol(x)), nv = min(nrow(x), p = ncol(x))) {
#   
#   res <- mclapply(1:ncol(x@matLeft), function(i) 
#     svd(x@matLeft[[i]] * x@matRight, nu = nu, nv = nv))
#   list(d = unlist(lapply(res,"[[","d")),
#        u = if(nu!=0) bdiag(lapply(res,"[[","u")) else NULL,
#        v = if(nv!=0) bdiag(lapply(res,"[[","v")) else NULL)
#   
# }
# 
# setMethod("svd", signature(x="rowtensorBlockMatrix"), svd.rowtensorBlockMatrix)