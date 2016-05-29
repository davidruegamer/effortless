evalMat <- function(x) if(class(x)=="kroneckersumBlockMatrix") 
  kronecker(x@matLeft, diag(ncol(x@matRight))) + kronecker(diag(ncol(x@matLeft)), x@matRight) else 
    stop("Function only defined for kroneckersumBlockMatrix class.")


setMethod("*", signature(e1="kroneckersumBlockMatrix", e2="numeric"),
          function(e1, e2){

            e1@matRight <- e1@matRight*e2
            e1@matLeft <- e1@matLeft*e2

            return(e1)
          }
)


setMethod("nrow", c("kroneckersumBlockMatrix"),
          function(x) nrow(x@matLeft) * nrow(x@matRight))