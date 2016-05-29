setMethod("crossprod", signature(x="bdMatrix", y="kroneckersumBlockMatrix"),
          function(x, y) {
            
            stopifnot(ncol(x) == nrow(y))
            
            # stop if block sizes are not the same
            if(any(dim(y@matRight)!=dim(x[[1]])) | !hasEqualBlocksizes(x))
              stop("crossprod for bdMatrix and kroneckersumBlockMatrix only possible if block sizes are equal.")
            
            diagPart <- split(rep(diag(y@matLeft), each=ncol(y@matRight)), 
                              rep(1:ncol(y@matRight), each=ncol(y@matLeft)))
            bdMatrix(lapply(1:length(x), function(i) crossprod(x[[i]], 
                                                               (Diagonal(diagPart[[i]]) + 
                                                                  y@matRight))))
            
          } )


