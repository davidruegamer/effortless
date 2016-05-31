### bdMatrix

setAs(from = "list", 
      to = "bdMatrix", 
      def = function(from) {
        bdMatrix(from) 
      }
)

setAs(from = "bdMatrix", 
      to = "vector", 
      def = function(from) {
        unlist(lapply(from@listOfBlocks,as.vector)) 
      }
)

#'  Coerce bdMatrix to vector
#' 
#' @export
as.vector.bdMatrix <- function(x, mode) as(x, "vector")

setAs(from = "Matrix", 
      to = "bdMatrix", 
      def = function(from) {
        
        stopifnot(from[1,1] != 0)
        
        end <- 0
        nr <- nrow(from)
        resL <- vector("list", ncol(from))
        
        for(i in 1:ncol(from)){
          
          sta <- end + 1 
          ei <- which(from[sta:nr, i] == 0)
          ei <- ifelse(length(ei) == 0, length(sta:nr), max(min(ei) - 1, 1))
          end <- (sta:nr)[ei]
          # check column i  
          if(sum(from[-1*(sta:end),i]) != 0) stop("Matrix not block-diagonal.")
          # check corresponding rows
          if(sum(as.vector(from[sta:end, -i])) != 0) stop("Matrix not block-diagonal.")
          
          resL[[i]] <- from[sta:end, i]
          
        }
        
        return(bdMatrix(resL))
        
      }
)

### kroneckersumBlockMatrix

setAs(from = "list", 
      to = "kroneckersumBlockMatrix", 
      def = function(from) {
        kroneckersumBlockMatrix(from) 
      }
)

### rowtensorBlockMatrix


setAs(from = "list", 
      to = "rowtensorBlockMatrix", 
      def = function(from) {
        rowtensorBlockMatrix(from) 
      }
)