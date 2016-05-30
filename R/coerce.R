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

as.vector.bdMatrix <- function(x, mode) as(x, "vector")

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