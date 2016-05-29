getBlockIndices <- function(x)
{
  
  stopifnot(class(x)=="bdMatrix")
  
  indRend <- cumsum(sapply(x@listOfBlocks, NROW))
  indCend <- cumsum(sapply(x@listOfBlocks, NCOL))
  indRsta <- c(1, indRend[-length(indRend)]+1)
  indCsta <- c(1, indCend[-length(indCend)]+1)
  
  indRsplit <- lapply(1:length(indRend), function(k) indRsta[k]:indRend[k])
  indCsplit <- lapply(1:length(indCend), function(k) indCsta[k]:indCend[k])
  
  return(list(rows=indRsplit, cols=indCsplit))
  
}


hasEqualBlocksizes <- function(x) if(class(x)=="bdMatrix") 
  (length(unique(sapply(x@listOfBlocks,NCOL)))==1 & length(unique(sapply(x@listOfBlocks,NROW)))==1) else
    stop("Method only implemented for class bdMatrix")


