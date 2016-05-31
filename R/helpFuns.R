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



isBlockDiag <- function(mat)
{
  
  stopifnot(mat[1,1] != 0)
  
  end <- 0
  nr <- nrow(mat)
  
  for(i in 1:ncol(mat)){
    
    sta <- end + 1 
    ei <- which(mat[sta:nr,i]==0)
    ei <- ifelse(length(ei)==0, length(sta:nr), max(min(ei) - 1, 1))
    end <- (sta:nr)[ei]
    # check column i  
    if(sum(mat[-1*(sta:end),i])!=0) return(FALSE)
    # check corresponding rows
    if(sum(as.vector(mat[sta:end, -i]))!=0) return(FALSE)
    
  }
  
  return(TRUE)
  
}