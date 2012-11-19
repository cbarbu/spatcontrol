".spam.elemul" <- function(e1,e2)
{
  if(is.vector(e1)) {
    if(length(e1) == 1){
      if(e1==0) return( spam(0,nrow(e2),ncol(e2)))
      else{  # just a scalar
        e2@entries <- e1*e2@entries
        return(e2)
      }
    }  else if(length(e1) == nrow(e2))
      return(diag.spam(e1) %*% e2)
    else # length(e1) == ncol(e2) is not required
      stop("e1 and e2 not conformable for efficient element-by-element multiplication")
  }
  else if(is.vector(e2)) {
    if(length(e2) == 1){
      if(e2==0)   return( spam(0,nrow(e1),ncol(e1)))
      else {
        e1@entries <- e2*e1@entries
        return(e1)
      }
    }
    else if(length(e2) == nrow(e1))
      return(diag.spam(e2) %*% e1)
    else
      stop("e1 and e2 not conformable for efficient element-by-element multiplication")
  }
  if(is.matrix(e1))
    e1 <- as.spam(e1)
  else if(is.matrix(e2))
    e2 <- as.spam(e2)
  if(!(is.spam(e1) && is.spam(e2)))
    stop("Arguments must be of class:  vector, matrix or spam")
  
  e1row <- e1@dimension[1]
  e1col <- e1@dimension[2]
  if(e1col != e2@dimension[2] | e1row != e2@dimension[1])
    stop("non-conformable matrices")
  nnzmax <- length(intersect(e1@colindices+e1col*(rep(1:e1row,diff(e1@rowpointers))-1),
                             e2@colindices+e2@dimension[2]*(rep(1:e2@dimension[1],diff(e2@rowpointers))-1)))+1
  z <- .Fortran("aemub",
                e1row,
                e1col,
                dcheck(e1@entries),
                e1@colindices,
                e1@rowpointers,
                dcheck(e2@entries),
                e2@colindices,
                e2@rowpointers,
                entries     = vector("double",nnzmax),
                colindices  = vector("integer",nnzmax),
                rowpointers = vector("integer",e1row+1),
                integer(e1col),
                double(e1col),
                as.integer(nnzmax),
                ierr = vector("integer",1),
                DUP = FALSE,
                PACKAGE = "spam"
                )
  if(z$ierr != 0)      stop("insufficient space for element-wise sparse matrix multiplication")
  nnz <- z$rowpointers[e1row+1]-1
  if(identical(z$entries,0)){#trap zero matrix
    z$colindices <- 1L
    z$rowpointers <- c(1L,rep(2L,e1row))
  }
  return(new("spam",entries=z$entries[1:nnz],colindices=z$colindices[1:nnz],rowpointers=z$rowpointers,dimension=c(e1row,e1col)))
}


setMethod("*",signature(e1="spam",e2="spam"), .spam.elemul)
setMethod("*",signature(e1="spam", e2="ANY"), .spam.elemul)
setMethod("*",signature(e1="ANY", e2="spam"), .spam.elemul)
