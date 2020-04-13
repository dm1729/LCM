Shuffle <- function(x,i){ #input vector of length 2^d, shuffles every 2^i things for i<d
  #e.g. (1,2,3,4,5,6,7,8) can do (2,1,4,3,6,5,8,7) or (3,4,1,2,7,8,5,6) or (5,6,7,8,1,2,3,4) (4,2 or 1 swap)
  d <- round(log2(length(x))) # can add error check that i<d
  y <- x
  for (k in c(1:2^(d-i))){
    y[ c( 1:(2^(i-1)) ) + (2^i)*(k-1) ] <- x[ c( 1:(2^(i-1)) ) + (2^i)*(k-1) + 2^(i-1) ]# y[1]<-x[2] if i=1, k=1
    y[ c( 1:(2^(i-1)) ) + (2^i)*(k-1) + 2^(i-1) ] <- x[ c( 1:(2^(i-1)) ) + (2^i)*(k-1) ]# y[2]<-x[1]
  }
  return(y)
}

isSignedMTP2 <- function(p){ #
  library(MLLPs)
  library(rje)
  library(IsingIPS)
d <- round(log2(length(p)))
# NEED TO FIT 2^(n-1) different blobs (because flipping everything in an MTP2 dist keeps is MTP2)
C <- IsingDesign(d)[1:2^(d-1),] #Only want half of them (never flipping)
D <- c(1:(d-1))
status=FALSE
for (i in 2^{d-1}){ #Choose sign flip perm
  if {status=FALSE}{
  q <- p #reset to original input before flipping
    for (j in D[C[,i]==-1]){ #Which things are we flipping
      #FLIP
      q <- Shuffle(q,j)
      status=all(MTP2ineqs(d)%*%subsetMatrix(d)%*%log(q)>=0) # check if flipped version is MTP2
    }
  }
}
return(status)    
}
