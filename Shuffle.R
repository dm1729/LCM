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