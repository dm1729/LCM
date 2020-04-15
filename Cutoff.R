Cutoff <- function(n){ #Gives cutoff points for different dimension strata in M_n
  #n to 2n+1 inclusive = n+2 cutoffs. Vector has last strata number of that dimension
C <- rep(0,n+2)
C[1]=1 # Last dimension 2n+1 strata is interior
for (l in c(1:(n+1))){ #Using Corollary 10, dimension k = 2n + 1 - l (so C[1] for k=0, C[k+1] for k=k)
  if (l==2){
    C[l+1] <- C[l] + choose(n,l)*2^l + n
  }
  else if (l==n){
    C[l+1] <- C[l] + choose(n,l)*2^l+ choose(n,2)
  }
  else if (l==n+1){
  C[l+1] <- C[l] + 1
  }
  else{
  C[l+1] <- C[l] + choose(n,l)*2^l #How many along from previous strata? Have extra if loop for adding in exceptional strata
  }
}
return(C) #e.g. n=4 says strata 1 is last of 9, strata 9 for 8, strata 37 for 7, strata 69 for 6, strata 91 for 5, strata 92 for 4
}