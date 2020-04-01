LCMineqs <- function (n) 
#Gives list of permissible simulatenous zeros of MTP2 constraints.
#WIll return allowed row numbers for rho matrix
U <- c(1,:choose(n,2)) # To choose argument of rho (subset of size 2, ordering from paper)
V <- c(1:2^{n-2}) # To choose sign pattern fo rho (binary order)
#N Total number of available strata. Form a matrix with N rows which are indicator vectors of length U*V "1 = exact"kl

#CALCULATE N
L <- c(2,n-1) #l values
K<-L #preallocate
for (i in L){
  K[i]=choose(n,i)*2*i
}
N <-  sum(K)+4*choose(n,2) +n + 2^n +choose(n,2) + 1 #add l=2, l=n, l=n+1
M <- matrix(0,nrow=N,ncol=choose(n,2)*2&{n-2}) # N rows for each strata, length(U)*length(V) many contraints
j<-1 # Counter to populate rows of M

#TYPE 1 STRATA
#NO CHANGE TO FIRST ROW of M AS NO EXACT EQUALITIES
j<-j+1

#TYPE 2 STRATA
#|I|=1 
A<-combn(n,2) #TO ALLOW US TO LOOK UP U COORDINATE WITH PAIR
A[1,]<-rev(A[1,])
A[2,]<-rev(A[2,]) #TO GIVE ORDERING OF PAPER

for ( i in c(1:n)){
  #CHOOSE ONE TO INCLUDE in I
  #GIVE U,V COORDINATES
  for (k in c(1,2)){
    
  }
    #s=1 case, conditional independence given X_i=1
    
}

#|I|=2 (different one)
for (i in c(1:combn(n,2))){
  a <- A[1,i]
  b <- A[2,i]
}

#|I|>2
for (I in c(3:n)){ # Set size of |I|
  for (i in combn(n,I)){# Subsets of size I (=|I|)
    for (k in c(1,2^i)){
      s <- combinat::hcube(rep(2,n))[k,] #Gives pattern for vals of entries of I (in label order)
      
    }
  }
}

#TYPE 3 STRATA

for (i in c(1:n))
{
  
}
  #CONDITIONAL INDEPENDECE GIVEN THIS VARIABLE

#TYPE 4 STRATA

for (i in c(1:combn(n,2))){
  a <- A[1,i]
  b <- A[2,i]
}
#TYPE 5 STRATA
M[j,] <- rep(1,length(U)*length(V)) # FULL INDEPENDENCE