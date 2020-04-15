LCMineqs <- function (n){ #based on https://arxiv.org/pdf/1710.01696.pdf and https://projecteuclid.org/download/pdf_1/euclid.aos/1015956713
  #Gives list of permissible simulatenous zeros of MTP2 constraints.
  #WIll return allowed row numbers for rho matrix
  U <- c(1:choose(n,2)) # To choose argument of rho (subset of size 2, ordering from paper)
  V <- c(1:2^{n-2}) # To choose sign pattern fo rho (binary order)
  #N Total number of available strata. Form a matrix with N rows which are indicator vectors of length U*V "1 = exact"kl
  
  #CALCULATE N
  L <- c(3:(n-1)) #l values
  K<-L #preallocate
  for (i in L){
    K[i-2]=choose(n,i)*(2^i)
  }
  N <-  sum(K)+2*n + 4*choose(n,2) + n + 2^n + choose(n,2) + 1 + 1 #add l=1, l=2, l=n, l=n+1, add also one for interior
  #N<-1000 # Extra allocation incase counter is wrong
  M <- matrix(0,nrow=N,ncol=choose(n,2)*2^(n-2)) # N rows for each strata, length(U)*length(V) many contraints
  j<-1 # Counter to populate rows of M
  
  #TYPE 1 STRATA
  #NO CHANGE TO FIRST ROW of M AS NO EXACT EQUALITIES
  j<-j+1
  
  #TYPE 2 STRATA of dimension 2n
  #|I|=1  ROWS 2 to 9
  A<-combn(n,2) #TO ALLOW US TO LOOK UP U COORDINATE WITH PAIR
  A[1,]<-rev(A[1,])
  A[2,]<-rev(A[2,]) #TO GIVE ORDERING OF PAPER
  B <- IsingIPS::IsingDesign(n-2) # to reference orderings for excluded nodes (Need IsingIPS)
  B <- (B+3)/2 # To replace -1 with 1 and with 2 (LOW/HIGH ordering)
  
  
  for ( i in c(1:n) ){ #APPEARS TO BE WORKING
    #CHOOSE ONE TO INCLUDE in I
    #GIVE U,V COORDINATES
    for (k in c(1,2)){
      for (u in U){
        if ( all(A[,u]!=i) ){ # If the choice 'u' of A column doesn't include the chosen 'conditioning node'
          #WHICH ONES TO INCLUDE. Update row correponding to (u,v) entries
          C <- which(c(1:n)[-A[,u]]==i) #Position if i in remaining n-2 points
          for (v in V){
            if (B[v,C]==k) { #Checks whether the conditioned n-2 tuple conditions on X_i=k
              M[j,2^(n-2)*(u-1)+v] <- 1 # As in this case we have independence
              #The entry of M maps the (u,v) pair to the row of the big constraint matrix
              #(Which corresponds to the rows in MTP2inequaliites, or the rho in B+F)
            } 
            
          }#CLOSE v in V
          
        }
      }# CLOSE u in U
      j <- j+1 #CHANGE of k means new strata
    }
  }
  
  #|I|=2 (different one) ROWS 10-33 inclusive (dimension 2n-1)
  for (i in c(1:choose(n,2))){
    a <- A[1,i]
    b <- A[2,i]
    for (k in c(1:4)){
      s <- combinat::hcube(rep(2,2))[k,] #Gives pattern for vals of entries of I (in label order)
      for (u in U){
        if (u!= i){ #Only pair to exclude is one with both (Had to pick a and b in reverse order to make work)
          #Do not update u=i as this corresponds to interaction between a,b which can be nonzero for theta_{ab} parameter
          if ( any(A[,u]==a) ){ #Does chosen pair (number u) contain a?
            C <- which(c(1:n)[-A[,u]]==b) # a excluded, gives position of b in remainder
            for (v in V){
              if (B[v,C]==s[1]){
                M[j,2^(n-2)*(u-1)+v] <- 1
              }
            } #END LOOP OVER V
          }
          else if ( any(A[,u]==b) ){ #Does chosen pair (number u) contain b?
            C <- which(c(1:n)[-A[,u]]==a) # b excluded (and other of log odds pair), gives position of a in remainder
            for (v in V){
              if (B[v,C]==s[2]){
                M[j,2^(n-2)*(u-1)+v] <- 1
              }
            }
          }
        }
      } #END LOOP OVER U
      j<- j+1 #MOVE TO NEXT STRATA BY ALTERING i or k
    }
  }
  
  
  #TYPE 3 STRATA #CONDITIONAL INDEPENDECE GIVEN THIS VARIABLE (dimension 2n-1)
  #ROWS 82-85
  
  for (i in c(1:n))
  {
    for (u in U){ #(Might be a cleaner way to do this, just want to set all 'columns' to 1 for all rows not including i)
      if ( all(A[,u]!=i) ){
        M[j,(2^(n-2)*(u-1)+1):(2^(n-2)*u)] <- 1
      }
      
    }
    j<-j+1
  }
  
  #TYPE 2 STRATA |I|>2 ROWS 34-81 inclusive - nothing populating! Now too much populating! (dimensions 2n-2 to n+1, descending)
  for (I in c(3:n)){ # Set size of |I|
    for (i in c(1:choose(n,I)) ){# Subsets of size I (=|I|)
      if (I==n){ #seemed to be getting issues because combn(n,n) is a vector not a matrix, so added this loop
        S<-c(1:n) #Gives only subset of size n
      }
      else{
      S <-combn(n,I)[,i] #Gives the subset
      }
      for (k in c(1:2^I)){
        s <- combinat::hcube(rep(2,I))[k,] #Gives pattern for vals of entries of I (in label order)
        for (u in U){
          C <- c(1:n)[-A[,u]]
          for (jj in (S)){
            if (any(C==jj)){ #If the index jj in I (S) is not in the chosen pair (so can be conditioned upon)
              for (v in V){
                if (B[v,which(C==jj)]==s[which(S==jj)]){ #Add in the one where the conditoning (n-2) tuple has X_j=s_j
                  M[j,(2^(n-2))*(u-1)+v] <- 1
                }
              }  #END LOOP OVER v
            }
          } #END LOOP OVER jj in S
          
        } #END LOOP OVER U
        j <- j+1  
      } #END LOOP OVER CHANGING CONDITIONING VALS
    }
  }

  
  #TYPE 4 STRATA #IMAGES OF THETA_ij (THETA_ab) (dimension n+1)
  #ROWS 86-91
  
  for (i in c(1:choose(n,2))){
    a <- A[1,i]
    b <- A[2,i]
    for (u in U){
      if (u!= i){ #Only pair to exclude is one with both UPDATED TO BE CORRECT
          M[j,(2^(n-2)*(u-1)+1):(2^(n-2)*u)] <- 1
        }
      
    } #END LOOP OVER U
    j<- j+1 #MOVE TO NEXT STRATA BY ALTERING i or k
  }
  #TYPE 5 STRATA ROW 92 (dimension n)
  M[j,] <- rep(1,length(U)*length(V)) # FULL INDEPENDENCE
  return(M)
}