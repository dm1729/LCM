##' Get MTP2 inequalities for binary distribution FROM MLLPs PACKAGE
##' 
##' @param d number of variables
##' 
##' @details Uses method in Bartolucci and Forcina (2000).
##' 
##' @import Matrix
MTP2ineqs <- function (d) {
  B0 <- combn(d, 2)
  B <- matrix(0, ncol(B0), d)
  for (i in seq_len(ncol(B0))) B[i,B0[,i]] <- 1
  c <- nrow(B)
  
  out <- Matrix(0, nrow=c*2^(d-2), ncol=2^d)
  
  ## basic operators
  e <- Matrix(c(0, 1), nrow=1, sparse=TRUE)
  E <- Matrix(c(1, 0, 1, 1), 2, 2, byrow=TRUE, sparse=TRUE)
  
  for (i in seq_len(c)) {
    A <- if (B[i,1] == 0) E
    else e
    
    for (j in seq_len(d)[-1]) {
      A <- if (B[i,j] == 0) kronecker(A, E)
      else kronecker(A, e)
    }
    
    out[(i-1)*2^(d-2) + seq_len(2^(d-2)), ] <- A
  }

  out  
}


##' Fit binary distribution to MTP2 constraints
##' 
##' Finds MLE for a vector of (presumed) counts.
##' 
##' @param y vector of counts to find MLE of
##' @param p starting point for \code{p}
##' @param control list of control parameters (see Details)
##' 
##' @details Finds inequality constraints using \code{MTP2ineqs}.
##' Control parameters include:
##' \itemize{
##' \item \code{max.it} maximum number of iterations;
##' \item \code{eps} sqrt machine precision (no effect currently);
##' \item \code{tol} convergence tolerance.
##' }
##' 
##' @importFrom MASS ginv
##' @importFrom osqp osqp
##' 
##' @export MTP2solve
FaceSolve <- function(y, z, p, control=list()) { #z==0 not allowed, for interior just use MTP2solve
  requireNamespace("CVXR", quietly = TRUE)
  
  n <- sum(y)
  if (missing(p)) p <- y/n #p <- (y+1)/sum(y+1)
  p0 <- y/n

  lo <- sum((y*log(p0))[y > 0])
  t <- length(y)
  d <- log2(t)
  
  # H <- MASS::ginv(diag(t)[,-1])
  
  D <- MTP2ineqs(d)[,-1]
  C <- subsetMatrix(d)[-1,]
  
  CM <- list(C=C, M=diag(t))
  
  # Get control parameters or use defaults
  con = list(max.it = 250, eps = .Machine$double.eps^0.5, tol = 1e-10)
  matches = match(names(control), names(con))
  con[matches] = control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ", 
                                   paste(names(control[is.na(matches)]), sep = ", "))
  
  theta <- log(p[-1])
  eta <- C %*% log(p)

  y1 <- y[-1]
  R <- C[,-1]
  Ri <- solve.default(R)
  Rit <- t(Ri)
  Rity <- Rit %*% y1
  
  ## for CVXR
  x <- Variable(length(theta))
  
  cont <- TRUE
  its <- 0
  
  while (cont) {
    p1 <- p[-1]
    Rp <- Rit %*% p1
    sc <- Rity - n*Rp

    DRp <- Rit %*% (c(p1)*Ri)
    FI <- n*(DRp - Rp %*% t(Rp))
    C <- chol(FI)
    
    theta0 <- eta + solve.default(FI, sc)
    # theta0 <- eta + backsolve(C, forwardsolve(C, sc, upper.tri = TRUE, transpose = TRUE))
        
    ## solve QP - try to speed this up!?
    obj <- Minimize(sum((C %*% (x-theta0))^2))
    prob <- Problem(obj, constraints = list( (D %*% x) [which(z==1),1] == 0)) #z 'face vector' from LCMineqs
    result <- psolve(prob, abstol = con$tol)
    tmp <- c(result$getValue(x))
    
    ## move
    dis <- sum(abs(tmp - eta))
    eta <- tmp
    
    its <- its + 1
    cont <- (dis > con$tol) && (its < con$max.it)
  }
  
  if (dis > con$tol) error = 1
  else error = 0
  
  ## compute log-likelihood
  p <- invMLLP(eta, CM, p)$p
  ll <- sum((y*log(p))[y > 0])
  dev <- 2*(lo - ll)
  rho <- MTP2ineqs(d)%*%subsetMatrix(d)%*%log(p)
  rho[abs(rho)<1e-8]<-0 #Gets rid of spurious nonzero
  boundary <- (rho==0)*1
  
  return(list(p=p, eta=eta, error=error, its=its, ll=ll, dev=dev,rho=rho,boundary=boundary))
}


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


LCMsolve <- function(y){
  library(CVXR)
  library(rje)
  library(MLLPs)
  library(IsingIPS)
  library(Matrix)
  d <- round(log2(length(y)))
  M <- LCMineqs(d)
  L <- length(M[,1])
  
  P <- matrix(0,nrow=(L),ncol=(2^d+2))
  
  l <- 1
  x <- isSignedMTP2(y/sum(y)) #Running facesolve unconstrained returns emperical dist.
  if (x==TRUE){
    P[l,1:2^d] <- y/sum(y) #Stores signed MTP2 dist
    P[l,(2^d+1)] <- FlatteningRank(P[l,1:2^d]) #Stores associated flatteningrank Technically can exit here if LCM.
    v <- ( abs( ( MTP2ineqs(d)%*%subsetMatrix(d)%*%log(P[l,1:2^d]) )[,1] ) <1e-8 ) #Gives sign pattern of p in form of M matrix
    P[l,(2^d+2)] <- any ( apply( t(M) == v , 2, all ) ) #Good pattern?
    }
 
  l <- 2 #Had some issues starting in interior and I think there's essentially no point.
 
  while (l <= L){
    z <- M[l,]
    p <- FaceSolve(y,z)$p
    x <- isSignedMTP2(p) #true/false
      if (x==TRUE){
        P[l,1:2^d] <- p #Stores signed MTP2 dist
        P[l,(2^d+1)] <- FlatteningRank(P[l,1:2^d]) #Stores associated flatteningrank
        v <- ( abs( ( MTP2ineqs(d)%*%subsetMatrix(d)%*%log(P[l,1:2^d]) )[,1] ) <1e-8 ) #Gives sign pattern of p in form of M matrix
        P[l,(2^d+2)] <- any ( apply( t(M) == v , 2, all ) ) #Good pattern?
        if (P[l,(2^d+1)] <= 2){#L to be the upper bound of this dimension if a FR 2 estimate has been found.
        L <- Cutoff(d)[which(Cutoff(d)>=l)[1]] # use >= because if we are already on last strata of dimension we stop.
        }
      }
    l <- l+1
  } #end while
  return(list(PMFs=P[1:L,],Strata=(l-1))) #Strata l-1 because the final one won't be looped over.
}
