library(MASS)
library(Matrix)
library(pracma)

WS <- function(x,lb,ub){
  if(is.null(x)){
    RR <- median(x,lb,ub)
  } else {
    RR <- apply(cbind(x,lb,ub),1,median)
  }
  return(RR)
}

expit <- function(x){
  exp(x)/(1+exp(x))
}

MM <- function(X.Input){
  if( is.null(dim(X.Input)) ){
    X.Output <- matrix(X.Input,length(X.Input),1)
  } else {
    X.Output <- matrix(0,dim(X.Input)[1],dim(X.Input)[2])
    for(tt in 1:dim(X.Input)[2]){
      X.Output[,tt] <- as.numeric( X.Input[,tt] )
    }
  }
  return(X.Output)
}

FT_RBF <- function(X,X.new=NULL,bw.median){
  
  X <- MM(X)
  
  if(is.null(X.new)){
    X.new <- X
  } else {
    X.new <- MM(X.new)
  }
  DM <- matrix(0,dim(X.new)[1],dim(X)[1])
  for(ii in 1:dim(X.new)[1]){
    for(jj in 1:dim(X)[1]){
      DM[ii,jj] <- exp(-sum( (X.new[ii,]-X[jj,])^2/bw.median ))
    }
  }
  return(DM)
}

FT_DISC <- function(X,X.new=NULL,bw){
  
  X <- MM(X)
  if(is.null(X.new)){
    X.new <- X
  } else {
    X.new <- MM(X.new)
  }
  DM <- matrix(0,dim(X.new)[1],dim(X)[1])
  for(ii in 1:dim(X.new)[1]){
    for(jj in 1:dim(X)[1]){
      INDICATOR <- as.numeric( sum((X.new[ii,]-X[jj,])^2)==0 )
      DM[ii,jj] <- INDICATOR + bw*(1-INDICATOR)
    }
  }
  return(DM)
}

FT_BWHeuristic <- function(X,p=0.5){
  X <- MM(X)
  c( quantile( c(dist(X))^2, p ) )
}

FT_PMMR <- function(Y,
                    Perturb,
                    Target,
                    Diagonal=NULL,
                    Perturb.bw,
                    Target.bw,
                    lambda,
                    NV=FALSE,
                    intercept=TRUE){
  
  n <- dim(Perturb)[1]
  
  if(NV==TRUE & n>200){
    
    n <- dim(Perturb)[1]
    if(is.null(Diagonal)){
      Diagonal <- rep(1,n)
    } 
    D <- diag(Diagonal)
    
    Intercept <- mean(Y)/mean(Diagonal)*as.numeric(intercept)
    
    K.PP <- FT_RBF(Perturb, Perturb, Perturb.bw)
    
    nv <- 200
    
    XX <- matrix(0,n,n)
    Total.nv <- max(round(n/200),1)
    
    for(nv.iter in 1:Total.nv){
      Index <- sample(n,nv,replace = F)
      V <- ginv(FT_RBF(Target[Index,],  Target[Index,],  Target.bw ))
      U <- FT_RBF(Target[Index,],  Target,  Target.bw )
      # print(dim(XX))
      XX <- XX + (1/lambda)*( (diag(rep(1,n)) - U%*% ginv( (t(U)%*%K.PP%*%U)/lambda + ginv(V) )%*%t(U)%*%K.PP/lambda ))%*%(U%*%V%*%t(U)) 
    }
    
    XX.mean <- XX/Total.nv
    alpha <- XX.mean%*%(Y-Diagonal*Intercept)
    
  } else {
    
    n <- dim(Perturb)[1]
    if(is.null(Diagonal)){
      Diagonal <- rep(1,n)
    } 
    D <- diag(Diagonal)
    
    Intercept <- mean(Y)/mean(Diagonal)*as.numeric(intercept)
    
    K.PP <- FT_RBF(Perturb, Perturb, Perturb.bw)
    K.TT <- FT_RBF(Target,  Target,  Target.bw )
    
    
    alpha <- ginv(K.TT%*%D%*%K.PP%*%D%*%K.TT + lambda*K.TT) %*% (K.TT%*%D%*%K.PP%*%(Y-Diagonal*Intercept))
    
  }
  
  Result <- list()
  Result$alpha     <- alpha
  Result$intercept <- Intercept
  
  
  return( Result )
  
}

FT_PMMR_CV <-function(Y.Train,
                      Perturb.Train,
                      Target.Train,
                      Diagonal.Train=NULL,
                      Y.Valid,
                      Perturb.Valid,
                      Target.Valid,
                      Diagonal.Valid=NULL,
                      Perturb.bw,
                      Target.bw,
                      lambda,
                      NV=FALSE,
                      intercept=TRUE){
  
  n.Train <- dim(Perturb.Train)[1]
  n.Valid <- dim(Perturb.Valid)[1]
  
  if(is.null(Diagonal.Train)){
    Diagonal.Train <- rep(1,n.Train)
  } 
  if(is.null(Diagonal.Valid)){
    Diagonal.Valid <- rep(1,n.Valid)
  } 
  D.Train <- diag(Diagonal.Train)
  D.Valid <- diag(Diagonal.Valid)
  
  Result <- FT_PMMR(Y               = Y.Train,
                    Perturb         = Perturb.Train,
                    Target          = Target.Train,
                    Diagonal        = Diagonal.Train,
                    Perturb.bw      = Perturb.bw,
                    Target.bw       = Target.bw,
                    lambda          = lambda,
                    NV              = NV,
                    intercept       = intercept ) 
  
  # K.PP.Train <- FT_RBF(Perturb.Train, Perturb.Train, Perturb.bw)
  K.PP.Valid <- FT_RBF(Perturb.Valid, Perturb.Valid, Perturb.bw)
  
  # K.TT.Train <- FT_RBF(Target.Train,  Target.Train,  Target.bw )
  K.TT.Cross <- FT_RBF(Target.Valid,  Target.Train,  Target.bw )
  
  h.Valid.hat <- t(K.TT.Cross)%*%Result$alpha + Result$intercept
  residual.Valid <- Y.Valid - D.Valid%*%c(h.Valid.hat)
  V.statistic <- t(residual.Valid)%*%K.PP.Valid%*%(residual.Valid)/n.Valid^2 
  U.statistic <- (V.statistic - sum((residual.Valid)^2*diag(K.PP.Valid))/n.Valid^2)*(n.Valid^2)/(n.Valid)/(n.Valid-1)
  
  RESULT <- list()
  RESULT$lambda <- lambda
  RESULT$Vstat  <- V.statistic
  RESULT$Ustat  <- U.statistic
  return(RESULT)
  
} 


################################################################################
# Auxiliary Functions
################################################################################

CrossVar <- function(PARAMETER,
                     ss,
                     subset=posA1D0,
                     response=Y*(1-D)*A,
                     target=WX,
                     perturb=ZX,
                     diagonal=(1-D)*A ){
  
  P.P <- rep(as.numeric(PARAMETER[1]),ncol(perturb))
  P.T <- rep(as.numeric(PARAMETER[2]),ncol(target))
  Lambda  <- as.numeric(PARAMETER[3])
  
  RISK <- rep(0,NumCVRep*NumCV)
  
  for(cv in 1:(NumCVRep*NumCV)){
    
    CV.Split.Index <- list()
    CV.Split.Index[[1]] <- intersect( SS.Index[[ss]][ -CV.Index[[ss]][[cv]] ],subset)
    CV.Split.Index[[2]] <- intersect( SS.Index[[ss]][ CV.Index[[ss]][[cv]] ] ,subset)
    
    response.CV  <- list()
    target.CV    <- list()
    perturb.CV   <- list()
    diagonal.CV   <- list()
    
    for(cvest in 1:2){
      response.CV[[cvest]] <- response[CV.Split.Index[[cvest]]]
      target.CV[[cvest]]   <- target[CV.Split.Index[[cvest]],]
      perturb.CV[[cvest]]  <- perturb[CV.Split.Index[[cvest]],]
      diagonal.CV[[cvest]] <- diagonal[CV.Split.Index[[cvest]]]
    }
    
    CV.result <- 
      FT_PMMR_CV(Y.Train       =response.CV[[1]],
                 Perturb.Train =perturb.CV[[1]],
                 Target.Train  =target.CV[[1]],
                 Diagonal.Train=diagonal.CV[[1]],
                 Y.Valid       =response.CV[[2]],
                 Perturb.Valid =perturb.CV[[2]],
                 Target.Valid  =target.CV[[2]],
                 Diagonal.Valid=diagonal.CV[[2]],
                 Perturb.bw    =exp(P.P),
                 Target.bw     =exp(P.T),
                 lambda        =exp(Lambda),
                 NV            =FALSE)
    
    RISK[c(cv)] <- c(CV.result$Vstat)
  }
  
  return( mean(RISK) )
}
