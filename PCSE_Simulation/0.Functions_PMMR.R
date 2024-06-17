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
  
  if(NV==TRUE){
    
    n <- dim(Perturb)[1]
    if(is.null(Diagonal)){
      Diagonal <- rep(1,n)
    } 
    D <- diag(Diagonal)
    
    Intercept <- mean(Y)/mean(Diagonal)*as.numeric(intercept)
    
    K.PP <- FT_RBF(Perturb, Perturb, Perturb.bw)
    nv <- ifelse(n>100,100,n)
    
    XX <- matrix(0,n^2,2)
    for(nv.iter in 1:dim(XX)[2]){
      Index <- sample(n,nv,replace = T)
      V <- ginv(FT_RBF(Target[Index,],  Target[Index,],  Target.bw ))
      U <- FT_RBF(Target[Index,],  Target,  Target.bw )
      XX[,nv.iter] <-
        as.vector( (1/lambda)*( (diag(rep(1,n)) - U%*% ginv( (t(U)%*%K.PP%*%U)/lambda + ginv(V) )%*%t(U)%*%K.PP/lambda ))%*%(U%*%V%*%t(U)) )
    }
    XX.mean <- matrix(apply(XX,2,mean),n,n)
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

CrossVal <- function(PARAMETER,
                     ss,
                     P.BW,
                     subgroup,
                     response,
                     diagonal,
                     target="WX"){
  
  Scale.P <- PARAMETER[1]
  Scale.T <- PARAMETER[2]
  Lambda  <- PARAMETER[3]
  
  RISK <- rep(0,NumCVRep*NumCV)
  
  for(cv in 1:(NumCVRep*NumCV)){
    
    CV.Split.Index <- list()
    CV.Split.Index[[1]] <- intersect( SS.Index[[ss]][ -CV.Index[[cv]] ] , subgroup)
    CV.Split.Index[[2]] <- intersect( SS.Index[[ss]][ CV.Index[[cv]] ]  , subgroup)
    
    # CV.Split.Index <- list()
    # CV.Split.Index[[1]] <- intersect( SS.Index[[ss]][ -CV.Index[[cv]] ] , SS.Index[[ss]])
    # CV.Split.Index[[2]] <- intersect( SS.Index[[ss]][ CV.Index[[cv]] ]  , SS.Index[[ss]])
    
    Y.CV  <- list()
    A.CV  <- list()
    D.CV  <- list()
    X.CV  <- list()
    W.CV  <- list()
    Z.CV  <- list()
    WX.CV <- list()
    ZX.CV <- list()
    
    Response.CV <- list()
    Diagonal.CV <- list()
    
    for(cvest in 1:2){
      Y.CV [[cvest]] <- Y[CV.Split.Index[[cvest]] ]
      A.CV [[cvest]] <- A[CV.Split.Index[[cvest]] ]
      D.CV [[cvest]] <- D[CV.Split.Index[[cvest]] ]
      X.CV [[cvest]] <- X[CV.Split.Index[[cvest]],]
      W.CV [[cvest]] <- W[CV.Split.Index[[cvest]] ]
      Z.CV [[cvest]] <- Z[CV.Split.Index[[cvest]] ]
      WX.CV[[cvest]] <- cbind(W.CV[[cvest]],X.CV[[cvest]])
      ZX.CV[[cvest]] <- cbind(Z.CV[[cvest]],X.CV[[cvest]])
      
      Response.CV[[cvest]] <- response[CV.Split.Index[[cvest]] ]
      Diagonal.CV[[cvest]] <- diagonal[CV.Split.Index[[cvest]] ]
    }
    
    if(target=="WX"){
      P.Mat <- ZX.CV
      T.Mat <- WX.CV
      P.BW.Fit <- as.numeric(ZX.bw*(P.BW)*10^Scale.P)
      T.BW.Fit <- as.numeric(WX.bw*10^Scale.T)
    } else {
      P.Mat <- WX.CV
      T.Mat <- ZX.CV
      P.BW.Fit <- as.numeric(WX.bw*(P.BW)*10^Scale.P)
      T.BW.Fit <- as.numeric(ZX.bw*10^Scale.T)
    }
    
    ## h11
    
    CV.Temp <- FT_PMMR_CV( Y.Train        = Response.CV[[1]],
                           Perturb.Train  = P.Mat[[1]],
                           Target.Train   = T.Mat[[1]],
                           Diagonal.Train = Diagonal.CV[[1]],
                           Y.Valid        = Response.CV[[2]],
                           Perturb.Valid  = P.Mat[[2]],
                           Target.Valid   = T.Mat[[2]],
                           Diagonal.Valid = Diagonal.CV[[2]],
                           Perturb.bw     = P.BW.Fit,
                           Target.bw      = T.BW.Fit,
                           lambda         = 10^Lambda,
                           NV             = NV.Type,
                           intercept      = Intercept.Type)
    
    RISK[c(cv)] <- c(CV.Temp$Ustat) + penalty(Scale.P , Scale.T , Lambda )
  }
  
  return( mean(RISK) )
}



GMM.Moment <- function(Y,
                       hat.matrix,
                       ss,
                       subgroup,
                       target="WX"){
  
  sub.index <- intersect( SS.Index[[ss]] , subgroup)
  
  Y.Temp  <- Y[sub.index]
  Z.Temp  <- Z[sub.index]
  W.Temp  <- W[sub.index]
  X1.Temp <- X1[sub.index]
  X2.Temp <- X2[sub.index]
  
  if(target=="WX"){
    coef <- as.matrix( cbind(1,
                             Z.Temp,X1.Temp,X2.Temp,
                             Z.Temp^2,X1.Temp^2,X2.Temp^2,
                             Z.Temp*X1.Temp,Z.Temp*X2.Temp,X1.Temp*X2.Temp,
                             as.matrix(model.matrix(~cut( Z.Temp,breaks=quantile(  Z.Temp,seq(0,1,0.25)),include.lowest = T))[,-1]),
                             as.matrix(model.matrix(~cut(X1.Temp,breaks=quantile( X1.Temp,seq(0,1,0.25)),include.lowest = T))[,-1]),
                             as.matrix(model.matrix(~cut(X2.Temp,breaks=quantile( X2.Temp,seq(0,1,0.25)),include.lowest = T))[,-1])) )
    
  } else {
    coef <- as.matrix( cbind(1,
                             W.Temp,X1.Temp,X2.Temp,
                             W.Temp^2,X1.Temp^2,X2.Temp^2,
                             W.Temp*X1.Temp,W.Temp*X2.Temp,X1.Temp*X2.Temp,
                             as.matrix(model.matrix(~cut( W.Temp,breaks=quantile(  W.Temp,seq(0,1,0.25)),include.lowest = T))[,-1]),
                             as.matrix(model.matrix(~cut(X1.Temp,breaks=quantile( X1.Temp,seq(0,1,0.25)),include.lowest = T))[,-1]),
                             as.matrix(model.matrix(~cut(X2.Temp,breaks=quantile( X2.Temp,seq(0,1,0.25)),include.lowest = T))[,-1])) )
  }
  
  return( 
    sapply(1:dim(hat.matrix)[2],
           function(tt){
             res <- as.numeric(Y.Temp - hat.matrix[sub.index,tt])
             max(apply(coef*matrix(res,length(res),dim(coef)[2]),
                       2,
                       function(v){ t.test(v)$statistic^2 }))
           })
    )
    
}


FT_Adversarial <- function(Y,
                           Perturb,
                           Target,
                           Diagonal=NULL,
                           Perturb.bw,
                           Target.bw,
                           Target.lambda,
                           Perturb.lambda,
                           NV=FALSE,
                           intercept=TRUE){
  
  
  
  n <- dim(Perturb)[1]
  if(is.null(Diagonal)){
    Diagonal <- rep(1,n)
  } 
  D <- diag(Diagonal)
  
  Intercept <- mean(Y)/mean(Diagonal)*as.numeric(intercept)
  
  K.PP <- FT_RBF(Perturb, Perturb, Perturb.bw)
  K.TT <- FT_RBF(Target,  Target,  Target.bw )
  
  Gamma <- 0.25*K.PP%*%ginv(K.PP/n + diag(rep(Perturb.lambda,n)))
  gamma <- -ginv(K.TT%*%D%*%Gamma%*%D%*%K.TT+n^2*Target.lambda*K.TT)%*%(K.TT%*%D%*%Gamma%*%matrix((Y-Diagonal*Intercept)))
  
  Result <- list()
  Result$alpha     <- gamma
  Result$intercept <- Intercept
  
  return(Result)
  
}

FT_Adversarial_CV <- function(Y.Train,
                              Perturb.Train,
                              Target.Train,
                              Diagonal.Train=NULL,
                              Y.Valid,
                              Perturb.Valid,
                              Target.Valid,
                              Diagonal.Valid=NULL,
                              Perturb.bw,
                              Target.bw,
                              Perturb.lambda,
                              Target.lambda,
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
  
  Result <- FT_Adversarial(Y               = Y.Train,
                           Perturb         = Perturb.Train,
                           Target          = Target.Train,
                           Diagonal        = Diagonal.Train,
                           Perturb.bw      = Perturb.bw,
                           Target.bw       = Target.bw,
                           Perturb.lambda  = Perturb.lambda,
                           Target.lambda   = Target.lambda,
                           NV              = NV,
                           intercept       = intercept) 
  
  # K.PP.Train <- FT_RBF(Perturb.Train, Perturb.Train, Perturb.bw)
  K.PP.Valid <- FT_RBF(Perturb.Valid, Perturb.Valid, Perturb.bw)
  
  # K.TT.Train <- FT_RBF(Target.Train,  Target.Train,  Target.bw )
  K.TT.Cross <- FT_RBF(Target.Valid,  Target.Train,  Target.bw )
  
  h.Valid.hat <- t(K.TT.Cross)%*%Result$alpha + Result$intercept
  
  Error.Proj <- D.Valid%*%h.Valid.hat
  
  RESULT <- list()
  RESULT$PL <- Perturb.lambda
  RESULT$TL <- Target.lambda
  RESULT$risk <- Error.Proj
  return(RESULT)
} 
