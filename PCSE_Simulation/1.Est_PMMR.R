################################################################################
# Cluster computing: BATCH = 1,...,16000
################################################################################

args <- (commandArgs(trailingOnly=TRUE))
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

################################################################################
# Basic Simulation Parameters
################################################################################

Sim.Setup.Obs <- (BATCH<=2000) 
# TRUE = Observational setting, FALSE = Experimental setting
seed.data <- (BATCH-1)%%2000+1
seed.CF   <- 1
N <- c(500,1000,1500,2000,
       500,1000,1500,2000)[(BATCH-1)%/%2000+1]
FOLDER <- "Result/"
set.seed( seed.data )

################################################################################
# Package and Source Files
################################################################################

source("0.Functions_PMMR.R")
source("0.MySL.R")
library(np)

################################################################################
# Parameters
################################################################################

## PMMR regularization para
PL <- -8-log(N/500,2)/2;      PU <- 1
BW.P.L <- -1;  BW.P.U <- 1
BW.T.L <- -1;  BW.T.U <- 1
Para.Grid <- expand.grid(seq(BW.P.L,BW.P.U,by=1),0,0,
                         seq(BW.T.L,BW.T.U,by=1),0,0,
                         seq(PL,PU,length=6)) 
Para.Grid[,2:3] <- Para.Grid[,1]
Para.Grid[,5:6] <- Para.Grid[,4]

## Superlearner parameters
SL.hpara <- list()
SL.hpara$SLL <- 1 # c(1,2,3,4,5,6,7,9,10)

# Superlearner basic learning algorithms:
# 1: GLM
# 2: lasso/ridge
# 3: earth
# 4: GAM
# 5: xgboost
# 6: polynomial spline
# 7: random forest
# 9: gbm
# 10: 1-layer MLP
SL.hpara$MLPL <- c(2,4)
SL.hpara$MTRY <- c(1,2)
SL.hpara$NMN <- min(N/50,25)
SL.hpara$MLPdecay <- 10^c(-1,-3)

## Cross-fitting/Cross-validation
CF       <- 2
NumCV    <- 5
NumCVRep <- 1


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
  
  P.P <- rep(as.numeric(PARAMETER[1]),3)
  P.T <- rep(as.numeric(PARAMETER[4]),3)
  Lambda  <- as.numeric(PARAMETER[7])
  
  RISK <- rep(0,NumCV)
  
  for(cv in 1:NumCV){
    
    CV.Split.Index <- list()
    CV.Split.Index[[1]] <- intersect( SS.Index[[ss]][ -CV.Index[[cv]] ],subset)
    CV.Split.Index[[2]] <- intersect( SS.Index[[ss]][ CV.Index[[cv]] ] ,subset)
    
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
    
    RISK[c(cv)] <- c(CV.result$Ustat)
  }
  
  return( mean(RISK) )
}

################################################################################
# Data Generation : Observational
################################################################################

if(Sim.Setup.Obs){
  
  #####################
  # DGP
  #####################
  
  betay0.1  <- 4
  betayx1.1 <- 0.35
  betayx2.1 <- 0.35
  betayw.1  <- 0.35
  betayu.1  <- 0.2
  
  betay0.0  <- 2
  betayx1.0 <- 0.35
  betayx2.0 <- 0.35
  betayw.0  <- 0.35
  betayu.0  <- 0.27
  
  betaz0  <- 0
  betazx1 <- 0.28
  betazx2 <- 0.28
  betazu  <- 0.45
  betaza  <- 0.53
  
  betaw0  <- 0
  betawx1 <- 0.17
  betawx2 <- 0.17
  betawu  <- 0.13
  
  betaa0  <- 0.05
  betaax1 <- 0.1
  betaax2 <- 0.1
  betaau  <- 0.07
  
  betad0.0  <- -0.5
  betadx1.0 <- 0.06
  betadx2.0 <- 0.06
  betadu.0  <- 0.02
  
  betad0.1  <- -0.7
  betadx1.1 <- 0.036
  betadx2.1 <- 0.036
  betadu.1  <- 0.015
  
  deltad0  <-  betad0.0  - betad0.1
  deltadx1 <-  betadx1.0 - betadx1.1
  deltadx2 <-  betadx2.0 - betadx2.1
  deltadu  <-  betadu.0  - betadu.1
  
  betaxu  <- 0.2
  
  #####################
  # Generate
  #####################
  
  U  <- rnorm(N)
  X1 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
  X2 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
  
  PrA1 <- expit( betaa0 + betaax1*X1 + betaax2*X2 + betaau*U )
  PrA0 <- 1-PrA1
  
  A <- rbinom(N,1,PrA1)
  
  Z  <- rnorm(N) + betaz0 + betazx1*X1 + betazx2*X2 + betazu*U + betaza*A
  W  <- rnorm(N) + betaw0 + betawx1*X1 + betawx2*X2 + betawu*U
  
  PrD1A1 <- 1-apply(cbind(0,exp(betad0.1+betadx1.1*X1+
                                  betadx2.1*X2+betadu.1*U),0.99),1,median)
  PrD1A0 <- 1-apply(cbind(0,exp(betad0.0+betadx1.0*X1+
                                  betadx2.0*X2+betadu.0*U),0.99),1,median)
  
  D.Ad1  <- rbinom(N,1,PrD1A1)      
  D.Ad0  <- rbinom(N,1,PrD1A0)      
  
  Y.Ay1.Ad1 <- rnorm( N, betay0.1 + betayx1.1*X1 + betayx2.1*X2 + 
                        betayw.1*W + betayu.1*U, 1  )*(1-D.Ad1)
  Y.Ay1.Ad0 <- rnorm( N, betay0.1 + betayx1.1*X1 + betayx2.1*X2 + 
                        betayw.1*W + betayu.1*U, 1  )*(1-D.Ad0)
  Y.Ay0.Ad1 <- rnorm( N, betay0.0 + betayx1.0*X1 + betayx2.0*X2 + 
                        betayw.0*W + betayu.0*U, 1  )*(1-D.Ad1)
  Y.Ay0.Ad0 <- rnorm( N, betay0.0 + betayx1.0*X1 + betayx2.0*X2 + 
                        betayw.0*W + betayu.0*U, 1  )*(1-D.Ad0)
  
  D <- A*(D.Ad1)+(1-A)*(D.Ad0)
  Y <- A*(Y.Ay1.Ad1)+(1-A)*(Y.Ay0.Ad0)
  
  Data <- data.frame( cbind(Y,A,D,W,Z,X1,X2,U) )
  
  #####################
  # True nuis fts
  #####################
  
  h11.true <- function(w,x1,x2){
    gamma0  <- (betay0.1 - betayu.1*betaw0/betawu)
    gammax1 <- (betayx1.1-betayu.1*betawx1/betawu)
    gammax2 <- (betayx2.1-betayu.1*betawx1/betawu)
    gammaw  <- (betayw.1+betayu.1/betawu)
    gamma0 + gammax1*x1 + gammax2*x2 + gammaw*w
  }
  
  h01.true <- function(w,x1,x2){
    gammaw  <- (betayu.1 + betayw.1*betawu)/betawu
    gammax1 <- betayx1.1 + betayw.1*betawx1 - gammaw*betawx1
    gammax2 <- betayx2.1 + betayw.1*betawx2 - gammaw*betawx2
    deltaw  <- betadu.0/betawu
    deltax1 <- betadx1.0 - deltaw*betawx1
    deltax2 <- betadx2.0 - deltaw*betawx2
    delta0  <- betad0.0 - deltaw*betaw0 - 0.5*deltaw^2
    gamma0  <- betay0.1 + betayw.1*betaw0 - gammaw*deltaw - gammaw*betaw0
    (gamma0+gammax1*X1+gammax2*X2+gammaw*W)*
      exp(delta0 + deltax1*X1 + deltax2*X2 + deltaw*W)
  }
  
  h00.true <- function(w,x1,x2){
    gammaw  <- (betayu.0 + betayw.0*betawu)/betawu
    gammax1 <- betayx1.0 + betayw.0*betawx1 - gammaw*betawx1
    gammax2 <- betayx2.0 + betayw.0*betawx2 - gammaw*betawx2
    deltaw  <- betadu.0/betawu
    deltax1 <- betadx1.0 - deltaw*betawx1
    deltax2 <- betadx2.0 - deltaw*betawx2
    delta0  <- betad0.0 - deltaw*betaw0 - 0.5*deltaw^2
    gamma0  <- betay0.0 + betayw.0*betaw0 - gammaw*deltaw - gammaw*betaw0
    (gamma0+gammax1*X1+gammax2*X2+gammaw*W)*
      exp(delta0 + deltax1*X1 + deltax2*X2 + deltaw*W)
  }
  
  h2.true <- function(w,x1,x2){
    deltaw  <- betadu.0/betawu
    deltax1 <- betadx1.0 - deltaw*betawx1
    deltax2 <- betadx2.0 - deltaw*betawx2
    delta0  <- betad0.0 - deltaw*betaw0 - 0.5*deltaw^2
    exp(delta0 + deltax1*X1 + deltax2*X2 + deltaw*W)
  }
  
  q0.true <- function(z,x1,x2){
    gammaz  <- betaau/betazu
    gammax1 <- betaax1 -gammaz*betazx1
    gammax2 <- betaax2 -gammaz*betazx2
    gamma0  <- betaa0 - 0.5*gammaz^2 - gammaz*betaz0
    1 + exp( gamma0 + gammax1*x1 + gammax2*x2 + gammaz*z )
  }
  
  q11.true <- function(z,x1,x2){
    gammaz  <- (deltadu - betaau)/betazu
    gammax1 <- deltadx1 - betaax1 - gammaz*betazx1
    gammax2 <- deltadx2 - betaax2 - gammaz*betazx2
    gamma0  <- deltad0  - betaa0  - 0.5*gammaz^2 - gammaz*betaz0 - gammaz*betaza
    
    sgammaz  <- (deltadu)/betazu
    sgammax1 <- deltadx1 - sgammaz*betazx1
    sgammax2 <- deltadx2 - sgammaz*betazx2
    sgamma0  <- deltad0  - 0.5*sgammaz^2   - sgammaz*betaz0 - sgammaz*betaza
    
    exp(gamma0 + gammax1*x1 + gammax2*x2 + gammaz*z) +
      exp(sgamma0 + sgammax1*x1 + sgammax2*x2 + sgammaz*z)
    
  }
  
}

################################################################################
# Data Generation : Experimental
################################################################################

if(!Sim.Setup.Obs){
  
  #####################
  # DGP
  #####################
  
  betay0.1  <- 4
  betayx1.1 <- 0.35
  betayx2.1 <- 0.35
  betayw.1  <- 0.35
  betayu.1  <- 0.2
  
  betay0.0  <- 2
  betayx1.0 <- 0.35
  betayx2.0 <- 0.35
  betayw.0  <- 0.35
  betayu.0  <- 0.27
  
  betaz0  <- 0
  betazx1 <- 0.28
  betazx2 <- 0.28
  betazu  <- 0.45
  
  betaw0  <- 0
  betawx1 <- 0.17
  betawx2 <- 0.17
  betawu  <- 0.13
  
  betaa0  <- 0.05
  betaax1 <- 0.1
  betaax2 <- 0.1
  betaau  <- 0.07
  
  betad0.0  <- -0.5
  betadx1.0 <- 0.03
  betadx2.0 <- 0.03
  betadu.0  <- 0.01
  betadz.0  <- 0.032
  
  betad0.1  <- -0.7
  betadx1.1 <- 0.014
  betadx2.1 <- 0.014
  betadu.1  <- 0.007
  betadz.1  <- 0.029
  
  
  deltad0  <-  betad0.0  - betad0.1
  deltadx1 <-  betadx1.0 - betadx1.1
  deltadx2 <-  betadx2.0 - betadx2.1
  deltadu  <-  betadu.0  - betadu.1
  
  betaxu  <- 0.2
  
  #####################
  # Generate
  #####################
  
  U  <- rnorm(N)
  X1 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
  X2 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
  Z  <- rnorm(N) + betaz0 + betazx1*X1 + betazx2*X2 + betazu*U
  W  <- rnorm(N) + betaw0 + betawx1*X1 + betawx2*X2 + betawu*U
  
  PrD1A1 <- 1-apply(cbind(0,exp(betad0.1+betadx1.1*X1+betadx2.1*X2 + 
                                  betadz.1*Z +betadu.1*U),0.99),1,median)
  PrD1A0 <- 1-apply(cbind(0,exp(betad0.0+betadx1.0*X1+betadx2.0*X2 + 
                                  betadz.0*Z +betadu.0*U),0.99),1,median)
  
  D.Ad1  <- rbinom(N,1,PrD1A1)
  D.Ad0  <- rbinom(N,1,PrD1A0)
  
  Y.Ay1.Ad1 <- rnorm( N, betay0.1 + betayx1.1*X1 + betayx2.1*X2 + 
                        betayw.1*W + betayu.1*U, 1  )*(1-D.Ad1)
  Y.Ay1.Ad0 <- rnorm( N, betay0.1 + betayx1.1*X1 + betayx2.1*X2 + 
                        betayw.1*W + betayu.1*U, 1  )*(1-D.Ad0)
  Y.Ay0.Ad1 <- rnorm( N, betay0.0 + betayx1.0*X1 + betayx2.0*X2 + 
                        betayw.0*W + betayu.0*U, 1  )*(1-D.Ad1)
  Y.Ay0.Ad0 <- rnorm( N, betay0.0 + betayx1.0*X1 + betayx2.0*X2 + 
                        betayw.0*W + betayu.0*U, 1  )*(1-D.Ad0)
  
  PrA1 <- 0.6
  PrA0 <- 1-PrA1
  
  A <- rbinom(N,1,PrA1)
  D <- A*(D.Ad1)+(1-A)*(D.Ad0)
  Y <- A*(Y.Ay1.Ad1)+(1-A)*(Y.Ay0.Ad0)
  
  Data <- data.frame( cbind(Y,A,D,W,Z,X1,X2,U) )
  
  #####################
  # True nuis fts
  #####################
  
  MeanUXWZ  <- matrix(c(0,0,0,betaw0,betaz0),5,1)
  SigmaUXWZ <- matrix(0,5,5)
  SigmaUXWZ[1,1] <- 1
  SigmaUXWZ[1,2] <- SigmaUXWZ[2,1] <- betaxu
  SigmaUXWZ[1,3] <- SigmaUXWZ[3,1] <- betaxu
  SigmaUXWZ[2,3] <- SigmaUXWZ[3,2] <- betaxu^2
  SigmaUXWZ[2,2] <- SigmaUXWZ[3,3] <- 1
  
  SigmaUXWZ[1:3,4] <- SigmaUXWZ[1:3,1:3]%*%c(betawu,betawx1,betawx2)
  SigmaUXWZ[1:3,5] <- SigmaUXWZ[1:3,1:3]%*%c(betazu,betazx1,betazx2)
  SigmaUXWZ[4,4] <- 1 + t(c(betawu,betawx1,betawx2))%*%
    SigmaUXWZ[1:3,1:3]%*%(c(betawu,betawx1,betawx2))
  
  SigmaUXWZ[4,1:3] <- t(SigmaUXWZ[1:3,4])
  SigmaUXWZ[5,1:3] <- t(SigmaUXWZ[1:3,5])
  SigmaUXWZ[5,5] <- 1 + t(c(betazu,betazx1,betazx2))%*%
    SigmaUXWZ[1:3,1:3]%*%(c(betazu,betazx1,betazx2))
  
  SigmaUXWZ[5,4] <- SigmaUXWZ[4,5] <- 
    (SigmaUXWZ[4:5,1:3]%*%solve(SigmaUXWZ[1:3,1:3])%*%SigmaUXWZ[1:3,4:5])[1,2]
  
  
  MGF <- function(tt,mu.vec,Sigma.mat){
    exp(sum(mu.vec*tt) + 0.5*t(tt)%*%Sigma.mat%*%(tt))
  }
  
  COND.ZU.XW.Mu <- function(w,x1,x2){
    MeanUXWZ[c(5,1)] + SigmaUXWZ[c(5,1),c(2,3,4)]%*%
      solve(SigmaUXWZ[c(2,3,4),c(2,3,4)])%*%( c(x1,x2,w) - MeanUXWZ[c(2,3,4)])
  }
  COND.ZU.XW.Var <- SigmaUXWZ[c(5,1),c(5,1)] - SigmaUXWZ[c(5,1),c(2,3,4)]%*%
    solve(SigmaUXWZ[c(2,3,4),c(2,3,4)])%*%SigmaUXWZ[c(2,3,4),c(5,1)]
  
  Smat <- solve(SigmaUXWZ[1:4,1:4])
  
  
  betau0.00wx <- Smat[1,1]^(-1)*betaw0*Smat[4,1] + Smat[1,1]*betadu.0
  betaux1.00wx <- -Smat[1,1]^(-1)*Smat[2,1]
  betaux2.00wx <- -Smat[1,1]^(-1)*Smat[3,1]
  betauw.00wx <- -Smat[1,1]^(-1)*Smat[4,1]
  
  
  h11.true <- function(w,x1,x2){
    gamma0  <- (betay0.1 - betayu.1*betaw0/betawu)
    gammax1 <- (betayx1.1-betayu.1*betawx1/betawu)
    gammax2 <- (betayx2.1-betayu.1*betawx1/betawu)
    gammaw  <- (betayw.1+betayu.1/betawu)
    gamma0 + gammax1*x1 + gammax2*x2 + gammaw*w
  }
  
  p0a0wx.true <- function(w,x1,x2){
    apply(cbind(w,x1,x2),
          1,
          function(v){
            exp(betad0.0+betadx1.0*v[2]+betadx2.0*v[3])*
              MGF(c(betadz.0,betadu.0),
                  COND.ZU.XW.Mu(v[1],v[2],v[3]),
                  COND.ZU.XW.Var )
          }
    )
    
  } 
  
  h01.true <- function(w,x1,x2){
    gamma0  <- betay0.1 - betayu.1*betaw0/betawu
    gammax1 <- betayx1.1 - betayu.1*betawx1/betawu
    gammax2 <- betayx2.1 - betayu.1*betawx2/betawu
    gammaw  <- betayw.1 + betayu.1/betawu
    
    p0a0wx.true(w,x1,x2)*(gamma0 + gammax1*x1 + gammax2*x2 + gammaw*w)
  }
  
  h00.true <- function(w,x1,x2){
    gamma0  <- betay0.0
    gammax1 <- betayx1.0
    gammax2 <- betayx2.0
    gammaw  <- betayw.0
    gammau0   <- betayu.0*(betau0.00wx)
    gammaux1  <- betayu.0*(betaux1.00wx)
    gammaux2  <- betayu.0*(betaux2.00wx)
    gammauw   <- betayu.0*(betauw.00wx)
    
    p0a0wx.true(w,x1,x2)*(gamma0 + gammax1*x1 + gammax2*x2 + gammaw*w + 
                            gammau0 + gammaux1*x1 + gammaux2*x2 + gammauw*w)
  }
  
  h2.true <- function(w,x1,x2){
    p0a0wx.true(w,x1,x2)
  }
  
  q0.true <- function(z,x1,x2){
    rep(1/PrA0,length(z))
  }
  
  q11.true <- function(z,x1,x2){
    gammaz  <- (deltadu)/betazu
    gammax1 <- deltadx1 - gammaz*betazx1
    gammax2 <- deltadx2 - gammaz*betazx2
    gamma0  <- deltad0 - log(PrA1)  - 0.5*gammaz^2   - gammaz*betaz0
    exp(gamma0 + gammax1*x1 + gammax2*x2 + gammaz*z)
  }
}

################################################################################
# Split Samples/Cross Validation Sets
################################################################################

set.A0D0 <- which(A==0 & D==0)
set.A0D1 <- which(A==0 & D==1)
set.A1D0 <- which(A==1 & D==0)
set.A1D1 <- which(A==1 & D==1)

Random.A0D0 <- sample( (set.A0D0),
                       length(set.A0D0) )
Random.A0D1 <- sample( (set.A0D1),
                       length(set.A0D1) )
Random.A1D0 <- sample( (set.A1D0),
                       length(set.A1D0) )
Random.A1D1 <- sample( (set.A1D1),
                       length(set.A1D1) )
q.A0D0 <- round(quantile(1:(length(Random.A0D0)+1),seq(0,1,length=11)))
q.A0D1 <- round(quantile(1:(length(Random.A0D1)+1),seq(0,1,length=11)))
q.A1D0 <- round(quantile(1:(length(Random.A1D0)+1),seq(0,1,length=11)))
q.A1D1 <- round(quantile(1:(length(Random.A1D1)+1),seq(0,1,length=11)))
Suffle <- NULL
for(tt in 1:10){
  Suffle <- c(Suffle,
              Random.A0D0[q.A0D0[tt]:(q.A0D0[tt+1]-1)],
              Random.A0D1[q.A0D1[tt]:(q.A0D1[tt+1]-1)],
              Random.A1D0[q.A1D0[tt]:(q.A1D0[tt+1]-1)],
              Random.A1D1[q.A1D1[tt]:(q.A1D1[tt+1]-1)])
}
Data <- Data[Suffle,]

SS.Index <- list()
for(ss in 1:CF){
  SS.Index[[ss]] <- ((seq(0,N,length=CF+1))[ss]+1):(seq(0,N,length=CF+1)[ss+1])
} 
CV.Index <- list()
CV.CUT <- seq(0,round(N/CF),length=NumCV+1)
for(cv in 1:NumCV){
  CV.Index[[cv]] <- ((CV.CUT)[cv]+1):(CV.CUT[cv+1])
}

if(NumCVRep>1){
  
  REORDER <- NULL
  for(cv in 1:NumCV){
    REORDER <- rbind(REORDER,((CV.CUT)[cv]+1):(CV.CUT[cv+1]))
  }
  for(tt in 2:NumCVRep){
    SFF <- sapply(1:(N/2/NumCV),function(t){sample(5,5)})
    for(cv in 1:NumCV){
      CV.Index[[cv + (tt-1)*NumCV]] <- 
        sapply(1:(N/2/NumCV),function(t){REORDER[SFF[cv,t],t]})
    }  
  }
}

X.pos <- which( substr(colnames(Data),1,1)=="X" )
W.pos <- which( substr(colnames(Data),1,1)=="W" )
Z.pos <- which( substr(colnames(Data),1,1)=="Z" )
A.pos <- which( substr(colnames(Data),1,1)=="A" )
U.pos <- which( substr(colnames(Data),1,1)=="U" )
Y.pos <- which( substr(colnames(Data),1,1)=="Y" )
D.pos <- which( substr(colnames(Data),1,1)=="D" )

Y <- Data[,Y.pos]
A <- Data[,A.pos]
D <- Data[,D.pos]
W <- Data[,W.pos]
Z <- Data[,Z.pos]
X <- Data[,X.pos]
X1 <- X[,1]
X2 <- X[,2]
U <- Data[,U.pos]
WX <- cbind(W,X)
ZX <- cbind(Z,X)
rY <- range(Y)

WX.bw <- rep(1,3)
ZX.bw <- rep(1,3)

Y.MM  <- list()
A.MM  <- list()
D.MM  <- list()
X.MM  <- list()
X1.MM  <- list()
X2.MM  <- list()
W.MM  <- list()
Z.MM  <- list()
WX.MM <- list()
ZX.MM <- list()


for(ss in 1:CF){
  
  Y.MM [[ss]]  <- Y[SS.Index[[ss]] ]
  A.MM [[ss]]  <- A[SS.Index[[ss]] ]
  D.MM [[ss]]  <- D[SS.Index[[ss]] ]
  X.MM [[ss]]  <- X[SS.Index[[ss]],]
  X1.MM [[ss]]  <- X[SS.Index[[ss]],1]
  X2.MM [[ss]]  <- X[SS.Index[[ss]],2]
  W.MM [[ss]]  <- W[SS.Index[[ss]] ]
  Z.MM [[ss]]  <- Z[SS.Index[[ss]] ]
  WX.MM[[ss]]  <- cbind(W.MM[[ss]],X.MM[[ss]])
  ZX.MM[[ss]]  <- cbind(Z.MM[[ss]],X.MM[[ss]])
  
}

posA1D0 <- which(A==1&D==0)
posA0D0 <- which(A==0&D==0)
posA0   <- which(A==0)
posD0   <- which(D==0)

SS.posA1D0 <- SS.posA0D0 <- SS.posA0 <- SS.posD0 <- list()
for(ss in 1:CF){
  SS.posA1D0[[ss]] <- which(A.MM[[ss]]==1&D.MM[[ss]]==0)
  SS.posA0D0[[ss]] <- which(A.MM[[ss]]==0&D.MM[[ss]]==0)
  SS.posA0[[ss]]   <- which(A.MM[[ss]]==0)
  SS.posD0[[ss]]   <- which(D.MM[[ss]]==0)
}

################################################################################
# ML for Experimental Settings: Pr(D=0|A=a_D,W,X)
################################################################################

if(!Sim.Setup.Obs){
  ModelD1gA0WX <- list()
  for(ss in 1:CF){
    ModelD1gA0WX[[ss]] <- 
      MySL(Data[ intersect(SS.Index[[ss]],which(A==0)), ],
           locY=D.pos,
           locX=c(W.pos,X.pos),
           Ydist=binomial(),
           SL.list=SL.hpara$SLL,
           MTRY=SL.hpara$MTRY,
           MLPL=SL.hpara$MLPL,
           NMN=SL.hpara$NMN,
           MLPdecay=SL.hpara$MLPdecay)
  }
  
  D0gA0WX.NoCF <- rep(0,N)
  for(ss in 1:CF){
    D0gA0WX.NoCF[ SS.Index[[ss]] ] <- 
      1-predict(ModelD1gA0WX[[ss]], 
                newdata=Data[ SS.Index[[ss]] , c(W.pos,X.pos)])$pred
  }
  
  D0gA0WX.CF <- rep(0,N)
  for(ss in 1:CF){
    D0gA0WX.CF[ SS.Index[[3-ss]] ] <- 
      1-predict(ModelD1gA0WX[[ss]], 
                newdata=Data[ SS.Index[[3-ss]] , c(W.pos,X.pos)])$pred
  }
  
  D0gA0WX.NoCF.MM <- list()
  for(ss in 1:CF){
    D0gA0WX.NoCF.MM[[ss]]  <- D0gA0WX.NoCF[ SS.Index[[ss]] ]
  }
  
}

################################################################################
# Cross Validation and Estimation
################################################################################

###################################################
# Bridge Function: h11
###################################################

Opt.Para.h11 <- list()
for(ss in 1:CF){
  
  CV.Curve <- apply(Para.Grid,
                    1,
                    FUN=function(vv){ CrossVar(vv,
                                               ss,
                                               subset=posA1D0,
                                               response=Y*(1-D)*A,
                                               target=WX,
                                               perturb=ZX,
                                               diagonal=(1-D)*A) })
  
  Opt.Para.h11[[ss]] <-
    as.numeric(Para.Grid[which.min(CV.Curve),])
}

pos.h11.MM <- list()
h11.MM <- list()
h11.predict <- list()

for(ss in 1:CF){
  
  pos.h11.MM[[ss]] <- SS.posA1D0[[ss]]
  
  h11.MM[[ss]] <- 
    FT_PMMR( Y         =(Y.MM[[ss]]*(A.MM[[ss]])*(1-D.MM[[ss]]))[pos.h11.MM[[ss]]],
             Perturb   =ZX.MM[[ss]][pos.h11.MM[[ss]],],
             Target    =WX.MM[[ss]][pos.h11.MM[[ss]],],
             Diagonal  =((A.MM[[ss]])*(1-D.MM[[ss]]))[pos.h11.MM[[ss]]],
             Perturb.bw=exp(Opt.Para.h11[[ss]][1:3]),
             Target.bw =exp(Opt.Para.h11[[ss]][4:6]),
             lambda    =exp(Opt.Para.h11[[ss]][7]),
             NV        =FALSE)
}

h11.predict[[1]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[1]][pos.h11.MM[[1]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h11[[1]][4:6]))%*%h11.MM[[1]]$alpha + h11.MM[[1]]$intercept
}

h11.predict[[2]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[2]][pos.h11.MM[[2]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h11[[2]][4:6]))%*%h11.MM[[2]]$alpha + h11.MM[[2]]$intercept
}

h11hat.NoCF <- h11hat.CF <- rep(0,N)
for(ss in 1:CF){
  h11hat.NoCF[SS.Index[[ss]]] <- h11.predict[[ss]](WX.MM[[ss]])
  h11hat.CF[SS.Index[[3-ss]]] <- h11.predict[[ss]](WX.MM[[3-ss]])
}

h11hat.NoCF.MM <- list()
for(ss in 1:CF){
  print(ss)
  h11hat.NoCF.MM[[ss]]  <- h11hat.NoCF[SS.Index[[ss]] ]
}

###################################################
# Bridge Function: h01
###################################################

if(Sim.Setup.Obs){
  
  Opt.Para.h01 <- list()
  for(ss in 1:CF){
    
    
    CV.Curve <- apply(Para.Grid,
                      1,
                      FUN=function(vv){ CrossVar(vv,
                                                 ss,
                                                 subset=posA0,
                                                 response=h11hat.NoCF*(1-D)*(1-A),
                                                 target=WX,
                                                 perturb=ZX,
                                                 diagonal=(1-A)) })
    
    Opt.Para.h01[[ss]] <-
      as.numeric(Para.Grid[which.min(CV.Curve),])
  }
  
  
  pos.h01.MM <- list()
  h01.MM <- list()
  h01.predict <- list()
  
  for(ss in 1:CF){
    
    pos.h01.MM[[ss]] <- SS.posA0[[ss]]
    Resp.Temp <- h11hat.NoCF.MM[[ss]]*(1-A.MM[[ss]])*(1-D.MM[[ss]])
    
    h01.MM[[ss]] <- 
      FT_PMMR( Y         =Resp.Temp[pos.h01.MM[[ss]]],
               Perturb   =ZX.MM[[ss]][pos.h01.MM[[ss]],],
               Target    =WX.MM[[ss]][pos.h01.MM[[ss]],],
               Diagonal  =(1-A.MM[[ss]])[pos.h01.MM[[ss]]],
               Perturb.bw=exp(Opt.Para.h01[[ss]][1:3]),
               Target.bw =exp(Opt.Para.h01[[ss]][4:6]),
               lambda    =exp(Opt.Para.h01[[ss]][7]),
               NV        =FALSE)
    
  }
  
  
  h01.predict[[1]] <- function(WX.New.Input){
    FT_RBF(X     = WX.MM[[1]][pos.h01.MM[[1]],],
           X.new = WX.New.Input,
           bw.median = exp(Opt.Para.h01[[1]][4:6]))%*%h01.MM[[1]]$alpha + h01.MM[[1]]$intercept
  }
  
  h01.predict[[2]] <- function(WX.New.Input){
    FT_RBF(X     = WX.MM[[2]][pos.h01.MM[[2]],],
           X.new = WX.New.Input,
           bw.median = exp(Opt.Para.h01[[2]][4:6]))%*%h01.MM[[2]]$alpha + h01.MM[[2]]$intercept
  }
  
  h01hat.NoCF <- h01hat.CF <- rep(0,N)
  for(ss in 1:CF){
    h01hat.NoCF[SS.Index[[ss]]] <- h01.predict[[ss]](WX.MM[[ss]])
    h01hat.CF[SS.Index[[3-ss]]] <- h01.predict[[ss]](WX.MM[[3-ss]])
  }
  
  h01hat.NoCF.MM <- list()
  for(ss in 1:CF){
    print(ss)
    h01hat.NoCF.MM[[ss]]  <- h01hat.NoCF[SS.Index[[ss]] ]
  }
} 

if(!Sim.Setup.Obs){
  h01hat.NoCF <- h11hat.NoCF*D0gA0WX.NoCF
  h01hat.CF   <- h11hat.CF*D0gA0WX.CF
}

###################################################
# Bridge Function: h00
###################################################

if(Sim.Setup.Obs){
  
  Opt.Para.h00 <- list()
  for(ss in 1:CF){
    
    CV.Curve <- apply(Para.Grid,
                      1,
                      FUN=function(vv){ CrossVar(vv,
                                                 ss,
                                                 subset=posA0,
                                                 response=Y*(1-D)*(1-A),
                                                 target=WX,
                                                 perturb=ZX,
                                                 diagonal=(1-A)) })
    
    Opt.Para.h00[[ss]] <-
      as.numeric(Para.Grid[which.min(CV.Curve),])
  }
  
  pos.h00.MM <- list()
  h00.MM <- list()
  h00.predict <- list()
  
  for(ss in 1:CF){
    
    pos.h00.MM[[ss]] <- SS.posA0[[ss]]
    Resp.Temp <- Y.MM[[ss]]*(1-A.MM[[ss]])*(1-D.MM[[ss]])
    
    h00.MM[[ss]] <- FT_PMMR( Y         =Resp.Temp[pos.h00.MM[[ss]]],
                             Perturb   =ZX.MM[[ss]][pos.h00.MM[[ss]],],
                             Target    =WX.MM[[ss]][pos.h00.MM[[ss]],],
                             Diagonal  =(1-A.MM[[ss]])[pos.h00.MM[[ss]]], 
                             Perturb.bw=exp(Opt.Para.h00[[ss]][1:3]),
                             Target.bw =exp(Opt.Para.h00[[ss]][4:6]),
                             lambda    =exp(Opt.Para.h00[[ss]][7]),
                             NV        =FALSE)
    
  }
  
  
  h00.predict[[1]] <- function(WX.New.Input){
    FT_RBF(X     = WX.MM[[1]][pos.h00.MM[[1]],],
           X.new = WX.New.Input,
           bw.median = exp(Opt.Para.h00[[1]][4:6]))%*%h00.MM[[1]]$alpha + h00.MM[[1]]$intercept
  }
  
  h00.predict[[2]] <- function(WX.New.Input){
    FT_RBF(X     = WX.MM[[2]][pos.h00.MM[[2]],],
           X.new = WX.New.Input,
           bw.median = exp(Opt.Para.h00[[2]][4:6]))%*%h00.MM[[2]]$alpha + h00.MM[[2]]$intercept
  }
  
  h00hat.NoCF <- h00hat.CF <- rep(0,N)
  for(ss in 1:CF){
    h00hat.NoCF[SS.Index[[ss]]] <- h00.predict[[ss]](WX.MM[[ss]])
    h00hat.CF[SS.Index[[3-ss]]] <- h00.predict[[ss]](WX.MM[[3-ss]])
  }
  
  h00hat.NoCF.MM <- list()
  for(ss in 1:CF){
    print(ss)
    h00hat.NoCF.MM[[ss]]  <- h00hat.NoCF[SS.Index[[ss]] ]
  }
}



###################################################
# Bridge Function: h2
###################################################

if(Sim.Setup.Obs){
  
  Opt.Para.h2 <- list()
  for(ss in 1:CF){
    
    CV.Curve <- apply(Para.Grid,
                      1,
                      FUN=function(vv){ CrossVar(vv,
                                                 ss,
                                                 subset=posA0,
                                                 response=(1-D)*(1-A),
                                                 target=WX,
                                                 perturb=ZX,
                                                 diagonal=(1-A)) })
    
    Opt.Para.h2[[ss]] <-
      as.numeric(Para.Grid[which.min(CV.Curve),])
  }
  
  pos.h2.MM <- list()
  h2.MM <- list()
  h2.predict <- list()
  
  for(ss in 1:CF){
    
    pos.h2.MM[[ss]] <- SS.posA0[[ss]]
    Resp.Temp <- (1-D.MM[[ss]])*(1-A.MM[[ss]])
    
    h2.MM[[ss]] <- FT_PMMR( Y         =Resp.Temp[pos.h2.MM[[ss]]],
                            Perturb   =ZX.MM[[ss]][pos.h2.MM[[ss]],],
                            Target    =WX.MM[[ss]][pos.h2.MM[[ss]],],
                            Diagonal  =(1-A.MM[[ss]])[pos.h2.MM[[ss]]], 
                            Perturb.bw=exp(Opt.Para.h2[[ss]][1:3]),
                            Target.bw =exp(Opt.Para.h2[[ss]][4:6]),
                            lambda    =exp(Opt.Para.h2[[ss]][7]),
                            NV        =FALSE)
    
  }
  
  h2.predict[[1]] <- function(WX.New.Input){
    FT_RBF(X     = WX.MM[[1]][pos.h2.MM[[1]],],
           X.new = WX.New.Input,
           bw.median = exp(Opt.Para.h2[[1]][4:6]))%*%h2.MM[[1]]$alpha + h2.MM[[1]]$intercept
  }
  
  h2.predict[[2]] <- function(WX.New.Input){
    FT_RBF(X     = WX.MM[[2]][pos.h2.MM[[2]],],
           X.new = WX.New.Input,
           bw.median = exp(Opt.Para.h2[[2]][4:6]))%*%h2.MM[[2]]$alpha + h2.MM[[2]]$intercept
  }
  
  h2hat.NoCF <- h2hat.CF <- rep(0,N)
  for(ss in 1:CF){
    h2hat.NoCF[SS.Index[[ss]]] <- h2.predict[[ss]](WX.MM[[ss]])
    h2hat.CF[SS.Index[[3-ss]]] <- h2.predict[[ss]](WX.MM[[3-ss]])
  }
  
  h2hat.NoCF.MM <- list()
  for(ss in 1:CF){
    print(ss)
    h2hat.NoCF.MM[[ss]]  <- h2hat.NoCF[SS.Index[[ss]] ]
  }
}

###################################################
# Bridge Function: q0
###################################################

if(Sim.Setup.Obs){
  
  Opt.Para.q0 <- list()
  for(ss in 1:CF){
    
    CV.Curve <- apply(Para.Grid,
                      1,
                      FUN=function(vv){ CrossVar(vv,
                                                 ss,
                                                 subset=(1:N),
                                                 response=rep(1,N),
                                                 target=ZX,
                                                 perturb=WX,
                                                 diagonal=(1-A)) })
    
    Opt.Para.q0[[ss]] <-
      as.numeric(Para.Grid[which.min(CV.Curve),])
  }
  
  pos.q0.MM <- list()
  q0.MM <- list()
  q0.predict <- list()
  
  for(ss in 1:CF){
    
    pos.q0.MM[[ss]] <- 1:(N/CF)
    
    Resp.Temp <- rep(1,N)
    
    q0.MM[[ss]] <- FT_PMMR( Y         =Resp.Temp[pos.q0.MM[[ss]]],
                            Perturb   =WX.MM[[ss]][pos.q0.MM[[ss]],],
                            Target    =ZX.MM[[ss]][pos.q0.MM[[ss]],],
                            Diagonal  =(1-A.MM[[ss]])[pos.q0.MM[[ss]]],
                            Perturb.bw=exp(Opt.Para.q0[[ss]][1:3]),
                            Target.bw =exp(Opt.Para.q0[[ss]][4:6]),
                            lambda    =exp(Opt.Para.q0[[ss]][7]),
                            NV        =FALSE)
  }
  
  
  q0.predict[[1]] <- function(ZX.New.Input){
    FT_RBF(X     = ZX.MM[[1]][pos.q0.MM[[1]],],
           X.new = ZX.New.Input,
           bw.median = exp(Opt.Para.q0[[1]][4:6]))%*%q0.MM[[1]]$alpha + q0.MM[[1]]$intercept 
  }
  
  q0.predict[[2]] <- function(ZX.New.Input){
    FT_RBF(X     = ZX.MM[[2]][pos.q0.MM[[2]],],
           X.new = ZX.New.Input,
           bw.median = exp(Opt.Para.q0[[2]][4:6]))%*%q0.MM[[2]]$alpha + q0.MM[[2]]$intercept
  }
  
  q0hat.NoCF <- q0hat.CF <- rep(0,N)
  for(ss in 1:CF){
    q0hat.NoCF[SS.Index[[ss]]] <- q0.predict[[ss]](ZX.MM[[ss]])
    q0hat.CF[SS.Index[[3-ss]]] <- q0.predict[[ss]](ZX.MM[[3-ss]])
  }
  
  q0hat.NoCF.MM <- list()
  for(ss in 1:CF){
    q0hat.NoCF.MM[[ss]] <- q0hat.NoCF[SS.Index[[ss]] ]
  }
}

if(!Sim.Setup.Obs){
  q0hat.NoCF <- q0hat.CF <- rep(1/PrA0,N)
  
  q0hat.NoCF.MM <- list()
  for(ss in 1:CF){
    print(ss)
    q0hat.NoCF.MM[[ss]] <- q0hat.NoCF[SS.Index[[ss]] ]
  }
}

###################################################
# Bridge Function: q11
###################################################

Opt.Para.q11 <- list()
for(ss in 1:CF){
  
  CV.Curve <- apply(Para.Grid,
                    1,
                    FUN=function(vv){ CrossVar(vv,
                                               ss,
                                               subset = posD0,
                                               response = (1-A)*q0hat.NoCF*(1-D),
                                               target=ZX,
                                               perturb=WX,
                                               diagonal = A*(1-D)) })
  
  Opt.Para.q11[[ss]] <-
    as.numeric(Para.Grid[which.min(CV.Curve),])
}

pos.q11.MM <- list()
q11.MM <- list()
q11.predict <- list()

for(ss in 1:CF){
  
  pos.q11.MM[[ss]] <- SS.posD0[[ss]]
  
  Resp.Temp <- (1-A.MM[[ss]])*(1-D.MM[[ss]])*(q0hat.NoCF.MM[[ss]])
  Diag.Temp <- A.MM[[ss]]*(1-D.MM[[ss]])
  
  q11.MM[[ss]] <- FT_PMMR( Y         =Resp.Temp[pos.q11.MM[[ss]]],
                           Perturb   =WX.MM[[ss]][pos.q11.MM[[ss]],],
                           Target    =ZX.MM[[ss]][pos.q11.MM[[ss]],],
                           Diagonal  =Diag.Temp[pos.q11.MM[[ss]]],
                           Perturb.bw=exp(Opt.Para.q11[[ss]][1:3]),
                           Target.bw =exp(Opt.Para.q11[[ss]][4:6]),
                           lambda    =exp(Opt.Para.q11[[ss]][7]),
                           NV        =FALSE)
}

q11.predict[[1]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX.MM[[1]][pos.q11.MM[[1]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q11[[1]][4:6]))%*%q11.MM[[1]]$alpha + q11.MM[[1]]$intercept
}

q11.predict[[2]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX.MM[[2]][pos.q11.MM[[2]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q11[[2]][4:6]))%*%q11.MM[[2]]$alpha + q11.MM[[2]]$intercept
}

q11hat.NoCF <- q11hat.CF <- rep(0,N)
for(ss in 1:CF){
  q11hat.NoCF[SS.Index[[ss]]] <- q11.predict[[ss]](ZX.MM[[ss]])
  q11hat.CF[SS.Index[[3-ss]]] <- q11.predict[[ss]](ZX.MM[[3-ss]])
}

q11hat.NoCF.MM <- list()
for(ss in 1:CF){
  print(ss)
  q11hat.NoCF.MM[[ss]]  <- q11hat.NoCF[SS.Index[[ss]] ]
}

################################################################################
# Estimand
################################################################################

if(Sim.Setup.Obs){
  IF.Numer.1 <- (A)*(1-D)*q11hat.CF*(Y-h11hat.CF)+(1-A)*q0hat.CF*((1-D)*h11hat.CF-h01hat.CF)+h01hat.CF
  IF.Numer.0 <- (1-A)*q0hat.CF*((1-D)*Y-h00hat.CF)+h00hat.CF
  IF.Denom   <- (1-A)*q0hat.CF*((1-D)-h2hat.CF)+h2hat.CF
  
  IF.Numer.1.True <-   
    (A)*(1-D)*q11.true(Z,X1,X2)*(Y-h11.true(W,X1,X2))+
    (1-A)*q0.true(Z,X1,X2)*((1-D)*h11.true(W,X1,X2)-h01.true(W,X1,X2))+h01.true(W,X1,X2)
  IF.Numer.0.True <- 
    (1-A)*q0.true(Z,X1,X2)*((1-D)*Y-h00.true(W,X1,X2))+h00.true(W,X1,X2)
  IF.Denom.True   <- 
    (1-A)*q0.true(Z,X1,X2)*((1-D)-h2.true(W,X1,X2))+h2.true(W,X1,X2)
  
  IF <- cbind(IF.Numer.1,
              IF.Numer.0,
              IF.Numer.1-IF.Numer.0,
              IF.Denom)
  Sigma <- as.matrix( var(IF)/N )
  
  Est.Numer.1 <- as.numeric( apply(IF,2,mean)[1] )
  Est.Numer.0 <- as.numeric( apply(IF,2,mean)[2] )
  Est.Numer   <- as.numeric( apply(IF,2,mean)[3] )
  Est.Denom   <- as.numeric( apply(IF,2,mean)[4] )
  Est.1       <- Est.Numer.1/Est.Denom
  Est.0       <- Est.Numer.0/Est.Denom
  Est         <- Est.Numer/Est.Denom
  
  IF.True <- cbind(IF.Numer.1.True,
                   IF.Numer.0.True,
                   IF.Numer.1.True-IF.Numer.0.True,
                   IF.Denom.True)
  Sigma.True <- as.matrix( var(IF.True)/N )
  
  Est.Numer.1.True <- as.numeric( apply(IF.True,2,mean)[1] )
  Est.Numer.0.True <- as.numeric( apply(IF.True,2,mean)[2] )
  Est.Numer.True   <- as.numeric( apply(IF.True,2,mean)[3] )
  Est.Denom.True   <- as.numeric( apply(IF.True,2,mean)[4] )
  Est.1.True       <- Est.Numer.1.True/Est.Denom.True
  Est.0.True       <- Est.Numer.0.True/Est.Denom.True
  Est.True         <- Est.Numer.True/Est.Denom.True
  
  contrast   <- matrix(c(0,0,1/Est.Denom,-Est/Est.Denom),1,4)
  SE         <- sqrt( contrast%*%Sigma%*%t(contrast) )
  
  contrast.True   <- matrix(c(0,0,1/Est.Denom.True,-Est.True/Est.Denom.True),1,4)
  SE.True    <- sqrt( contrast.True%*%Sigma.True%*%t(contrast.True) )
  
  IF.Est      <- (IF[,3]-Est*IF[,4])/Est.Denom
  IF.Est.True <- (IF.True[,3]-Est.True*IF.True[,4])/Est.Denom.True
  
  BOOT <- BOOT.True <- BOOT.GMM <- rep(0,10000)
  for(bb in 1:10000){
    vv <- mean( rnorm(N)*IF.Est )
    BOOT[bb] <- vv
    vv <- mean( rnorm(N)*IF.Est.True )
    BOOT.True[bb] <- vv
  }
  
  Result1 <- c(Est, SE, sd(BOOT), 
               Est.Numer.1, Est.Numer.0, Est.Denom,
               
               Est.True, SE.True, sd(BOOT.True), 
               Est.Numer.1.True, Est.Numer.0.True, Est.Denom.True,
               
               c( sd(h11hat.CF), sd(h01hat.CF), sd(h00hat.CF), 
                  sd(h2hat.CF), sd(q0hat.CF), sd(q11hat.CF) ) , 
               c( sd(h11.true(W,X1,X2)), sd(h01.true(W,X1,X2)), 
                  sd(h00.true(W,X1,X2)), sd(h2.true(W,X1,X2)), 
                  sd(q0.true(Z,X1,X2)), sd(q11.true(Z,X1,X2)) )
               
  )
  
  Result1 <- matrix(Result1,1,length(Result1))
  
  colnames(Result1) <- c("Est",    "SE",    "SE.B",
                         "Est.Numer.1", "Est.Numer.0", "Est.Denom",
                         
                         "Est.True",    "SE.True",    "SE.B.True",
                         "Est.Numer.1.True", "Est.Numer.0.True", "Est.Denom.True",
                         
                         "SE.h11.PMMR", "SE.h01.PMMR", "SE.h00.PMMR", 
                         "SE.h2.PMMR", "SE.q0.PMMR", "SE.q11.PMMR",
                         "SE.h11.True", "SE.h01.True", "SE.h00.True", 
                         "SE.h2.True", "SE.q0.True", "SE.q11.True"
                         
  )
  
  write.csv(Result1,
            sprintf("%sResult_PMMR_obs_N%0.4d_B%0.5d.csv",
                    FOLDER,
                    N, 
                    seed.data),
            row.names=F)
}

if(!Sim.Setup.Obs){
  IF.Numer.1 <-   (A)*(1-D)*q11hat.CF*(Y-h11hat.CF)+
    (1-A)*q0hat.CF*((1-D)*h11hat.CF-h01hat.CF)+h01hat.CF
  IF.Numer.0 <- (1-A)*q0hat.CF*(1-D)*Y
  IF.Denom   <- (1-A)*q0hat.CF*(1-D)
  
  IF.Numer.1.True <-   
    (A)*(1-D)*q11.true(Z,X1,X2)*(Y-h11.true(W,X1,X2))+
    (1-A)*q0.true(Z,X1,X2)*((1-D)*h11.true(W,X1,X2)-h01.true(W,X1,X2))+h01.true(W,X1,X2)
  IF.Numer.0.True <- 
    (1-A)*q0.true(Z,X1,X2)*((1-D)*Y)
  IF.Denom.True   <- 
    (1-A)*q0.true(Z,X1,X2)*(1-D)
  
  IF <- cbind(IF.Numer.1,
              IF.Numer.0,
              IF.Numer.1-IF.Numer.0,
              IF.Denom)
  Sigma <- as.matrix( var(IF)/N )
  
  Est.Numer.1 <- as.numeric( apply(IF,2,mean)[1] )
  Est.Numer.0 <- as.numeric( apply(IF,2,mean)[2] )
  Est.Numer   <- as.numeric( apply(IF,2,mean)[3] )
  Est.Denom   <- as.numeric( apply(IF,2,mean)[4] )
  Est.1       <- Est.Numer.1/Est.Denom
  Est.0       <- Est.Numer.0/Est.Denom
  Est         <- Est.Numer/Est.Denom
  
  IF.True <- cbind(IF.Numer.1.True,
                   IF.Numer.0.True,
                   IF.Numer.1.True-IF.Numer.0.True,
                   IF.Denom.True)
  Sigma.True <- as.matrix( var(IF.True)/N )
  
  Est.Numer.1.True <- as.numeric( apply(IF.True,2,mean)[1] )
  Est.Numer.0.True <- as.numeric( apply(IF.True,2,mean)[2] )
  Est.Numer.True   <- as.numeric( apply(IF.True,2,mean)[3] )
  Est.Denom.True   <- as.numeric( apply(IF.True,2,mean)[4] )
  Est.1.True       <- Est.Numer.1.True/Est.Denom.True
  Est.0.True       <- Est.Numer.0.True/Est.Denom.True
  Est.True         <- Est.Numer.True/Est.Denom.True
  
  contrast   <- matrix(c(0,0,1/Est.Denom,-Est/Est.Denom),1,4)
  SE         <- sqrt( contrast%*%Sigma%*%t(contrast) )
  
  contrast.True   <- matrix(c(0,0,1/Est.Denom.True,-Est.True/Est.Denom.True),1,4)
  SE.True    <- sqrt( contrast.True%*%Sigma.True%*%t(contrast.True) )
  
  IF.Est      <- (IF[,3]-Est*IF[,4])/Est.Denom
  IF.Est.True <- (IF.True[,3]-Est.True*IF.True[,4])/Est.Denom.True
  
  BOOT <- BOOT.True <- BOOT.GMM <- rep(0,10000)
  for(bb in 1:10000){
    vv <- mean( rnorm(N)*IF.Est )
    BOOT[bb] <- vv
    vv <- mean( rnorm(N)*IF.Est.True )
    BOOT.True[bb] <- vv
  }
  
  Result1 <- c(Est, SE, sd(BOOT), 
               Est.Numer.1, Est.Numer.0, Est.Denom,
               
               Est.True, SE.True, sd(BOOT.True), 
               Est.Numer.1.True, Est.Numer.0.True, Est.Denom.True,
               
               c( sd(h11hat.CF),  sd(q11hat.CF) ) , 
               c( sd(h11.true(W,X1,X2)),  sd(q11.true(Z,X1,X2)) )
               
  )
  
  Result1 <- matrix(Result1,1,length(Result1))
  
  colnames(Result1) <- c("Est",    "SE",    "SE.B",
                         "Est.Numer.1", "Est.Numer.0", "Est.Denom",
                         
                         "Est.True",    "SE.True",    "SE.B.True",
                         "Est.Numer.1.True", "Est.Numer.0.True", "Est.Denom.True",
                         
                         "SE.h11.PMMR", "SE.q11.PMMR",
                         "SE.h11.True", "SE.q11.True")
  
  write.csv(Result1,
            sprintf("%sResult_PMMR_exp_N%0.4d_B%0.5d.csv",
                    FOLDER,
                    N, 
                    seed.data),
            row.names=F)  
}

