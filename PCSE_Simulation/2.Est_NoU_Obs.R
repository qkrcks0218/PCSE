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

Sim.Setup.Obs <- T
# TRUE = Observational setting, FALSE = Experimental setting
seed.data <- BATCH
seed.CF   <- 1
N <- 500
FOLDER <- "WithNoU/"
set.seed( seed.data )

################################################################################
# Package and Source Files
################################################################################

source("0.Functions_PMMR.R")
source("0.MySL.R") 
source("0.DGP.R") 

################################################################################
# Parameters
################################################################################

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
# Data Generation : Observational
################################################################################

U    <- rnorm(N) 
X1 <- data.frame(matrix(rnorm(dimX*N),N,dimX))
X2 <- data.frame(matrix(rnorm(dimX*N),N,dimX))

PrA1 <- expit( betaa0 + 
                 apply(betaax1*X1,1,sum) + 
                 apply(betaax2*X2,1,sum) + 
                 betaau*U )
PrA0 <- 1-PrA1

A <- rbinom(N,1,PrA1)

Z  <- ( rnorm(N) + 
          betaz0 + 
          apply(betazx1*X1,1,sum) + 
          apply(betazx2*X2,1,sum) + 
          betazu*U + 
          betaza*A )
W  <- ( rnorm(N) + 
          betaw0 + 
          apply(betawx1*X1,1,sum) + 
          apply(betawx2*X2,1,sum) + 
          betawu*U )

PrD1A1 <- 1-apply(cbind(0,exp(betad0.1+
                                apply(betadx1.1*X1,1,sum) + 
                                apply(betadx2.1*X2,1,sum) + 
                                betadu.1*U),0.99),1,median)
PrD1A0 <- 1-apply(cbind(0,exp(betad0.0+
                                apply(betadx1.0*X1,1,sum) + 
                                apply(betadx2.0*X2,1,sum) + 
                                betadu.0*U),0.99),1,median)

D.Ad1  <- rbinom(N,1,PrD1A1)      
D.Ad0  <- rbinom(N,1,PrD1A0)      

Y.Ay1.Ad1 <- rnorm( N, betay0.1 + 
                      apply(betayx1.1*X1,1,sum) + 
                      apply(betayx2.1*X2,1,sum) + 
                      betayw.1*W + 
                      betayu.1*U, 1  )*(1-D.Ad1)
Y.Ay1.Ad0 <- rnorm( N, betay0.1 + 
                      apply(betayx1.1*X1,1,sum) + 
                      apply(betayx2.1*X2,1,sum) + 
                      betayw.1*W + 
                      betayu.1*U, 1  )*(1-D.Ad0)
Y.Ay0.Ad1 <- rnorm( N, betay0.0 + 
                      apply(betayx1.0*X1,1,sum) + 
                      apply(betayx2.0*X2,1,sum) + 
                      betayw.0*W + 
                      betayu.0*U, 1  )*(1-D.Ad1)
Y.Ay0.Ad0 <- rnorm( N, betay0.0 + 
                      apply(betayx1.0*X1,1,sum) + 
                      apply(betayx2.0*X2,1,sum) + 
                      betayw.0*W + 
                      betayu.0*U, 1  )*(1-D.Ad0)

D <- A*(D.Ad1)+(1-A)*(D.Ad0)
Y <- A*(Y.Ay1.Ad1)+(1-A)*(Y.Ay0.Ad0)

DataOriginal <- Data <- data.frame( cbind(Y,A,D,W,Z,X1,X2,U) )

#####################
# True nuis fts
#####################

h11.true <- function(w,x1,x2){
  gamma0  <- (betay0.1 - betayu.1*betaw0/betawu)
  gammax1 <- (betayx1.1-betayu.1*betawx1/betawu)
  gammax2 <- (betayx2.1-betayu.1*betawx2/betawu)
  gammaw  <- (betayw.1+betayu.1/betawu)
  
  return( 
    gamma0 + 
      apply(gammax1*x1,1,sum) + 
      apply(gammax2*x2,1,sum) + 
      gammaw*w
  )
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
  
  return( 
    (gamma0 + 
       apply(gammax1*x1,1,sum) + 
       apply(gammax2*x2,1,sum) + 
       gammaw*W)*
      exp(delta0 + 
            apply(deltax1*x1,1,sum) + 
            apply(deltax2*x2,1,sum) + 
            deltaw*W)
  )
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
  return( 
    (gamma0 + 
       apply(gammax1*x1,1,sum) + 
       apply(gammax2*x2,1,sum) + 
       gammaw*W)*
      exp(delta0 + 
            apply(deltax1*x1,1,sum) + 
            apply(deltax2*x2,1,sum) + 
            deltaw*W)
  )
}

h2.true <- function(w,x1,x2){
  deltaw  <- betadu.0/betawu
  deltax1 <- betadx1.0 - deltaw*betawx1
  deltax2 <- betadx2.0 - deltaw*betawx2
  delta0  <- betad0.0 - deltaw*betaw0 - 0.5*deltaw^2
  
  return( 
    exp(delta0 + 
          apply(deltax1*x1,1,sum) + 
          apply(deltax2*x2,1,sum) + 
          deltaw*W)
  )
}

q0.true <- function(z,x1,x2){
  gammaz  <- betaau/betazu
  gammax1 <- betaax1 -gammaz*betazx1
  gammax2 <- betaax2 -gammaz*betazx2
  gamma0  <- betaa0 - 0.5*gammaz^2 - gammaz*betaz0
  
  return(
    1 + exp( gamma0 + 
               apply(gammax1*x1,1,sum) + 
               apply(gammax2*x2,1,sum) + 
               gammaz*z )
    
  )
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
  
  return(
    exp(gamma0 + 
          apply(gammax1*x1,1,sum) + 
          apply(gammax2*x2,1,sum) + 
          gammaz*z) +
      exp(sgamma0 + 
            apply(sgammax1*x1,1,sum) + 
            apply(sgammax2*x2,1,sum) + 
            sgammaz*z)
  )
  
}

h11.true.vec <- h11.true(W,X1,X2)
h01.true.vec <- h01.true(W,X1,X2)
h00.true.vec <- h00.true(W,X1,X2)
h2.true.vec  <- h2.true(W,X1,X2)
q11.true.vec <- q11.true(Z,X1,X2)
q0.true.vec  <- q0.true(Z,X1,X2)

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
DataOriginal <- DataOriginal[Suffle,]

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

MEANMAT <- matrix(apply(Data[,c(Y.pos,Z.pos,W.pos,X.pos)],2,mean),N,length(c(Y.pos,Z.pos,W.pos,X.pos)),byrow=T)
SDMAT <- matrix(apply(Data[,c(Y.pos,Z.pos,W.pos,X.pos)],2,sd),N,length(c(Y.pos,Z.pos,W.pos,X.pos)),byrow=T)

Data[,c(Y.pos,Z.pos,W.pos,X.pos)] <- 
  (Data[,c(Y.pos,Z.pos,W.pos,X.pos)]-MEANMAT)/SDMAT

Y <- Data[,Y.pos]
A <- Data[,A.pos]
D <- Data[,D.pos]
W <- Data[,W.pos]
Z <- Data[,Z.pos]
X <- Data[,X.pos]
X1 <- X[,1:dimX]
X2 <- X[,dimX+1:dimX]
U <- Data[,U.pos]
WX <- cbind(W,X)
ZX <- cbind(Z,X)
rY <- range(Y)

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
  X1.MM [[ss]]  <- X[SS.Index[[ss]],1:dimX]
  X2.MM [[ss]]  <- X[SS.Index[[ss]],dimX+1:dimX]
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
# Cross Validation and Estimation
################################################################################

###################################################
# No U: Pr(A|X)
###################################################

ModelAgL <- ModelDgAL <- ModelYgADL <- list()
for(ss in 1:CF){
  ModelAgL[[ss]] <- 
    MySL(Data[ SS.Index[[ss]] , ],
         locY=A.pos,
         locX=c(W.pos,Z.pos,X.pos),
         Ydist=binomial(),
         SL.list=SL.hpara$SLL,
         MTRY=SL.hpara$MTRY,
         MLPL=SL.hpara$MLPL,
         NMN=SL.hpara$NMN,
         MLPdecay=SL.hpara$MLPdecay)
  
  ModelDgAL[[ss]] <- 
    MySL(Data[ SS.Index[[ss]] , ],
         locY=D.pos,
         locX=c(A.pos,W.pos,Z.pos,X.pos),
         Ydist=binomial(),
         SL.list=SL.hpara$SLL,
         MTRY=SL.hpara$MTRY,
         MLPL=SL.hpara$MLPL,
         NMN=SL.hpara$NMN,
         MLPdecay=SL.hpara$MLPdecay)
  
  ModelYgADL[[ss]] <-
    MySL(Data[ intersect(SS.Index[[ss]],which(D==0)) , ],
         locY=Y.pos,
         locX=c(A.pos,W.pos,Z.pos,X.pos),
         Ydist=gaussian(),
         SL.list=SL.hpara$SLL,
         MTRY=SL.hpara$MTRY,
         MLPL=SL.hpara$MLPL,
         NMN=SL.hpara$NMN,
         MLPdecay=SL.hpara$MLPdecay)
}

##

A1gL.NoCF <- rep(0,N)
for(ss in 1:CF){
  A1gL.NoCF[ SS.Index[[ss]] ] <- 
    predict(ModelAgL[[ss]], 
            newdata=Data[ SS.Index[[ss]] , c(W.pos,Z.pos,X.pos)])$pred
}

A1gL.CF <- rep(0,N)
for(ss in 1:CF){
  A1gL.CF[ SS.Index[[3-ss]] ] <- 
    predict(ModelAgL[[ss]], 
            newdata=Data[ SS.Index[[ss]] , c(W.pos,Z.pos,X.pos)])$pred
}

A1gL.NoCF.MM <- list()
for(ss in 1:CF){
  A1gL.NoCF.MM[[ss]]  <- A1gL.NoCF[ SS.Index[[ss]] ]
}

##

D0gA0L.NoCF <- D0gA1L.NoCF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[ss]] , c(A.pos,W.pos,Z.pos,X.pos)]
  
  Data.Temp[,1] <- 1
  D0gA1L.NoCF[ SS.Index[[ss]] ] <- 
    1 - predict(ModelDgAL[[ss]], newdata=Data.Temp)$pred
  
  Data.Temp[,1] <- 0
  D0gA0L.NoCF[ SS.Index[[ss]] ] <- 
    1 - predict(ModelDgAL[[ss]], newdata=Data.Temp)$pred
}

D0gA0L.CF <- D0gA1L.CF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[ss]] , c(A.pos,W.pos,Z.pos,X.pos)]
  
  Data.Temp[,1] <- 1
  D0gA1L.CF[ SS.Index[[3-ss]] ] <- 
    1 - predict(ModelDgAL[[ss]], newdata=Data.Temp)$pred
  
  Data.Temp[,1] <- 0
  D0gA0L.CF[ SS.Index[[3-ss]] ] <- 
    1 - predict(ModelDgAL[[ss]], newdata=Data.Temp)$pred
  
}

D0gA0L.NoCF.MM <- D0gA1L.NoCF.MM <- list()
for(ss in 1:CF){
  D0gA0L.NoCF.MM[[ss]]  <- D0gA0L.NoCF[ SS.Index[[ss]] ]
  D0gA1L.NoCF.MM[[ss]]  <- D0gA1L.NoCF[ SS.Index[[ss]] ]
}

##

YgA1D0L.NoCF <- YgA0D0L.NoCF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[ss]] , c(A.pos,W.pos,Z.pos,X.pos)]
  
  Data.Temp[,1] <- 1
  YgA1D0L.NoCF[ SS.Index[[ss]] ] <- 
    predict(ModelYgADL[[ss]], newdata=Data.Temp)$pred
  
  Data.Temp[,1] <- 0
  YgA0D0L.NoCF[ SS.Index[[ss]] ] <- 
    predict(ModelYgADL[[ss]], newdata=Data.Temp)$pred
}

YgA1D0L.CF <- YgA0D0L.CF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[ss]] , c(A.pos,W.pos,Z.pos,X.pos)]
  
  Data.Temp[,1] <- 1
  YgA1D0L.CF[ SS.Index[[3-ss]] ] <- 
    predict(ModelYgADL[[ss]], newdata=Data.Temp)$pred
  
  Data.Temp[,1] <- 0
  YgA0D0L.CF[ SS.Index[[3-ss]] ] <- 
    predict(ModelYgADL[[ss]], newdata=Data.Temp)$pred
  
}

YgA1D0L.NoCF.MM <- YgA0D0L.NoCF.MM <- list()
for(ss in 1:CF){
  YgA1D0L.NoCF.MM[[ss]]  <- YgA1D0L.NoCF[ SS.Index[[ss]] ]
  YgA0D0L.NoCF.MM[[ss]]  <- YgA0D0L.NoCF[ SS.Index[[ss]] ]
}

##

q11.NoU <- (1/A1gL.CF)*(D0gA0L.CF/D0gA1L.CF)
q0.NoU  <- (1/(1-A1gL.CF))
h11.NoU <- YgA1D0L.CF
h01.NoU <- YgA1D0L.CF*D0gA0L.CF
h00.NoU <- YgA0D0L.CF*D0gA0L.CF
h2.NoU  <- D0gA0L.CF

################################################################################
# Estimand
################################################################################

## Estimated

# IF.Numer.1 <- ((A)*(1-D)*q11hat.CF*(Y-h11hat.CF)+(1-A)*q0hat.CF*((1-D)*h11hat.CF-h01hat.CF)+h01hat.CF)*SDMAT[1,1]
# IF.Numer.0 <- ((1-A)*q0hat.CF*((1-D)*Y-h00hat.CF)+h00hat.CF)*SDMAT[1,1]
# IF.Denom   <- (1-A)*q0hat.CF*((1-D)-h2hat.CF)+h2hat.CF
# 
# IF <- cbind(IF.Numer.1,
#             IF.Numer.0,
#             IF.Numer.1-IF.Numer.0,
#             IF.Denom)
# Sigma <- as.matrix( var(IF)/N )
# 
# Est.Numer.1 <- as.numeric( apply(IF,2,mean)[1] )
# Est.Numer.0 <- as.numeric( apply(IF,2,mean)[2] )
# Est.Numer   <- as.numeric( apply(IF,2,mean)[3] )
# Est.Denom   <- as.numeric( apply(IF,2,mean)[4] )
# Est.1       <- Est.Numer.1/Est.Denom
# Est.0       <- Est.Numer.0/Est.Denom
# Est         <- Est.Numer/Est.Denom
# 
# IF.Est      <- (IF[,3]-Est*IF[,4])/Est.Denom
# 
# contrast   <- matrix(c(0,0,1/Est.Denom,-Est/Est.Denom),1,4)
# SE         <- sqrt( contrast%*%Sigma%*%t(contrast) )

## No U

IF.Numer.1.NoU <-   
  ((A)*(1-D)*q11.NoU*(Y-h11.NoU)+
  (1-A)*q0.NoU*((1-D)*h11.NoU-h01.NoU)+h01.NoU)*SDMAT[1,1]
IF.Numer.0.NoU <- 
  ((1-A)*q0.NoU*((1-D)*Y-h00.NoU)+h00.NoU)*SDMAT[1,1]
IF.Denom.NoU   <- 
  (1-A)*q0.NoU*((1-D)-h2.NoU)+h2.NoU

IF.NoU <- cbind(IF.Numer.1.NoU,
                IF.Numer.0.NoU,
                IF.Numer.1.NoU-IF.Numer.0.NoU,
                IF.Denom.NoU)
Sigma.NoU <- as.matrix( var(IF.NoU)/N )

Est.Numer.1.NoU <- as.numeric( apply(IF.NoU,2,mean)[1] )
Est.Numer.0.NoU <- as.numeric( apply(IF.NoU,2,mean)[2] )
Est.Numer.NoU   <- as.numeric( apply(IF.NoU,2,mean)[3] )
Est.Denom.NoU   <- as.numeric( apply(IF.NoU,2,mean)[4] )
Est.1.NoU       <- Est.Numer.1.NoU/Est.Denom.NoU
Est.0.NoU       <- Est.Numer.0.NoU/Est.Denom.NoU
Est.NoU         <- Est.Numer.NoU/Est.Denom.NoU

IF.Est.NoU <- (IF.NoU[,3]-Est.NoU*IF.NoU[,4])/Est.Denom.NoU

contrast.NoU   <- matrix(c(0,0,1/Est.Denom.NoU,-Est.NoU/Est.Denom.NoU),1,4)
SE.NoU    <- sqrt( contrast.NoU%*%Sigma.NoU%*%t(contrast.NoU) )

## True

IF.Numer.1.True <-   
  ((A)*(1-D)*q11.true.vec*(DataOriginal$Y-h11.true.vec)+
     (1-A)*q0.true.vec*((1-D)*h11.true.vec-h01.true.vec)+h01.true.vec)
IF.Numer.0.True <- 
  ((1-A)*q0.true.vec*((1-D)*DataOriginal$Y-h00.true.vec)+h00.true.vec)
IF.Denom.True   <- 
  (1-A)*q0.true.vec*((1-D)-h2.true.vec)+h2.true.vec

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


IF.Est.True <- (IF.True[,3]-Est.True*IF.True[,4])/Est.Denom.True

contrast.True   <- matrix(c(0,0,1/Est.Denom.True,-Est.True/Est.Denom.True),1,4)
SE.True    <- sqrt( contrast.True%*%Sigma.True%*%t(contrast.True) )


BOOT <- BOOT.True <- BOOT.NoU <- rep(0,10000)
for(bb in 1:10000){
  # vv <- mean( rnorm(N)*IF.Est )
  # BOOT[bb] <- vv
  vv <- mean( rnorm(N)*IF.Est.NoU )
  BOOT.NoU[bb] <- vv
  vv <- mean( rnorm(N)*IF.Est.True )
  BOOT.True[bb] <- vv
}

Result1 <- c(
  seed.data,
  # Est, SE, sd(BOOT), 
  # Est.Numer.1, Est.Numer.0, Est.Denom,
  
  Est.NoU, SE.NoU, sd(BOOT.NoU), 
  Est.Numer.1.NoU, Est.Numer.0.NoU, Est.Denom.NoU,
  
  Est.True, SE.True, sd(BOOT.True), 
  Est.Numer.1.True, Est.Numer.0.True, Est.Denom.True
  )

Result1 <- matrix(Result1,1,length(Result1))

colnames(Result1) <- 
  c("Seed",
    # "Est",    "SE",    "SE.B",
    # "Est.Numer.1", "Est.Numer.0", "Est.Denom",
    
    "Est.No",    "SE.No",    "SE.B.No",
    "Est.Numer.1.No", "Est.Numer.0.No", "Est.Denom.No",
    
    "Est.True",    "SE.True",    "SE.B.True",
    "Est.Numer.1.True", "Est.Numer.0.True", "Est.Denom.True"
    )

write.csv(Result1,
          sprintf("%sResult_NoU_obs_N%s_B%0.5d.csv",
                  FOLDER,
                  N, 
                  seed.data),
          row.names=F)


