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

Each.Data <- 1000
Sim.Para <- expand.grid(1:Each.Data,c(1,0),(1:4)*500)
sd.U <- Sim.Para[BATCH,2]
N <- Sim.Para[BATCH,3]
seed.data <- Sim.Para[BATCH,1]

seed.CF   <- 1

FOLDER <- ""
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
SL.hpara$SLL <- c(1,2,3,4,
                  5,
                  6,7,9,10)

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
NumCVRep <- 5

################################################################################
# Data Generation : Observational
################################################################################

U  <- rnorm(N)*sd.U
X0 <- matrix(rnorm(N*dimX0),N,dimX0)

PrA1 <- expit( betaa0 + 
                 apply(betaax0*X0,1,sum) +
                 betaau*U )
PrA0 <- 1-PrA1

A <- rbinom(N,1,PrA1)

Z  <- ( rnorm(N) + 
          betaz0 + 
          apply(betazx0*X0,1,sum) + 
          betazu*U + 
          betaza*A )
W  <- ( rnorm(N) + 
          betaw0 + 
          apply(betawx0*X0,1,sum) +
          betawu*U )

X1A1 <- betax10.1+betax1x0.1*X0+betax1u.1*matrix(U,N,dimX0)+rnorm(N*dimX0)
X1A0 <- betax10.0+betax1x0.0*X0+betax1u.0*matrix(U,N,dimX0)+rnorm(N*dimX0)

PrD1A1 <- 1-apply(cbind(0,exp(betad0.1+
                                apply(betadx0.1*X0,1,sum) + 
                                apply(betadx1.1*X1A1,1,sum) + 
                                betadu.1*U),0.99),1,median)
PrD1A0 <- 1-apply(cbind(0,exp(betad0.0+
                                apply(betadx0.0*X0,1,sum) +
                                apply(betadx1.0*X1A0,1,sum) +
                                betadu.0*U),0.99),1,median)

D.Ad1  <- rbinom(N,1,PrD1A1)      
D.Ad0  <- rbinom(N,1,PrD1A0)      

Y.Ay1.Ad1 <- rnorm( N, betay0.1 + 
                      apply(betayx0.1*X0,1,sum) +
                      apply(betayx1.1*X1A1,1,sum) +
                      betayw.1*W + 
                      betayu.1*U, 1  )*(1-D.Ad1)
Y.Ay1.Ad0 <- rnorm( N, betay0.1 + 
                      apply(betayx0.1*X0,1,sum) +
                      apply(betayx1.1*X1A0,1,sum) +
                      betayw.1*W + 
                      betayu.1*U, 1  )*(1-D.Ad0)
Y.Ay0.Ad1 <- rnorm( N, betay0.0 + 
                      apply(betayx0.0*X0,1,sum) +
                      apply(betayx1.0*X1A1,1,sum) +
                      betayw.0*W + 
                      betayu.0*U, 1  )*(1-D.Ad1)
Y.Ay0.Ad0 <- rnorm( N, betay0.0 + 
                      apply(betayx0.0*X0,1,sum) +
                      apply(betayx1.0*X1A0,1,sum) +
                      betayw.0*W + 
                      betayu.0*U, 1  )*(1-D.Ad0)

X1 <- A*X1A1 + (1-A)*X1A0
D <- A*(D.Ad1)+(1-A)*(D.Ad0)
Y <- A*(Y.Ay1.Ad1)+(1-A)*(Y.Ay0.Ad0)

DataOriginal <- Data <- data.frame( cbind(Y,A,D,W,Z,X0,X1) )
colnames(DataOriginal) <- colnames(Data) <- 
  c("Y","A","D","W","Z",
    sprintf("X0.%d",1:dimX0),
    sprintf("X1.%d",1:dimX0))

#####################
# True nuis fts
#####################

h11.true.vec <- h11.true(W,X0,X1)
h01.true.vec <- h01.true(W,X0)
h00.true.vec <- h00.true(W,X0)
h2.true.vec  <- h2.true(W,X0)
q0.true.vec  <- q0.true(Z,X0)
q11.true.vec <- q11.true(Z,X0,X1)

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
q.A0D0 <- round(quantile(1:(length(Random.A0D0)+1),seq(0,1,length=1+2*NumCV)))
q.A0D1 <- round(quantile(1:(length(Random.A0D1)+1),seq(0,1,length=1+2*NumCV)))
q.A1D0 <- round(quantile(1:(length(Random.A1D0)+1),seq(0,1,length=1+2*NumCV)))
q.A1D1 <- round(quantile(1:(length(Random.A1D1)+1),seq(0,1,length=1+2*NumCV)))
Suffle <- NULL
for(tt in 1:(2*NumCV)){
  Suffle <- c(Suffle,
              Random.A0D0[q.A0D0[tt]:(q.A0D0[tt+1]-1)],
              Random.A0D1[q.A0D1[tt]:(q.A0D1[tt+1]-1)],
              Random.A1D0[q.A1D0[tt]:(q.A1D0[tt+1]-1)],
              Random.A1D1[q.A1D1[tt]:(q.A1D1[tt+1]-1)])
}
Data <- Data[Suffle,]
DataOriginal <- DataOriginal[Suffle,]

h11.true.vec <- h11.true.vec [Suffle]
h01.true.vec <- h01.true.vec [Suffle]
h00.true.vec <- h00.true.vec [Suffle]
h2.true.vec  <- h2.true.vec  [Suffle]
q11.true.vec <- q11.true.vec [Suffle]
q0.true.vec  <- q0.true.vec  [Suffle]

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
    SFF <- sapply(1:(N/2/NumCV),function(t){sample(NumCV,NumCV)})
    for(cv in 1:NumCV){
      CV.Index[[cv + (tt-1)*NumCV]] <- 
        sapply(1:(N/2/NumCV),function(t){REORDER[SFF[cv,t],t]})
    }  
  }
}

X0.pos <- which( substr(colnames(Data),1,2)=="X0" )
X1.pos <- which( substr(colnames(Data),1,2)=="X1" )
X.pos <- which( substr(colnames(Data),1,1)=="X" )
W.pos <- which( substr(colnames(Data),1,1)=="W" )
Z.pos <- which( substr(colnames(Data),1,1)=="Z" )
A.pos <- which( substr(colnames(Data),1,1)=="A" )
U.pos <- which( substr(colnames(Data),1,1)=="U" )
Y.pos <- which( substr(colnames(Data),1,1)=="Y" )
D.pos <- which( substr(colnames(Data),1,1)=="D" )

MEANMAT <- matrix(apply(Data[,c(Y.pos,Z.pos,W.pos,X.pos)],2,mean),
                  N,length(c(Y.pos,Z.pos,W.pos,X.pos)),byrow=T)
SDMAT <- matrix(apply(Data[,c(Y.pos,Z.pos,W.pos,X.pos)],2,sd),
                N,length(c(Y.pos,Z.pos,W.pos,X.pos)),byrow=T)

Data[,c(Y.pos,Z.pos,W.pos,X.pos)] <- 
  (Data[,c(Y.pos,Z.pos,W.pos,X.pos)]-MEANMAT)/SDMAT

Y <- Data[,Y.pos]
A <- Data[,A.pos]
D <- Data[,D.pos]
W <- Data[,W.pos]
Z <- Data[,Z.pos]
X <- Data[,X.pos]
X0 <- Data[,X0.pos]
X1 <- Data[,X1.pos]
U <- Data[,U.pos]
WX01 <- cbind(W,X0,X1)
ZX01 <- cbind(Z,X0,X1)
WX0 <- cbind(W,X0)
ZX0 <- cbind(Z,X0)
rY <- range(Y)

Y.MM  <- list()
A.MM  <- list()
D.MM  <- list()
WX01.MM <- list()
ZX01.MM <- list()
WX0.MM <- list()
ZX0.MM <- list()

for(ss in 1:CF){
  
  Y.MM [[ss]]  <- Y[SS.Index[[ss]] ]
  A.MM [[ss]]  <- A[SS.Index[[ss]] ]
  D.MM [[ss]]  <- D[SS.Index[[ss]] ]
  WX01.MM[[ss]]  <- WX01[SS.Index[[ss]], ]
  ZX01.MM[[ss]]  <- ZX01[SS.Index[[ss]], ]
  WX0.MM[[ss]]  <- WX0[SS.Index[[ss]], ]
  ZX0.MM[[ss]]  <- ZX0[SS.Index[[ss]], ]
  
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

ModelAgL0 <- ModelAgD0L01 <- 
  ModelYgA1D0L01 <- 
  ModelYgA0D0L0 <- ModelDgA0L0 <- list()
for(ss in 1:CF){
  ModelAgL0[[ss]] <- 
    MySL(Data[ SS.Index[[ss]] , ],
         locY=A.pos,
         locX=c(W.pos,Z.pos,X0.pos),
         Ydist=binomial(),
         SL.list=SL.hpara$SLL,
         MTRY=SL.hpara$MTRY,
         MLPL=SL.hpara$MLPL,
         NMN=SL.hpara$NMN,
         MLPdecay=SL.hpara$MLPdecay)
  
  ModelAgD0L01[[ss]] <- 
    MySL(Data[ intersect(SS.Index[[ss]],which(D==0)) , ],
         locY=A.pos,
         locX=c(W.pos,Z.pos,X0.pos,X1.pos),
         Ydist=binomial(),
         SL.list=SL.hpara$SLL,
         MTRY=SL.hpara$MTRY,
         MLPL=SL.hpara$MLPL,
         NMN=SL.hpara$NMN,
         MLPdecay=SL.hpara$MLPdecay)
  
  ModelYgA1D0L01[[ss]] <-
    MySL(Data[ intersect(SS.Index[[ss]],which(D==0 & A==1)) , ],
         locY=Y.pos,
         locX=c(W.pos,Z.pos,X0.pos,X1.pos),
         Ydist=gaussian(),
         SL.list=SL.hpara$SLL,
         MTRY=SL.hpara$MTRY,
         MLPL=SL.hpara$MLPL,
         NMN=SL.hpara$NMN,
         MLPdecay=SL.hpara$MLPdecay)
  
  ModelYgA0D0L0[[ss]] <-
    MySL(Data[ intersect(SS.Index[[ss]],which(D==0 & A==0)) , ],
         locY=Y.pos,
         locX=c(W.pos,Z.pos,X0.pos),
         Ydist=gaussian(),
         SL.list=SL.hpara$SLL,
         MTRY=SL.hpara$MTRY,
         MLPL=SL.hpara$MLPL,
         NMN=SL.hpara$NMN,
         MLPdecay=SL.hpara$MLPdecay)
  
  ModelDgA0L0[[ss]] <- 
    MySL(Data[ intersect(SS.Index[[ss]],which(A==0)) , ],
         locY=D.pos,
         locX=c(W.pos,Z.pos,X0.pos),
         Ydist=binomial(),
         SL.list=SL.hpara$SLL,
         MTRY=SL.hpara$MTRY,
         MLPL=SL.hpara$MLPL,
         NMN=SL.hpara$NMN,
         MLPdecay=SL.hpara$MLPdecay)
  
}

##

A1gL0.NoCF <- rep(0,N)
for(ss in 1:CF){
  A1gL0.NoCF[ SS.Index[[ss]] ] <- 
    predict(ModelAgL0[[ss]], 
            newdata=Data[ SS.Index[[ss]] , c(W.pos,Z.pos,X0.pos)])$pred
}

A1gL0.CF <- rep(0,N)
for(ss in 1:CF){
  A1gL0.CF[ SS.Index[[3-ss]] ] <- 
    predict(ModelAgL0[[ss]], 
            newdata=Data[ SS.Index[[3-ss]] , c(W.pos,Z.pos,X0.pos)])$pred
}

A1gL0.NoCF.MM <- list()
for(ss in 1:CF){
  A1gL0.NoCF.MM[[ss]]  <- A1gL0.NoCF[ SS.Index[[ss]] ]
}

##

A1gD0L01.NoCF <- rep(0,N)
for(ss in 1:CF){
  A1gD0L01.NoCF[ SS.Index[[ss]] ] <- 
    predict(ModelAgD0L01[[ss]], 
            newdata=Data[ SS.Index[[ss]] , c(W.pos,Z.pos,X0.pos,X1.pos)])$pred
}

A1gD0L01.CF <- rep(0,N)
for(ss in 1:CF){
  A1gD0L01.CF[ SS.Index[[3-ss]] ] <- 
    predict(ModelAgD0L01[[ss]], 
            newdata=Data[ SS.Index[[3-ss]] , c(W.pos,Z.pos,X0.pos,X1.pos)])$pred
}

A1gD0L01.NoCF.MM <- list()
for(ss in 1:CF){
  A1gD0L01.NoCF.MM[[ss]]  <- A1gD0L01.NoCF[ SS.Index[[ss]] ]
}

##

YgA1D0L01.NoCF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[ss]] , c(W.pos,Z.pos,X0.pos,X1.pos)]
  YgA1D0L01.NoCF[ SS.Index[[ss]] ] <- 
    predict(ModelYgA1D0L01[[ss]], newdata=Data.Temp)$pred 
}

YgA1D0L01.CF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[3-ss]] , c(W.pos,Z.pos,X0.pos,X1.pos)]
  YgA1D0L01.CF[ SS.Index[[3-ss]] ] <- 
    predict(ModelYgA1D0L01[[ss]], newdata=Data.Temp)$pred 
  
}

YgA1D0L01.NoCF.MM <- list()
for(ss in 1:CF){
  YgA1D0L01.NoCF.MM[[ss]]  <- YgA1D0L01.NoCF[ SS.Index[[ss]] ] 
}

##

YgA0D0L0.NoCF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[ss]] , c(W.pos,Z.pos,X0.pos)]
  YgA0D0L0.NoCF[ SS.Index[[ss]] ] <- 
    predict(ModelYgA0D0L0[[ss]], newdata=Data.Temp)$pred
}

YgA0D0L0.CF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[3-ss]] , c(W.pos,Z.pos,X0.pos)]
  YgA0D0L0.CF[ SS.Index[[3-ss]] ] <- 
    predict(ModelYgA0D0L0[[ss]], newdata=Data.Temp)$pred
  
}

YgA0D0L0.NoCF.MM <- list()
for(ss in 1:CF){
  YgA0D0L0.NoCF.MM[[ss]]  <- YgA0D0L0.NoCF[ SS.Index[[ss]] ]
}

##

D0gA0L0.NoCF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[ss]] , c(W.pos,Z.pos,X0.pos)]
  
  D0gA0L0.NoCF[ SS.Index[[ss]] ] <- 
    1 - predict(ModelDgA0L0[[ss]], newdata=Data.Temp)$pred
}

D0gA0L0.CF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[3-ss]] , c(W.pos,Z.pos,X0.pos)]
  D0gA0L0.CF[ SS.Index[[3-ss]] ] <- 
    1 - predict(ModelDgA0L0[[ss]], newdata=Data.Temp)$pred
  
}

D0gA0L0.NoCF.MM <- list()
for(ss in 1:CF){
  D0gA0L0.NoCF.MM[[ss]]  <- D0gA0L0.NoCF[ SS.Index[[ss]] ]
}

##

Data.Extend <- cbind(Data,(1-D)*YgA1D0L01.NoCF)
colnames(Data.Extend) <- c(colnames(Data),"V")
V.pos <- which(colnames(Data.Extend)=="V")

Modelh01 <- list()

for(ss in 1:CF){
  
  Modelh01[[ss]] <-
    MySL(Data.Extend[ intersect(SS.Index[[ss]],which(A==0)) , ],
         locY=V.pos,
         locX=c(W.pos,Z.pos,X0.pos),
         Ydist=gaussian(),
         SL.list=SL.hpara$SLL,
         MTRY=SL.hpara$MTRY,
         MLPL=SL.hpara$MLPL,
         NMN=SL.hpara$NMN,
         MLPdecay=SL.hpara$MLPdecay)
  
}

##

h01.NoCF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[ss]] , c(W.pos,Z.pos,X0.pos)]
  
  h01.NoCF[ SS.Index[[ss]] ] <- 
    predict(Modelh01[[ss]], newdata=Data.Temp)$pred
}

h01.CF <- rep(0,N)
for(ss in 1:CF){
  
  Data.Temp <- Data[ SS.Index[[3-ss]] , c(W.pos,Z.pos,X0.pos)]
  h01.CF[ SS.Index[[3-ss]] ] <- 
    predict(Modelh01[[ss]], newdata=Data.Temp)$pred
  
}

h01.NoCF.MM <- list()
for(ss in 1:CF){
  h01.NoCF.MM[[ss]]  <- h01.NoCF[ SS.Index[[ss]] ]
}







##

q11.NoU <- (1/(1-A1gL0.CF))*((1-A1gD0L01.CF)/A1gD0L01.CF)
q0.NoU  <- (1/(1-A1gL0.CF))
h11.NoU <- YgA1D0L01.CF
h01.NoU <- h01.CF
h00.NoU <- YgA0D0L0.CF*D0gA0L0.CF
h2.NoU  <- D0gA0L0.CF

################################################################################
# Estimand
################################################################################

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

BOOT <- BOOT.NoU <- rep(0,10000)
for(bb in 1:10000){
  vv <- mean( rnorm(N)*IF.Est.NoU )
  BOOT.NoU[bb] <- vv
}

Result1 <- c(Est.NoU, SE.NoU, sd(BOOT.NoU), 
             Est.Numer.1.NoU, Est.Numer.0.NoU, Est.Denom.NoU)

Result1 <- matrix(Result1,1,length(Result1))

colnames(Result1) <- 
  c("Est.No",    "SE.No",    "SE.B.No",
    "Est.Numer.1.No", "Est.Numer.0.No", "Est.Denom.No")

write.csv(Result1,
          sprintf("%sResult_Ign_obs_N%s_sdU%d_B%0.5d.csv",
                  FOLDER,
                  N, 
                  sd.U,
                  seed.data),
          row.names=F)


