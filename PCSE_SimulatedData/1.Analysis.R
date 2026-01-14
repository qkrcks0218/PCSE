################################################################################
# Cluster computing
################################################################################

args <- (commandArgs(trailingOnly=TRUE))
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

################################################################################
# Generate a Simulated Dataset
################################################################################

set.seed(1) # Fix a random seed for generating a simulated dataset
N <- 487
A <- rep(c(0,1),c(258,229))
PrA0 <- mean(1-A)
X01 <- rbinom(N,1,0.5)
X02 <- rbinom(N,1,0.5)
X03 <- rbinom(N,1,0.5)
X04 <- rbinom(N,1,0.5)
X05 <- rbinom(N,1,0.5)
U   <- rnorm(N)
X1  <- rbinom(N,1,0.2+0.05*A*(X01+X02+X03+X04+X05)+0.05*U)
W <- sqrt(1-0.2^2)*rnorm(N) + 0.2*U
Z <- sqrt(1-0.2^2)*rnorm(N) + 0.2*U
D <- rbinom(N,1,0.15+0.01*(X01+X02+X03+X04+X05)+0.1*X1+0.02*U+0.1*A)
Y <- rep(N)
Y <- -4.5-2*A - 0.5*(X01+X02+X03+X04+X05)-1*X1-2*U-rnorm(N)*20
Y <- Y*(1-D)

Data <- data.frame(cbind(Y,A,D,Z,W,X01,X02,X03,X04,X05,X1))

################################################################################
# Random Seed
################################################################################

seed.data <- BATCH
set.seed( seed.data ) 
N <- nrow(Data)

################################################################################
# Package and Source Files
################################################################################

source("0.Functions_PMMR.R")
source("0.MySL.R")  

###############################################################################
# Parameters
################################################################################

## PMMR regularization para
PL <- -8;      PU <- -2
BW.L <- 0;     BW.U <- 4

Para.Grid <- expand.grid(0, 
                         seq(BW.L,BW.U,by=0.5),
                         seq(PL,PU,by=1))

## Cross-fitting/Cross-validation
CF       <- 2
NumCV    <- 5
NumCVRep <- 5

## Superlearner parameters
SL.hpara <- list()
SL.hpara$SLL <- c(1,2,3,4,5,6,7,9,10)

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

N <- dim(Data)[1]

X0.pos  <- which( substr(colnames(Data),1,2)=="X0" )
X1.pos  <- which( substr(colnames(Data),1,2)=="X1" )
X.pos   <- which( substr(colnames(Data),1,1)=="X" ) 
W.pos   <- which( substr(colnames(Data),1,1)=="W" )
Z.pos   <- which( substr(colnames(Data),1,1)=="Z" )
A.pos   <- which( substr(colnames(Data),1,1)=="A" )
U.pos   <- which( substr(colnames(Data),1,1)=="U" )
Y.pos   <- which( substr(colnames(Data),1,1)=="Y" )
D.pos   <- which( substr(colnames(Data),1,1)=="D" )

Y    <- Data[,Y.pos]
A    <- Data[,A.pos]
D    <- Data[,D.pos]
W    <- Data[,W.pos]
Z    <- Data[,Z.pos]
X    <- Data[,X.pos]
X0   <- Data[,X0.pos]
X1   <- Data[,X1.pos] 
WX01 <- cbind(W,X0,X1)
ZX01 <- cbind(Z,X0,X1)
WX0  <- cbind(W,X0)
ZX0  <- cbind(Z,X0)

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
q.A0D0 <- round(quantile(1:(length(Random.A0D0)+1),seq(0,1,length=CF*NumCV+1)))
q.A0D1 <- round(quantile(1:(length(Random.A0D1)+1),seq(0,1,length=CF*NumCV+1)))
q.A1D0 <- round(quantile(1:(length(Random.A1D0)+1),seq(0,1,length=CF*NumCV+1)))
q.A1D1 <- round(quantile(1:(length(Random.A1D1)+1),seq(0,1,length=CF*NumCV+1)))
Suffle <- NULL
for(tt in 1:(2*NumCV)){
  Suffle <- c(Suffle,
              Random.A0D0[q.A0D0[tt]:(q.A0D0[tt+1]-1)],
              Random.A0D1[q.A0D1[tt]:(q.A0D1[tt+1]-1)],
              Random.A1D0[q.A1D0[tt]:(q.A1D0[tt+1]-1)],
              Random.A1D1[q.A1D1[tt]:(q.A1D1[tt+1]-1)])
}
Data <- Data[Suffle,]

SS.Index <- list()
SS.Cut   <- round((seq(0,N,length=CF+1)))
for(ss in 1:CF){
  SS.Index[[ss]] <- (SS.Cut[ss]+1):(SS.Cut[ss+1])
} 
CV.Index <- list()
CV.Index[[1]] <- list()
CV.Index[[2]] <- list()
for(ss in 1:CF){
  CV.CUT               <- round( seq(0,length(SS.Index[[ss]]),length=NumCV+1) )
  for(cv in 1:NumCV){
    CV.Index[[ss]][[cv]] <- ((CV.CUT)[cv]+1):(CV.CUT[cv+1])
  }
}

if(NumCVRep>1){
  
  REORDER <- NULL
  for(ss in 1:CF){
    CV.CUT         <- round( seq(0,length(SS.Index[[ss]]),length=NumCV+1) )
    for(tt in 2:NumCVRep){
      SS.Suffle.Index <- sample(1:length(SS.Index[[ss]]), length(SS.Index[[ss]]))
      for(cv in 1:NumCV){
        CV.Index[[ss]][[cv + (tt-1)*NumCV]] <- sort( SS.Suffle.Index[((CV.CUT)[cv]+1):(CV.CUT[cv+1])] )
      }  
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
# PMMR
################################################################################

###################################################
# Hyper-Parameters
###################################################

MEDDIST.WX0  <- log(median(dist(cbind(X0,W)))^2)
MEDDIST.ZX0  <- log(median(dist(cbind(X0,Z)))^2)
MEDDIST.WX01 <- log(median(dist(cbind(X0,X1,W)))^2)
MEDDIST.ZX01 <- log(median(dist(cbind(X0,X1,Z)))^2)

Para.Grid.WX0.Target <- Para.Grid.ZX0.Target <- Para.Grid
Para.Grid.WX0.Target[,2] <- Para.Grid.WX0.Target[,2] + MEDDIST.WX0
Para.Grid.ZX0.Target[,1] <- Para.Grid.ZX0.Target[,1] + MEDDIST.WX0
Para.Grid.WX0.Target[,1] <- Para.Grid.WX0.Target[,1] + MEDDIST.ZX0
Para.Grid.ZX0.Target[,2] <- Para.Grid.ZX0.Target[,2] + MEDDIST.ZX0

Para.Grid.WX01.Target <- Para.Grid.ZX01.Target <- Para.Grid
Para.Grid.WX01.Target[,2] <- Para.Grid.WX01.Target[,2] + MEDDIST.WX01
Para.Grid.ZX01.Target[,1] <- Para.Grid.ZX01.Target[,1] + MEDDIST.WX01
Para.Grid.WX01.Target[,1] <- Para.Grid.WX01.Target[,1] + MEDDIST.ZX01
Para.Grid.ZX01.Target[,2] <- Para.Grid.ZX01.Target[,2] + MEDDIST.ZX01

###################################################
# Bridge Function: h11
###################################################

Opt.Para.h11 <- list()
for(ss in 1:CF){
  
  CV.Curve <- apply(Para.Grid.WX01.Target,
                    1,
                    FUN=function(vv){ CrossVar(vv,
                                               ss,
                                               subset=posA1D0,
                                               response=Y*(1-D)*A,
                                               target=WX01,
                                               perturb=ZX01,
                                               diagonal=(1-D)*A) })
  
  Opt.Para.h11[[ss]] <-
    as.numeric(Para.Grid.WX01.Target[which.min(CV.Curve),])
}

pos.h11.MM <- list()
h11.MM <- list()
h11.predict <- list()

for(ss in 1:CF){
  
  pos.h11.MM[[ss]] <- SS.posA1D0[[ss]]
  
  h11.MM[[ss]] <- 
    FT_PMMR( Y         =(Y.MM[[ss]]*(A.MM[[ss]])*(1-D.MM[[ss]]))[pos.h11.MM[[ss]]],
             Perturb   =ZX01.MM[[ss]][pos.h11.MM[[ss]],],
             Target    =WX01.MM[[ss]][pos.h11.MM[[ss]],],
             Diagonal  =((A.MM[[ss]])*(1-D.MM[[ss]]))[pos.h11.MM[[ss]]],
             Perturb.bw=exp(Opt.Para.h11[[ss]][1]),
             Target.bw =exp(Opt.Para.h11[[ss]][2]),
             lambda    =exp(Opt.Para.h11[[ss]][3]),
             NV        =FALSE)
}

h11.predict[[1]] <- function(WX.New.Input){
  FT_RBF(X     = WX01.MM[[1]][pos.h11.MM[[1]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h11[[1]][2]))%*%h11.MM[[1]]$alpha + h11.MM[[1]]$intercept
}

h11.predict[[2]] <- function(WX.New.Input){
  FT_RBF(X     = WX01.MM[[2]][pos.h11.MM[[2]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h11[[2]][2]))%*%h11.MM[[2]]$alpha + h11.MM[[2]]$intercept
}

h11hat.NoCF <- h11hat.CF <- rep(0,N)
for(ss in 1:CF){
  h11hat.NoCF[SS.Index[[ss]]] <- h11.predict[[ss]](WX01.MM[[ss]])
  h11hat.CF[SS.Index[[3-ss]]] <- h11.predict[[ss]](WX01.MM[[3-ss]])
}

h11hat.NoCF.MM <- list()
for(ss in 1:CF){
  print(ss)
  h11hat.NoCF.MM[[ss]]  <- h11hat.NoCF[SS.Index[[ss]] ]
}

###################################################
# Bridge Function: q0
###################################################

q0hat.NoCF <- q0hat.CF <- rep(1/PrA0,N)

q0hat.NoCF.MM <- list()
for(ss in 1:CF){
  print(ss)
  q0hat.NoCF.MM[[ss]] <- q0hat.NoCF[SS.Index[[ss]] ]
}

###################################################
# Bridge Function: q11
###################################################

Opt.Para.q11 <- list()
for(ss in 1:CF){
  
  CV.Curve <- apply(Para.Grid.ZX01.Target,
                    1,
                    FUN=function(vv){ CrossVar(vv,
                                               ss,
                                               subset = posD0,
                                               response = (1-A)*q0hat.NoCF*(1-D),
                                               target=ZX01,
                                               perturb=WX01,
                                               diagonal = A*(1-D)) })
  
  Opt.Para.q11[[ss]] <-
    as.numeric(Para.Grid.ZX01.Target[which.min(CV.Curve),])
}

pos.q11.MM <- list()
q11.MM <- list()
q11.predict <- list()

for(ss in 1:CF){
  
  pos.q11.MM[[ss]] <- SS.posD0[[ss]]
  
  Resp.Temp <- (1-A.MM[[ss]])*(1-D.MM[[ss]])*(q0hat.NoCF.MM[[ss]])
  Diag.Temp <- A.MM[[ss]]*(1-D.MM[[ss]])
  
  q11.MM[[ss]] <- FT_PMMR( Y         =Resp.Temp[pos.q11.MM[[ss]]],
                           Perturb   =WX01.MM[[ss]][pos.q11.MM[[ss]],],
                           Target    =ZX01.MM[[ss]][pos.q11.MM[[ss]],],
                           Diagonal  =Diag.Temp[pos.q11.MM[[ss]]],
                           Perturb.bw=exp(Opt.Para.q11[[ss]][1]),
                           Target.bw =exp(Opt.Para.q11[[ss]][2]),
                           lambda    =exp(Opt.Para.q11[[ss]][3]),
                           NV        =FALSE)
}

q11.predict[[1]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX01.MM[[1]][pos.q11.MM[[1]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q11[[1]][2]))%*%q11.MM[[1]]$alpha + q11.MM[[1]]$intercept
}

q11.predict[[2]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX01.MM[[2]][pos.q11.MM[[2]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q11[[2]][2]))%*%q11.MM[[2]]$alpha + q11.MM[[2]]$intercept
}

q11hat.NoCF <- q11hat.CF <- rep(0,N)
for(ss in 1:CF){
  q11hat.NoCF[SS.Index[[ss]]] <- q11.predict[[ss]](ZX01.MM[[ss]])
  q11hat.CF[SS.Index[[3-ss]]] <- q11.predict[[ss]](ZX01.MM[[3-ss]])
}

q11hat.NoCF.MM <- list()
for(ss in 1:CF){
  print(ss)
  q11hat.NoCF.MM[[ss]]  <- q11hat.NoCF[SS.Index[[ss]] ]
}

################################################################################
# Estimand
################################################################################

## Estimated

IF.Numer.1 <- ((A)*(1-D)*q11hat.CF*(Y-h11hat.CF)+
                 (1-A)*q0hat.CF*(1-D)*h11hat.CF)
IF.Numer.0 <- ((1-A)*q0hat.CF*(1-D)*Y)
IF.Denom   <- (1-A)*q0hat.CF*(1-D)

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

IF.Est      <- (IF[,3]-Est*IF[,4])/Est.Denom

contrast   <- matrix(c(0,0,1/Est.Denom,-Est/Est.Denom,
                       1/Est.Denom,0,0,-Est.1/Est.Denom,
                       0,1/Est.Denom,0,-Est.0/Est.Denom),3,4,byrow=T)
SE         <- sqrt( diag(contrast%*%Sigma%*%t(contrast)) )


Result1 <- c(seed.data,
             Est, SE[1],
             Est.1, SE[2],
             Est.0, SE[3])

Result1 <- matrix(Result1,1,length(Result1))

colnames(Result1) <- 
  c("Seed",
    "Est",   "SE",
    "Est.1", "SE.1",
    "Est.0", "SE.0") 

write.csv(Result1,
          sprintf("Result_PMMR_B%0.5d.csv",
                  seed.data), 
          row.names=F)

################################################################################
# Ignorability
################################################################################

# No variation in X will cause an error for ML approaches, so add very tiny noise
Data[,X.pos] <- Data[,X.pos] + 
  matrix(rnorm(nrow(Data)*5)*0.0001,nrow(Data),6)

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

###################################################
# No U: Pr(A|X)
###################################################

ModelAgD0L01 <- 
  ModelYgA1D0L01 <- list()
for(ss in 1:CF){
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

q11.NoU <- (1/(PrA0))*((1-A1gD0L01.CF)/A1gD0L01.CF)
q0.NoU  <- (1/PrA0)
h11.NoU <- YgA1D0L01.CF

################################################################################
# Estimand
################################################################################

## No U

IF.Numer.1.NoU <-   
  ((A)*(1-D)*q11.NoU*(Y-h11.NoU)+
     (1-A)*q0.NoU*((1-D)*h11.NoU)) 
IF.Numer.0.NoU <- 
  ((1-A)*q0.NoU*((1-D)*Y)) 
IF.Denom.NoU   <- 
  (1-A)*q0.NoU*((1-D))

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

IF.Est.NoU      <- (IF.NoU[,3]-Est.NoU*IF.NoU[,4])/Est.Denom.NoU

contrast.NoU   <- matrix(c(0,0,1/Est.Denom.NoU,-Est.NoU/Est.Denom.NoU,
                           1/Est.Denom.NoU,0,0,-Est.1.NoU/Est.Denom.NoU,
                           0,1/Est.Denom.NoU,0,-Est.0.NoU/Est.Denom.NoU),3,4,byrow=T)
SE.NoU         <- sqrt( diag(contrast.NoU%*%Sigma.NoU%*%t(contrast.NoU)) )


Result2 <- c(seed.data,
             Est.NoU, SE.NoU[1],
             Est.1.NoU, SE.NoU[2],
             Est.0.NoU, SE.NoU[3])

Result2 <- matrix(Result2,1,length(Result2))

colnames(Result2) <- 
  c("Seed",
    "Est",   "SE",
    "Est.1", "SE.1",
    "Est.0", "SE.0")  

write.csv(Result2,
          sprintf("Result_Ign_B%0.5d.csv", 
                  seed.data),
          row.names=F)

