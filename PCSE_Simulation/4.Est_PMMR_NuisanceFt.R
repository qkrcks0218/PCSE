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

Each.Data <- 520
Sim.Setup.Obs <- T

Sim.Para <- rbind(expand.grid(1:Each.Data,seq(8,16,by=1)))

nn <- Sim.Para[BATCH,2]
N <- round( (50*2^(nn/2))/2 )*2
N.0 <- 10^(4)

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

## PMMR regularization para
PL <- -8;      PU <- -2
BW.L <- 0;     BW.U <- 4

Para.Grid <- expand.grid(0, 
                         seq(BW.L,BW.U,by=2),
                         seq(PL,PU,by=3)) 

## Cross-fitting/Cross-validation
CF       <- 2
NumCV    <- 5
NumCVRep <- 3

################################################################################
# Data Generation : Observational
################################################################################

U  <- rnorm(N)
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
q.A0D0 <- round(quantile(1:(length(Random.A0D0)+1),seq(0,1,length=(2*NumCV)+1)))
q.A0D1 <- round(quantile(1:(length(Random.A0D1)+1),seq(0,1,length=(2*NumCV)+1)))
q.A1D0 <- round(quantile(1:(length(Random.A1D0)+1),seq(0,1,length=(2*NumCV)+1)))
q.A1D1 <- round(quantile(1:(length(Random.A1D1)+1),seq(0,1,length=(2*NumCV)+1)))
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
for(ss in 1:1){
  
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

for(ss in 1:1){
  
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
for(ss in 1:1){
  h11hat.NoCF[SS.Index[[ss]]] <- h11.predict[[ss]](WX01.MM[[ss]])
  h11hat.CF[SS.Index[[3-ss]]] <- h11.predict[[ss]](WX01.MM[[3-ss]])
}

h11hat.NoCF.MM <- list()
for(ss in 1:1){
  print(ss)
  h11hat.NoCF.MM[[ss]]  <- h11hat.NoCF[SS.Index[[ss]] ]
}

###################################################
# Bridge Function: h01
###################################################


Opt.Para.h01 <- list()
for(ss in 1:1){
  
  
  CV.Curve <- apply(Para.Grid.WX0.Target,
                    1,
                    FUN=function(vv){ CrossVar(vv,
                                               ss,
                                               subset=posA0,
                                               response=h11hat.NoCF*(1-D)*(1-A),
                                               target=WX0,
                                               perturb=ZX0,
                                               diagonal=(1-A)) })
  
  Opt.Para.h01[[ss]] <-
    as.numeric(Para.Grid.WX0.Target[which.min(CV.Curve),])
}


pos.h01.MM <- list()
h01.MM <- list()
h01.predict <- list()

for(ss in 1:1){
  
  pos.h01.MM[[ss]] <- SS.posA0[[ss]]
  Resp.Temp <- h11hat.NoCF.MM[[ss]]*(1-A.MM[[ss]])*(1-D.MM[[ss]])
  
  h01.MM[[ss]] <- 
    FT_PMMR( Y         =Resp.Temp[pos.h01.MM[[ss]]],
             Perturb   =ZX0.MM[[ss]][pos.h01.MM[[ss]],],
             Target    =WX0.MM[[ss]][pos.h01.MM[[ss]],],
             Diagonal  =(1-A.MM[[ss]])[pos.h01.MM[[ss]]],
             Perturb.bw=exp(Opt.Para.h01[[ss]][1]),
             Target.bw =exp(Opt.Para.h01[[ss]][2]),
             lambda    =exp(Opt.Para.h01[[ss]][3]),
             NV        =FALSE)
  
}


h01.predict[[1]] <- function(WX.New.Input){
  FT_RBF(X     = WX0.MM[[1]][pos.h01.MM[[1]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h01[[1]][2]))%*%h01.MM[[1]]$alpha + h01.MM[[1]]$intercept
}

h01.predict[[2]] <- function(WX.New.Input){
  FT_RBF(X     = WX0.MM[[2]][pos.h01.MM[[2]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h01[[2]][2]))%*%h01.MM[[2]]$alpha + h01.MM[[2]]$intercept
}

h01hat.NoCF <- h01hat.CF <- rep(0,N)
for(ss in 1:1){
  h01hat.NoCF[SS.Index[[ss]]] <- h01.predict[[ss]](WX0.MM[[ss]])
  h01hat.CF[SS.Index[[3-ss]]] <- h01.predict[[ss]](WX0.MM[[3-ss]])
}

h01hat.NoCF.MM <- list()
for(ss in 1:1){
  print(ss)
  h01hat.NoCF.MM[[ss]]  <- h01hat.NoCF[SS.Index[[ss]] ]
}

###################################################
# Bridge Function: h00
###################################################


Opt.Para.h00 <- list()
for(ss in 1:1){
  
  CV.Curve <- apply(Para.Grid.WX0.Target,
                    1,
                    FUN=function(vv){ CrossVar(vv,
                                               ss,
                                               subset=posA0,
                                               response=Y*(1-D)*(1-A),
                                               target=WX0,
                                               perturb=ZX0,
                                               diagonal=(1-A)) })
  
  Opt.Para.h00[[ss]] <-
    as.numeric(Para.Grid.WX0.Target[which.min(CV.Curve),])
}

pos.h00.MM <- list()
h00.MM <- list()
h00.predict <- list()

for(ss in 1:1){
  
  pos.h00.MM[[ss]] <- SS.posA0[[ss]]
  Resp.Temp <- Y.MM[[ss]]*(1-A.MM[[ss]])*(1-D.MM[[ss]])
  
  h00.MM[[ss]] <- FT_PMMR( Y         =Resp.Temp[pos.h00.MM[[ss]]],
                           Perturb   =ZX0.MM[[ss]][pos.h00.MM[[ss]],],
                           Target    =WX0.MM[[ss]][pos.h00.MM[[ss]],],
                           Diagonal  =(1-A.MM[[ss]])[pos.h00.MM[[ss]]], 
                           Perturb.bw=exp(Opt.Para.h00[[ss]][1]),
                           Target.bw =exp(Opt.Para.h00[[ss]][2]),
                           lambda    =exp(Opt.Para.h00[[ss]][3]),
                           NV        =FALSE)
  
}


h00.predict[[1]] <- function(WX.New.Input){
  FT_RBF(X     = WX0.MM[[1]][pos.h00.MM[[1]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h00[[1]][2]))%*%h00.MM[[1]]$alpha + h00.MM[[1]]$intercept
}

h00.predict[[2]] <- function(WX.New.Input){
  FT_RBF(X     = WX0.MM[[2]][pos.h00.MM[[2]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h00[[2]][2]))%*%h00.MM[[2]]$alpha + h00.MM[[2]]$intercept
}

h00hat.NoCF <- h00hat.CF <- rep(0,N)
for(ss in 1:1){
  h00hat.NoCF[SS.Index[[ss]]] <- h00.predict[[ss]](WX0.MM[[ss]])
  h00hat.CF[SS.Index[[3-ss]]] <- h00.predict[[ss]](WX0.MM[[3-ss]])
}

h00hat.NoCF.MM <- list()
for(ss in 1:1){
  print(ss)
  h00hat.NoCF.MM[[ss]]  <- h00hat.NoCF[SS.Index[[ss]] ]
}

###################################################
# Bridge Function: h2
###################################################

Opt.Para.h2 <- list()
for(ss in 1:1){
  
  CV.Curve <- apply(Para.Grid.WX0.Target,
                    1,
                    FUN=function(vv){ CrossVar(vv,
                                               ss,
                                               subset=posA0,
                                               response=(1-D)*(1-A),
                                               target=WX0,
                                               perturb=ZX0,
                                               diagonal=(1-A)) })
  
  Opt.Para.h2[[ss]] <-
    as.numeric(Para.Grid.WX0.Target[which.min(CV.Curve),])
}

pos.h2.MM <- list()
h2.MM <- list()
h2.predict <- list()

for(ss in 1:1){
  
  pos.h2.MM[[ss]] <- SS.posA0[[ss]]
  Resp.Temp <- (1-D.MM[[ss]])*(1-A.MM[[ss]])
  
  h2.MM[[ss]] <- FT_PMMR( Y         =Resp.Temp[pos.h2.MM[[ss]]],
                          Perturb   =ZX0.MM[[ss]][pos.h2.MM[[ss]],],
                          Target    =WX0.MM[[ss]][pos.h2.MM[[ss]],],
                          Diagonal  =(1-A.MM[[ss]])[pos.h2.MM[[ss]]], 
                          Perturb.bw=exp(Opt.Para.h2[[ss]][1]),
                          Target.bw =exp(Opt.Para.h2[[ss]][2]),
                          lambda    =exp(Opt.Para.h2[[ss]][3]),
                          NV        =FALSE)
  
}

h2.predict[[1]] <- function(WX.New.Input){
  FT_RBF(X     = WX0.MM[[1]][pos.h2.MM[[1]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h2[[1]][2]))%*%h2.MM[[1]]$alpha + h2.MM[[1]]$intercept
}

h2.predict[[2]] <- function(WX.New.Input){
  FT_RBF(X     = WX0.MM[[2]][pos.h2.MM[[2]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h2[[2]][2]))%*%h2.MM[[2]]$alpha + h2.MM[[2]]$intercept
}

h2hat.NoCF <- h2hat.CF <- rep(0,N)
for(ss in 1:1){
  h2hat.NoCF[SS.Index[[ss]]] <- h2.predict[[ss]](WX0.MM[[ss]])
  h2hat.CF[SS.Index[[3-ss]]] <- h2.predict[[ss]](WX0.MM[[3-ss]])
}

h2hat.NoCF.MM <- list()
for(ss in 1:1){
  print(ss)
  h2hat.NoCF.MM[[ss]]  <- h2hat.NoCF[SS.Index[[ss]] ]
}

###################################################
# Bridge Function: q0
###################################################

Opt.Para.q0 <- list()
for(ss in 1:1){
  
  CV.Curve <- apply(Para.Grid.ZX0.Target,
                    1,
                    FUN=function(vv){ CrossVar(vv,
                                               ss,
                                               subset=(1:N),
                                               response=rep(1,N),
                                               target=ZX0,
                                               perturb=WX0,
                                               diagonal=(1-A)) })
  
  Opt.Para.q0[[ss]] <-
    as.numeric(Para.Grid.ZX0.Target[which.min(CV.Curve),])
}

pos.q0.MM <- list()
q0.MM <- list()
q0.predict <- list()

for(ss in 1:1){
  
  pos.q0.MM[[ss]] <- 1:(N/CF)
  
  Resp.Temp <- rep(1,N)
  
  q0.MM[[ss]] <- FT_PMMR( Y         =Resp.Temp[pos.q0.MM[[ss]]],
                          Perturb   =WX0.MM[[ss]][pos.q0.MM[[ss]],],
                          Target    =ZX0.MM[[ss]][pos.q0.MM[[ss]],],
                          Diagonal  =(1-A.MM[[ss]])[pos.q0.MM[[ss]]],
                          Perturb.bw=exp(Opt.Para.q0[[ss]][1]),
                          Target.bw =exp(Opt.Para.q0[[ss]][2]),
                          lambda    =exp(Opt.Para.q0[[ss]][3]),
                          NV        =FALSE)
}


q0.predict[[1]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX0.MM[[1]][pos.q0.MM[[1]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q0[[1]][2]))%*%q0.MM[[1]]$alpha + q0.MM[[1]]$intercept 
}

q0.predict[[2]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX0.MM[[2]][pos.q0.MM[[2]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q0[[2]][2]))%*%q0.MM[[2]]$alpha + q0.MM[[2]]$intercept
}

q0hat.NoCF <- q0hat.CF <- rep(0,N)
for(ss in 1:1){
  q0hat.NoCF[SS.Index[[ss]]] <- q0.predict[[ss]](ZX0.MM[[ss]])
  q0hat.CF[SS.Index[[3-ss]]] <- q0.predict[[ss]](ZX0.MM[[3-ss]])
}

q0hat.NoCF.MM <- list()
for(ss in 1:1){
  q0hat.NoCF.MM[[ss]] <- q0hat.NoCF[SS.Index[[ss]] ]
}

###################################################
# Bridge Function: q11
###################################################

Opt.Para.q11 <- list()
for(ss in 1:1){
  
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

for(ss in 1:1){
  
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
for(ss in 1:1){
  q11hat.NoCF[SS.Index[[ss]]] <- q11.predict[[ss]](ZX01.MM[[ss]])
  q11hat.CF[SS.Index[[3-ss]]] <- q11.predict[[ss]](ZX01.MM[[3-ss]])
}

q11hat.NoCF.MM <- list()
for(ss in 1:1){
  print(ss)
  q11hat.NoCF.MM[[ss]]  <- q11hat.NoCF[SS.Index[[ss]] ]
}

















#####################
# Generate
#####################

U.0  <- rnorm(N.0)
X0.0 <- matrix(rnorm(N.0*dimX0),N.0,dimX0)

PrA1.0 <- expit( betaa0 + 
                   betaax0*X0.0 +
                   betaau*U.0 )
PrA0.0 <- 1-PrA1.0

A.0 <- rbinom(N.0,1,PrA1.0)

Z.0  <- ( rnorm(N.0) + 
            betaz0 + 
            apply(betazx0*X0.0,1,sum) + 
            betazu*U.0 + 
            betaza*A.0 )
W.0  <- ( rnorm(N.0) + 
            betaw0 + 
            apply(betawx0*X0.0,1,sum) +
            betawu*U.0 )

X1A1.0 <- betax10.1+betax1x0.1*X0.0+betax1u.1*matrix(U.0,N.0,dimX0)+rnorm(N.0)
X1A0.0 <- betax10.0+betax1x0.0*X0.0+betax1u.0*matrix(U.0,N.0,dimX0)+rnorm(N.0) 

X1.0 <- A.0*X1A1.0 + (1-A.0)*X1A0.0

h11.true.vec <- h11.true(W.0,X0.0,X1.0)
h01.true.vec <- h01.true(W.0,X0.0)
h00.true.vec <- h00.true(W.0,X0.0)
h2.true.vec  <- h2.true(W.0,X0.0)
q0.true.vec  <- q0.true(Z.0,X0.0)
q11.true.vec <- q11.true(Z.0,X0.0,X1.0)

################################################################################
# Summary
################################################################################

W.0.Scale <- (W.0-MEANMAT[1,3])/SDMAT[1,3]
Z.0.Scale <- (Z.0-MEANMAT[1,2])/SDMAT[1,2]
X0.0.Scale <- (X0.0-matrix(MEANMAT[1,3+1:dimX0],N.0,dimX0,byrow=T))/
  matrix(SDMAT[1,3+1:dimX0],N.0,dimX0,byrow=T)
X1.0.Scale <- (X1.0-matrix(MEANMAT[1,3+dimX0+1:dimX0],N.0,dimX0,byrow=T))/
  matrix(SDMAT[1,3+dimX0+1:dimX0],N.0,dimX0,byrow=T)

h11hat.CF.Full <- h11.predict[[1]]( cbind(W.0.Scale,X0.0.Scale,X1.0.Scale) )
h01hat.CF.Full <- h01.predict[[1]]( cbind(W.0.Scale,X0.0.Scale) )
h00hat.CF.Full <- h00.predict[[1]]( cbind(W.0.Scale,X0.0.Scale) )
h2hat.CF.Full  <- h2.predict[[1]](  cbind(W.0.Scale,X0.0.Scale) )
q0hat.CF.Full  <- q0.predict[[1]](  cbind(Z.0.Scale,X0.0.Scale) )
q11hat.CF.Full <- q11.predict[[1]]( cbind(Z.0.Scale,X0.0.Scale,X1.0.Scale) )

h11.true.vec <- (h11.true.vec - MEANMAT[1,1])/SDMAT[1,1]
h01.true.vec <- (h01.true.vec - MEANMAT[1,1]*h2.true.vec)/SDMAT[1,1]
h00.true.vec <- (h00.true.vec - MEANMAT[1,1]*h2.true.vec)/SDMAT[1,1]

MSE.h11 <- mean(((h11hat.CF.Full- h11.true.vec)^2))
MSE.h01 <- mean(((h01hat.CF.Full- h01.true.vec)^2))
MSE.h00 <- mean(((h00hat.CF.Full- h00.true.vec)^2))
MSE.h2 <- mean(((h2hat.CF.Full- h2.true.vec)^2))
MSE.q0 <- mean(((q0hat.CF.Full- q0.true.vec)^2))
MSE.q11 <- mean(((q11hat.CF.Full- q11.true.vec)^2))

RESULT <- data.frame(matrix(c(N,seed.data,
                              MSE.h11,MSE.h01,MSE.h00,MSE.h2,MSE.q0,MSE.q11,
                              c(Opt.Para.h11[[1]],
                                Opt.Para.h01[[1]],
                                Opt.Para.h00[[1]],
                                Opt.Para.h2[[1]],
                                Opt.Para.q0[[1]],
                                Opt.Para.q11[[1]])),1,2+6+3*6))
colnames(RESULT) <- c("N","Seed",
                      "MSE.h11", "MSE.h01", "MSE.h00", "MSE.h2", "MSE.q0", "MSE.q11",
                      sprintf("Para.h11.%s",c("Pert","Main","Lambda")),
                      sprintf("Para.h01.%s",c("Pert","Main","Lambda")),
                      sprintf("Para.h00.%s",c("Pert","Main","Lambda")),
                      sprintf("Para.h2.%s",c("Pert","Main","Lambda")),
                      sprintf("Para.q0.%s",c("Pert","Main","Lambda")),
                      sprintf("Para.q11.%s",c("Pert","Main","Lambda")))

write.csv(RESULT,
          sprintf("PSE_Nuisance_N%0.2d_B%0.5d.csv",
                  nn,seed.data),
          row.names=F)

