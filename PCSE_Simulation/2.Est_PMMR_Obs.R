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
FOLDER <- "WithU/"
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
PL <- -8;      PU <- 0
BW.L <- 3;     BW.U <- 4

Para.Grid <- expand.grid(seq(BW.L,BW.U,by=0.25),
                         seq(PL,PU,by=2))
Para.Grid <- cbind(Para.Grid[,1],Para.Grid)



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
  
  P.P <- rep(as.numeric(PARAMETER[1]),1+2*dimX)
  P.T <- rep(as.numeric(PARAMETER[2]),1+2*dimX)
  Lambda  <- as.numeric(PARAMETER[3])
  
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
 
#####################
# Generate
#####################

U    <- rnorm(N)
# X1.1 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U
# X1.2 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
# X1.3 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U
# X1.4 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
# X1.5 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
# X2.1 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U
# X2.2 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
# X2.3 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
# X2.4 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
# X2.5 <- sqrt(1-betaxu^2)*rnorm(N) + betaxu*U 
# X1   <- cbind(X1.1,X1.2,X1.3,X1.4,X1.5)
# X2   <- cbind(X2.1,X2.2,X2.3,X2.4,X2.5)
X1 <- data.frame(matrix(rnorm(dimX*N),N,dimX))
X2 <- data.frame(matrix(rnorm(dimX*N),N,dimX))

colnames(X1) <- sprintf("X1.%0.2d",1:dimX)
colnames(X2) <- sprintf("X1.%0.2d",1:dimX)

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
if(dimX==1){
  X1 <- matrix(X1,N,dimX)
  X2 <- matrix(X2,N,dimX)
}
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
  X2.MM [[ss]]  <- X[SS.Index[[ss]],1+1:dimX]
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
             Perturb.bw=exp(Opt.Para.h11[[ss]][1]),
             Target.bw =exp(Opt.Para.h11[[ss]][2]),
             lambda    =exp(Opt.Para.h11[[ss]][3]),
             NV        =FALSE)
}

h11.predict[[1]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[1]][pos.h11.MM[[1]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h11[[1]][2]))%*%h11.MM[[1]]$alpha + h11.MM[[1]]$intercept
}

h11.predict[[2]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[2]][pos.h11.MM[[2]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h11[[2]][2]))%*%h11.MM[[2]]$alpha + h11.MM[[2]]$intercept
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
             Perturb.bw=exp(Opt.Para.h01[[ss]][1]),
             Target.bw =exp(Opt.Para.h01[[ss]][2]),
             lambda    =exp(Opt.Para.h01[[ss]][3]),
             NV        =FALSE)
  
}


h01.predict[[1]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[1]][pos.h01.MM[[1]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h01[[1]][2]))%*%h01.MM[[1]]$alpha + h01.MM[[1]]$intercept
}

h01.predict[[2]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[2]][pos.h01.MM[[2]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h01[[2]][2]))%*%h01.MM[[2]]$alpha + h01.MM[[2]]$intercept
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

###################################################
# Bridge Function: h00
###################################################


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
                           Perturb.bw=exp(Opt.Para.h00[[ss]][1]),
                           Target.bw =exp(Opt.Para.h00[[ss]][2]),
                           lambda    =exp(Opt.Para.h00[[ss]][3]),
                           NV        =FALSE)
  
}


h00.predict[[1]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[1]][pos.h00.MM[[1]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h00[[1]][2]))%*%h00.MM[[1]]$alpha + h00.MM[[1]]$intercept
}

h00.predict[[2]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[2]][pos.h00.MM[[2]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h00[[2]][2]))%*%h00.MM[[2]]$alpha + h00.MM[[2]]$intercept
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

###################################################
# Bridge Function: h2
###################################################

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
                          Perturb.bw=exp(Opt.Para.h2[[ss]][1]),
                          Target.bw =exp(Opt.Para.h2[[ss]][2]),
                          lambda    =exp(Opt.Para.h2[[ss]][3]),
                          NV        =FALSE)
  
}

h2.predict[[1]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[1]][pos.h2.MM[[1]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h2[[1]][2]))%*%h2.MM[[1]]$alpha + h2.MM[[1]]$intercept
}

h2.predict[[2]] <- function(WX.New.Input){
  FT_RBF(X     = WX.MM[[2]][pos.h2.MM[[2]],],
         X.new = WX.New.Input,
         bw.median = exp(Opt.Para.h2[[2]][2]))%*%h2.MM[[2]]$alpha + h2.MM[[2]]$intercept
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

###################################################
# Bridge Function: q0
###################################################

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
                          Perturb.bw=exp(Opt.Para.q0[[ss]][1]),
                          Target.bw =exp(Opt.Para.q0[[ss]][2]),
                          lambda    =exp(Opt.Para.q0[[ss]][3]),
                          NV        =FALSE)
}


q0.predict[[1]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX.MM[[1]][pos.q0.MM[[1]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q0[[1]][2]))%*%q0.MM[[1]]$alpha + q0.MM[[1]]$intercept 
}

q0.predict[[2]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX.MM[[2]][pos.q0.MM[[2]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q0[[2]][2]))%*%q0.MM[[2]]$alpha + q0.MM[[2]]$intercept
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
                           Perturb.bw=exp(Opt.Para.q11[[ss]][1]),
                           Target.bw =exp(Opt.Para.q11[[ss]][2]),
                           lambda    =exp(Opt.Para.q11[[ss]][3]),
                           NV        =FALSE)
}

q11.predict[[1]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX.MM[[1]][pos.q11.MM[[1]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q11[[1]][2]))%*%q11.MM[[1]]$alpha + q11.MM[[1]]$intercept
}

q11.predict[[2]] <- function(ZX.New.Input){
  FT_RBF(X     = ZX.MM[[2]][pos.q11.MM[[2]],],
         X.new = ZX.New.Input,
         bw.median = exp(Opt.Para.q11[[2]][2]))%*%q11.MM[[2]]$alpha + q11.MM[[2]]$intercept
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

## Estimated

IF.Numer.1 <- ((A)*(1-D)*q11hat.CF*(Y-h11hat.CF)+(1-A)*q0hat.CF*((1-D)*h11hat.CF-h01hat.CF)+h01hat.CF)*SDMAT[1,1]
IF.Numer.0 <- ((1-A)*q0hat.CF*((1-D)*Y-h00hat.CF)+h00hat.CF)*SDMAT[1,1]
IF.Denom   <- (1-A)*q0hat.CF*((1-D)-h2hat.CF)+h2hat.CF

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

contrast   <- matrix(c(0,0,1/Est.Denom,-Est/Est.Denom),1,4)
SE         <- sqrt( contrast%*%Sigma%*%t(contrast) )


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
  vv <- mean( rnorm(N)*IF.Est )
  BOOT[bb] <- vv
  # vv <- mean( rnorm(N)*IF.Est.NoU )
  # BOOT.NoU[bb] <- vv
  vv <- mean( rnorm(N)*IF.Est.True )
  BOOT.True[bb] <- vv
}

Result1 <- c(
  seed.data,
  Est, SE, sd(BOOT),
  Est.Numer.1, Est.Numer.0, Est.Denom,
  
  # Est.NoU, SE.NoU, sd(BOOT.NoU), 
  # Est.Numer.1.NoU, Est.Numer.0.NoU, Est.Denom.NoU,
  
  Est.True, SE.True, sd(BOOT.True), 
  Est.Numer.1.True, Est.Numer.0.True, Est.Denom.True
)

Result1 <- matrix(Result1,1,length(Result1))

colnames(Result1) <- 
  c("Seed",
    "Est",    "SE",    "SE.B",
    "Est.Numer.1", "Est.Numer.0", "Est.Denom",
    
    # "Est.No",    "SE.No",    "SE.B.No",
    # "Est.Numer.1.No", "Est.Numer.0.No", "Est.Denom.No",
    
    "Est.True",    "SE.True",    "SE.B.True",
    "Est.Numer.1.True", "Est.Numer.0.True", "Est.Denom.True"
  )


write.csv(Result1,
          sprintf("%sResult_PMMR_obs_N%s_B%0.5d.csv",
                  FOLDER,
                  N, 
                  seed.data),
          row.names=F)

