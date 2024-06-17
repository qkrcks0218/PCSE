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

N <- 10^7
source("0.DGP.R")

Calculate.Estimand <- function(Sim.Setup.Obs){
  
  ################################################################################
  # Data Generation : Observational
  ################################################################################
  
  if(Sim.Setup.Obs){
    
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
    
    Data <- data.frame( cbind(Y,A,D,W,Z,X1,X2,U) )
    
    
  }
  
  ################################################################################
  # Data Generation : Experimental
  ################################################################################
  
  if(!Sim.Setup.Obs){
    
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
    X1 <- matrix(rnorm(dimX*N),N,dimX)
    X2 <- matrix(rnorm(dimX*N),N,dimX)
    
    Z  <- ( rnorm(N) + 
              betaz0 + 
              apply(betazx1*X1,1,sum) + 
              apply(betazx2*X2,1,sum) + 
              betazu*U )
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
    
    
    A <- rbinom(N,1,PrA1)
    D <- A*(D.Ad1)+(1-A)*(D.Ad0)
    Y <- A*(Y.Ay1.Ad1)+(1-A)*(Y.Ay0.Ad0)
    
    Data <- data.frame( cbind(Y,A,D,W,Z,X1,X2,U) )
    
    
  }
  
  
  
  return( mean( (Y.Ay1.Ad0 - Y.Ay0.Ad0)[D.Ad0==0] ) )
}

Eff.Obs <- Calculate.Estimand(TRUE)
Eff.Exp <- Calculate.Estimand(FALSE)

RESULT <- matrix(c(BATCH,Eff.Obs,Eff.Exp),1,3)
RESULT <- as.data.frame(RESULT)
colnames(RESULT) <- c("Batch","Eff_Obs","Eff_Exp")

write.csv(RESULT,
          sprintf("Estimand/Effect_B%0.5d.csv",BATCH),
                  row.names=F)


## After Merge All Files

RESULT <- read.csv("Estimand.csv")
apply(RESULT,2,mean)
# 1.982401    1.982148 


