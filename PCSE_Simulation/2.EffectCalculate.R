#####################
# Cluster computing: BATCH=1,...,1000
#####################

args <- (commandArgs(trailingOnly=TRUE))
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

expit <- function(v){exp(v)/(1+exp(v))}

N <- 10^6

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

Eff.Obs <- mean((Y.Ay1.Ad0 - Y.Ay0.Ad0)[D.Ad0==0])








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

Eff.Exp <- mean((Y.Ay1.Ad0 - Y.Ay0.Ad0)[D.Ad0==0])

write.csv(matrix(c(Eff.Obs,Eff.Exp),1,2),
          sprintf("Effect/Effect_B%0.5d.csv",BATCH),
          row.names=F)
 