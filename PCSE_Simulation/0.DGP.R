expit <- function(v){1/(1+exp(-v))}

################################################################################
# Parameters
################################################################################

#####################
# DGP
#####################

dimX <- 3

betay0.1  <- 4
betayx1.1 <- 0.3
betayx2.1 <- 0.4
betayw.1  <- 0.35
betayu.1  <- 0.5

betay0.0  <- 2
betayx1.0 <- 0.5
betayx2.0 <- 0.2
betayw.0  <- 0.5
betayu.0  <- 0.5

betaz0  <- 0
betazx1 <- 0.4
betazx2 <- 0.2
betazu  <- -1.5
betaza  <- -0.5

betaw0  <- 0
betawx1 <- 0.2
betawx2 <- 0.4
betawu  <- 1

betaa0  <- 0
betaax1 <- 0.2
betaax2 <- 0.1
betaau  <- 0.2

betad0.0  <- -0.6  # lower -> PrD0|A0 becomes smaller
betadx1.0 <- 0.01
betadx2.0 <- 0.01
betadu.0  <- 0.1
betadz.0  <- 0.01

betad0.1  <- -0.8  # lower -> PrD0|A1 becomes smaller
betadx1.1 <- 0.01
betadx2.1 <- 0.01
betadu.1  <- -0.05
betadz.1  <- 0.01

deltad0  <-  betad0.0  - betad0.1
deltadx1 <-  betadx1.0 - betadx1.1
deltadx2 <-  betadx2.0 - betadx2.1
deltadu  <-  betadu.0  - betadu.1

betaxu  <- 0

PrA1 <- 0.5
PrA0 <- 1-PrA1
