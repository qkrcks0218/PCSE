################################################################################
# For Paper
################################################################################

Est.1 <- Est.0 <- Est <- SE.1 <- SE.0 <- SE <- rep(0,2)
Est.1.NoU <- Est.0.NoU <- Est.NoU <- SE.1.NoU <- SE.0.NoU <- SE.NoU <- rep(0,2)

RRR1.1 <- read.csv("Result/RESULT_PMMR_Merge.csv")
RRR1.0 <- read.csv("Result/RESULT_Ign_Merge.csv")

nrow(RRR1.1);nrow(RRR1.0)

SUM <- function(RRR){
  EST <- apply(RRR[,1+c(3,5,1)],2,median)
  SE <- sqrt( apply((RRR[,1+c(3,5,1)] - 
                       matrix(apply(RRR[,1+c(3,5,1)],2,median),
                              nrow(RRR),3,byrow=T))^2 + 
                      RRR[,1+c(4,6,2)]^2,2,median) )
  LIST <- list()
  LIST$EST <- round(EST,2)
  LIST$SE  <- round(SE,2)
  LIST$LB <- round(EST-1.96*SE,2)
  LIST$UB <- round(EST+1.96*SE,2)
  RR <- cbind(LIST$EST,LIST$SE,LIST$LB,LIST$UB)
  
  return(RR)
}

SUM(RRR1.1)
SUM(RRR1.0)


