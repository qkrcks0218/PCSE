################################################################################
# For Paper
################################################################################

Est.1 <- Est.0 <- Est <- SE.1 <- SE.0 <- SE <- rep(0,2)
Est.1.NoU <- Est.0.NoU <- Est.NoU <- SE.1.NoU <- SE.0.NoU <- SE.NoU <- rep(0,2)

RRR1.1 <- read.csv("Result_F/Result_PMMR_MICE.csv")
RRR1.0 <- read.csv("Result_F/Result_Ign_MICE.csv")
RRR0.1 <- read.csv("Result_F/Result_PMMR_FULL.csv")
RRR0.0 <- read.csv("Result_F/Result_Ign_FULL.csv")

# IND1 <- intersect( RRR1.1$Seed[!is.na(apply(RRR1.1,1,sum))] , 
#                    RRR1.0$Seed[!is.na(apply(RRR1.0,1,sum))] )
#   
# IND0 <- intersect( RRR0.1$Seed[!is.na(apply(RRR0.1,1,sum))] , 
#                    RRR0.0$Seed[!is.na(apply(RRR0.0,1,sum))] )
# 
# RRR1.1 <- RRR1.1[which( RRR1.1$Seed %in% sample(IND1,500) ),-1]
# RRR1.0 <- RRR1.0[which( RRR1.0$Seed %in% sample(IND1,500) ),-1]
# RRR0.1 <- RRR0.1[which( RRR0.1$Seed %in% sample(IND0,200) ),-1]
# RRR0.0 <- RRR0.0[which( RRR0.0$Seed %in% sample(IND0,200) ),-1]
# 
# write.csv(RRR1.1,"Result_F/Result_PMMR_MICE.csv",row.names=F)
# write.csv(RRR1.0,"Result_F/Result_Ign_MICE.csv",row.names=F)
# write.csv(RRR0.1,"Result_F/Result_PMMR_FULL.csv",row.names=F)
# write.csv(RRR0.0,"Result_F/Result_Ign_FULL.csv",row.names=F)

SUM <- function(RRR){
  EST <- apply(RRR[,c(3,5,1)],2,median)
  SE <- sqrt( apply((RRR[,c(3,5,1)] - 
                       matrix(apply(RRR[,c(3,5,1)],2,median),
                              nrow(RRR),3,byrow=T))^2 + 
                      RRR[,c(4,6,2)]^2,2,median) )
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

SUM(RRR0.1)
SUM(RRR0.0)

