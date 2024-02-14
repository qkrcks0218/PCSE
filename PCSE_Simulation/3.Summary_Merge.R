################################################################################
# Aggregate files in Effect and Result Folders
################################################################################

LIST <- list.files("./Effect")
RRR <- read.csv( sprintf("./Effect/%s",LIST[1]) )
for(tt in 2:length(LIST)){
  RRR <- rbind(RRR, read.csv( sprintf("./Effect/%s",LIST[tt]) ))
}
# write.csv(RRR,"FinalData/TrueEffect.csv",row.names=F)

LIST <- list.files("./Result")
for(type in c("exp","obs")){
  for(nn in c("0500","1000","1500","2000")){
    POS <- which(substr(LIST,13,15)==type & substr(LIST,18,21)==nn)
    LIST.SHORT <- LIST[POS]
    RRR <- read.csv( sprintf("./Result/%s",LIST.SHORT[1]) )
    for(tt in 2:length(LIST.SHORT)){
      RRR <- rbind(RRR, read.csv( sprintf("./Result/%s",LIST.SHORT[tt]) ))
    }
    # write.csv(RRR,sprintf("FinalData/Result_PMMR_%s_N%s.csv",type,nn),row.names=F)
  }
}

################################################################################
# Aggregate files in Effect and Result Folders
################################################################################

png(file="Simulation.png",
    height=4*0.9,
    width=10*0.9,
    unit="in",
    res=500)

MAR    <- c(1,0,0.5,0.5)
par(mar=MAR*c(0,1,0,1))
layout(rbind(c(5,6),matrix(1:4,2,2)),heights=c(2,10,8))

Eff.Obs <- mean(read.csv("FinalData/TrueEffect.csv")[,1])
Eff.Exp <- mean(read.csv("FinalData/TrueEffect.csv")[,2])
AXIS2 <- -7.5

WANTTOSEE <- "Eff"

N.Grid <- c(500,1000,1500,2000)

Est <- SE <- SE.B <- list()
FL   <- "FinalData"


FILE <- c(sprintf("%s/Result_PMMR_obs_N%0.4d.csv", FL,N.Grid[1]),
          sprintf("%s/Result_PMMR_obs_N%0.4d.csv", FL,N.Grid[2]),
          sprintf("%s/Result_PMMR_obs_N%0.4d.csv", FL,N.Grid[3]),
          sprintf("%s/Result_PMMR_obs_N%0.4d.csv", FL,N.Grid[4]))
TT <- 1:4

N.Text <- c(sprintf("N=%s",N.Grid[TT]))

for(Est.T in c("Eff")){
  
  for(tt in TT){
    
    Result.temp <- read.csv(FILE[tt])
    Result.Est <- Result.temp$Est
    Result.SE <- Result.temp$SE
    Result.SE.B <- Result.temp$SE.B
    Est[[tt]] <- Result.Est
    SE[[tt]]  <- Result.SE
    SE.B[[tt]]  <- Result.SE.B
    
  }
  
  par(mar=MAR*c(0,1,1,1))
  for(tt in TT[1]){
    boxplot(Est[[tt]]-(Eff.Obs),at=TT[1],
            xlab="",ylab="",axes=FALSE,
            xlim=c(-0.5,max(TT)+0.5),
            ylim=c(-2.5,2.5),
            medlwd = 1,
            pch=19,cex=0.2)
    
  }
  if(length(TT)>1){
    for(tt in TT[-1]){
      boxplot(Est[[tt]]-(Eff.Obs),at=tt,
              xlab="",ylab="",axes=FALSE,add=TRUE,
              xlim=c(-0.5,max(TT)+0.5),
              ylim=c(-2.5,2.5),
              medlwd = 1,
              pch=19,cex=0.2)
    }
    
  }
  axis(2,line=AXIS2,cex.axis=0.8)
  title(ylab="Bias",line=-4.5)
  segments(0.5,0,max(TT)+0.5,0,col=c(2),lty=2)
  axis(3,col=1,at=(1:8)[TT],labels=N.Text,line=-0.5,tick=F)
  
  
  Summary <- function(tt){
    c(mean( Est[[tt]] - (Eff.Obs) )*10,
      sd( Est[[tt]] )*10,
      mean( SE[[tt]] )*10,
      mean( SE.B[[tt]] )*10,
      mean( abs( Est[[tt]] - (Eff.Obs) )/SE[[tt]] <= qnorm(0.975) ),
      mean( abs( Est[[tt]] - (Eff.Obs) )/SE.B[[tt]] <= qnorm(0.975) ),
      length(Est[[tt]])
    )
  }
  
  YP <- seq(1,0,length=8)[2:7]
  
  boxplot(cbind(rep(0,2),rep(0,2),rep(0,2),rep(0,2)),
          ylim=c(10,11),bty="n",xlim=c(-0.5,max(TT)+0.5),
          axes=F)
  text(0,10+YP[1],"Bias (x10)")
  text(0,10+YP[2],"ESE (x10)")
  text(0,10+YP[3],"ASE (x10)")
  text(0,10+YP[4],"BSE (x10)")
  text(0,10+YP[5],"Coverage (ASE)")
  text(0,10+YP[6],"Coverage (BSE)")
  # text(0,10+YP[7],"Repetition")
  for(tt in TT){
    vvv <- sapply(tt:tt,Summary)
    text(tt,10+YP[1],sprintf("%0.3f",vvv[1]))
    text(tt,10+YP[2],sprintf("%0.3f",vvv[2]))
    text(tt,10+YP[3],sprintf("%0.3f",vvv[3]))
    text(tt,10+YP[4],sprintf("%0.3f",vvv[4]))
    text(tt,10+YP[5],sprintf("%0.3f",vvv[5]))
    text(tt,10+YP[6],sprintf("%0.3f",vvv[6]))
    # text(tt,10+YP[7],sprintf("%d",vvv[7]))
  }
  
  
  
}


FILE <- c(sprintf("%s/Result_PMMR_exp_N%0.4d.csv", FL,N.Grid[1]),
          sprintf("%s/Result_PMMR_exp_N%0.4d.csv", FL,N.Grid[2]),
          sprintf("%s/Result_PMMR_exp_N%0.4d.csv", FL,N.Grid[3]),
          sprintf("%s/Result_PMMR_exp_N%0.4d.csv", FL,N.Grid[4]))
TT <- 1:4

N.Text <- c(sprintf("N=%s",N.Grid[TT]))

for(Est.T in c("Eff")){
  
  for(tt in TT){
    
    Result.temp <- read.csv(FILE[tt])
    Result.Est <- Result.temp$Est
    Result.SE <- Result.temp$SE
    Result.SE.B <- Result.temp$SE.B
    Est[[tt]] <- Result.Est
    SE[[tt]]  <- Result.SE
    SE.B[[tt]]  <- Result.SE.B
    
  }
  
  par(mar=MAR*c(0,1,1,1))
  for(tt in TT[1]){
    boxplot(Est[[tt]]-(Eff.Exp),at=TT[1],
            xlab="",ylab="",axes=FALSE,
            xlim=c(-0.5,max(TT)+0.5),
            ylim=c(-2.5,2.5),
            medlwd = 1,
            pch=19,cex=0.2)
    
  }
  if(length(TT)>1){
    for(tt in TT[-1]){
      boxplot(Est[[tt]]-(Eff.Exp),at=tt,
              xlab="",ylab="",axes=FALSE,add=TRUE,
              xlim=c(-0.5,max(TT)+0.5),
              ylim=c(-2.5,2.5),
              medlwd = 1,
              pch=19,cex=0.2)
    }
    
  }
  axis(2,line=AXIS2,cex.axis=0.8)
  title(ylab="Bias",line=-4.5)
  segments(0.5,0,max(TT)+0.5,0,col=c(2),lty=2)
  axis(3,col=1,at=(1:8)[TT],labels=N.Text,line=-0.5,tick=F)
  
  
  Summary <- function(tt){
    c(mean( Est[[tt]] - (Eff.Exp) )*10,
      sd( Est[[tt]] )*10,
      mean( SE[[tt]] )*10,
      mean( SE.B[[tt]] )*10,
      mean( abs( Est[[tt]] - (Eff.Exp) )/SE[[tt]] <= qnorm(0.975) ),
      mean( abs( Est[[tt]] - (Eff.Exp) )/SE.B[[tt]] <= qnorm(0.975) ),
      length(Est[[tt]])
    )
  }
  
  YP <- seq(1,0,length=8)[2:7]
  
  boxplot(cbind(rep(0,2),rep(0,2),rep(0,2),rep(0,2)),
          ylim=c(10,11),bty="n",xlim=c(-0.5,max(TT)+0.5),
          axes=F)
  text(0,10+YP[1],"Bias (x10)")
  text(0,10+YP[2],"ESE (x10)")
  text(0,10+YP[3],"ASE (x10)")
  text(0,10+YP[4],"BSE (x10)")
  text(0,10+YP[5],"Coverage (ASE)")
  text(0,10+YP[6],"Coverage (BSE)")
  # text(0,10+YP[7],"Repetition")
  for(tt in TT){
    vvv <- sapply(tt:tt,Summary)
    text(tt,10+YP[1],sprintf("%0.3f",vvv[1]))
    text(tt,10+YP[2],sprintf("%0.3f",vvv[2]))
    text(tt,10+YP[3],sprintf("%0.3f",vvv[3]))
    text(tt,10+YP[4],sprintf("%0.3f",vvv[4]))
    text(tt,10+YP[5],sprintf("%0.3f",vvv[5]))
    text(tt,10+YP[6],sprintf("%0.3f",vvv[6]))
    # text(tt,10+YP[7],sprintf("%d",vvv[7]))
  }
  
  
  
}




MAR    <- c(1,0,1,0.5)
par(mar=MAR*c(0,1,0,1))
plot.new()
text(0.5,0.8,"Observational Setting",cex=1.4,font=1)
plot.new()
text(0.5,0.8,"Experimental Setting",cex=1.4,font=1)

dev.off()

