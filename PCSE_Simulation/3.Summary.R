library(png)


addImg <- function( obj,  x = NULL,  y = NULL,  width = NULL,  interpolate = TRUE ){ 
  USR <- par()$usr 
  PIN <- par()$pin 
  DIM <- dim(obj) 
  ARp <- DIM[1]/DIM[2] 
  WIDi <- width/(USR[2]-USR[1])*PIN[1] 
  HEIi <- WIDi * ARp 
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) 
  rasterImage(image = obj, 
              xleft = x, xright = x+(width),
              ybottom = y-0.5*(HEIu), ytop = y+0.5*(HEIu), 
              interpolate = interpolate)
}

addImg2 <- function( obj,  x = NULL,  y = NULL,  height = NULL,  interpolate = TRUE ){ 
  USR <- par()$usr 
  PIN <- par()$pin 
  DIM <- dim(obj) 
  ARp <- DIM[1]/DIM[2] 
  HEIi <- height/(USR[4]-USR[3])*PIN[2]
  WIDi <- HEIi / ARp 
  WIDi <- WIDi/PIN[1]*(USR[2]-USR[1]) 
  rasterImage(image = obj, 
              xleft = x, xright = x+(WIDi),
              ybottom = y-0.5*(height), ytop = y+0.5*(height), 
              interpolate = interpolate)
}

# setwd("D:/Dropbox/Chan/Research/Postdoc2022/PSE_Simluation_May/DGPs/May6")
setwd("/Users/chanpark/Library/CloudStorage/Dropbox/Chan/Research/Postdoc2022/PSE_Simluation_May/NewSim_Final")

TAU1 <- readPNG("TAU1.png")
TAU2 <- readPNG("TAU2.png")
C1 <- readPNG("Caption1.png")
C2 <- readPNG("Caption2.png")

Estimand <- mean( apply(read.csv("Estimand.csv"),2,mean)[2:3] )

png("Simulation_May.png",height=6,width=8,unit="in",res=500)

MAR    <- c(1,0,0.5,0.5)
par(mar=MAR*c(0,1,0,1))
layout(rbind(c(7,7),c(5,6),c(1,3),c(8,8),c(2,4)),heights=c(1.5,2,10,1.5,8))

Eff.Obs <- Estimand
Eff.Exp <- Estimand
AXIS2 <- -7.5

WANTTOSEE <- "Eff"

N.Grid <- c(500,1000,1500,2000)

Est <- SE <- SE.B <- Est.NoU <- list()
FL   <- "Result" 

FILE <- c(sprintf("%s/Result_U_obs_N%s.csv", FL,N.Grid[1]),
          sprintf("%s/Result_U_obs_N%s.csv", FL,N.Grid[2]),
          sprintf("%s/Result_U_obs_N%s.csv", FL,N.Grid[3]),
          sprintf("%s/Result_U_obs_N%s.csv", FL,N.Grid[4]))

FILE.NoU <- c(sprintf("%s/Result_NoU_obs_N%s.csv", FL,N.Grid[1]),
              sprintf("%s/Result_NoU_obs_N%s.csv", FL,N.Grid[2]),
              sprintf("%s/Result_NoU_obs_N%s.csv", FL,N.Grid[3]),
              sprintf("%s/Result_NoU_obs_N%s.csv", FL,N.Grid[4]))

TT <- 1:4

N.Text <- c(sprintf("N=%s",N.Grid[TT]))

XL <- c(-0.5,max(TT)+0.5)
YL <- c(-1.25,1.25)
shift <- 0.25
boxwidth <- 0.25
COL <- list()
COL[[1]] <- rgb(0,0,0,0.5)
COL[[2]] <- rgb(1,0,0,0.5)
XP1 <- 1.2
XP2 <- 3.0
TW  <- 1.2

for(Est.T in c("Eff")){
  
  for(tt in TT){
    
    Result.temp <- read.csv(FILE[tt])
    Result.Est <- Result.temp$Est
    Result.SE <- Result.temp$SE
    Result.SE.B <- Result.temp$SE.B
    Est[[tt]] <- Result.Est
    SE[[tt]]  <- Result.SE
    SE.B[[tt]]  <- Result.SE.B
    
    Result.temp.NoU <- read.csv(FILE.NoU[tt])
    Est.NoU[[tt]] <- Result.temp.NoU$Est.No
    
  }
  
  par(mar=MAR*c(0,1,1,1))
  for(tt in TT[1]){
    boxplot(Est[[tt]]-(Eff.Obs),at=TT[1]-shift/2,
            xlab="",ylab="",axes=FALSE,
            xlim=XL,
            ylim=YL,
            medlwd = 1,
            boxwex = boxwidth,
            col=COL[[1]],
            pch=19,cex=0.2)
    boxplot(Est.NoU[[tt]]-(Eff.Obs),at=TT[1]+shift/2,
            xlab="",ylab="",axes=FALSE,add=TRUE,
            xlim=XL,
            ylim=YL,
            medlwd = 1,
            boxwex = boxwidth,
            col=COL[[2]],
            pch=19,cex=0.2)
    
  }
  if(length(TT)>1){
    for(tt in TT[-1]){
      boxplot(Est[[tt]]-(Eff.Obs),at=tt-shift/2,
              xlab="",ylab="",axes=FALSE,add=TRUE,
              xlim=XL,
              ylim=YL,
              medlwd = 1,
              boxwex = boxwidth,
              col=COL[[1]],
              pch=19,cex=0.2)
      boxplot(Est.NoU[[tt]]-(Eff.Obs),at=tt+shift/2,
              xlab="",ylab="",axes=FALSE,add=TRUE,
              xlim=XL,
              ylim=YL,
              medlwd = 1,
              boxwex = boxwidth,
              col=COL[[2]],
              pch=19,cex=0.2)
    }
    
  }
  # points(XP1,-1,col=COL[[1]],pch=15,cex=2)
  # points(XP2,-1,col=COL[[2]],pch=15,cex=2)
  # addImg(TAU1, 
  #        x = XP1+0.2, 
  #        y = -1, 
  #        width = TW)
  # addImg(TAU2, 
  #        x = XP2+0.2, 
  #        y = -1, 
  #        width = TW)
  
  axis(2,line=AXIS2,cex.axis=0.8)
  title(ylab="Bias",line=-4.5)
  segments(0.5,0,max(TT)+0.5,0,col=c(2),lty=2)
  axis(3,col=1,at=(1:8)[TT],labels=N.Text,line=-1,tick=F)
  
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
          ylim=c(10,11),bty="n",xlim=XL,
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



FILE <- c(sprintf("%s/Result_U_exp_N%s.csv", FL,N.Grid[1]),
          sprintf("%s/Result_U_exp_N%s.csv", FL,N.Grid[2]),
          sprintf("%s/Result_U_exp_N%s.csv", FL,N.Grid[3]),
          sprintf("%s/Result_U_exp_N%s.csv", FL,N.Grid[4]))

FILE.NoU <- c(sprintf("%s/Result_NoU_exp_N%s.csv", FL,N.Grid[1]),
              sprintf("%s/Result_NoU_exp_N%s.csv", FL,N.Grid[2]),
              sprintf("%s/Result_NoU_exp_N%s.csv", FL,N.Grid[3]),
              sprintf("%s/Result_NoU_exp_N%s.csv", FL,N.Grid[4]))

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
    
    Result.temp.NoU <- read.csv(FILE.NoU[tt])
    Est.NoU[[tt]] <- Result.temp.NoU$Est.No
    
  }
  
  par(mar=MAR*c(0,1,1,1))
  for(tt in TT[1]){
    boxplot(Est[[tt]]-(Eff.Obs),at=TT[1]-shift/2,
            xlab="",ylab="",axes=FALSE,
            xlim=XL,
            ylim=YL,
            medlwd = 1,
            boxwex = boxwidth,
            col=COL[[1]],
            pch=19,cex=0.2)
    boxplot(Est.NoU[[tt]]-(Eff.Obs),at=TT[1]+shift/2,
            xlab="",ylab="",axes=FALSE,add=TRUE,
            xlim=XL,
            ylim=YL,
            medlwd = 1,
            boxwex = boxwidth,
            col=COL[[2]],
            pch=19,cex=0.2)
    
  }
  if(length(TT)>1){
    for(tt in TT[-1]){
      boxplot(Est[[tt]]-(Eff.Obs),at=tt-shift/2,
              xlab="",ylab="",axes=FALSE,add=TRUE,
              xlim=XL,
              ylim=YL,
              medlwd = 1,
              boxwex = boxwidth,
              col=COL[[1]],
              pch=19,cex=0.2)
      boxplot(Est.NoU[[tt]]-(Eff.Obs),at=tt+shift/2,
              xlab="",ylab="",axes=FALSE,add=TRUE,
              xlim=XL,
              ylim=YL,
              medlwd = 1,
              boxwex = boxwidth,
              col=COL[[2]],
              pch=19,cex=0.2)
    }
    
  }
  # points(XP1,-1,col=COL[[1]],pch=15,cex=2)
  # points(XP2,-1,col=COL[[2]],pch=15,cex=2)
  # addImg(TAU1, 
  #        x = XP1+0.2, 
  #        y = -1, 
  #        width = TW)
  # addImg(TAU2, 
  #        x = XP2+0.2, 
  #        y = -1, 
  #        width = TW)
  
  axis(2,line=AXIS2,cex.axis=0.8)
  title(ylab="Bias",line=-4.5)
  segments(0.5,0,max(TT)+0.5,0,col=c(2),lty=2)
  axis(3,col=1,at=(1:8)[TT],labels=N.Text,line=-1,tick=F)
  
  axis(2,line=AXIS2,cex.axis=0.8)
  title(ylab="Bias",line=-4.5)
  segments(0.5,0,max(TT)+0.5,0,col=c(2),lty=2)
  axis(3,col=1,at=(1:8)[TT],labels=N.Text,line=-1,tick=F)
  
  
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
          ylim=c(10,11),bty="n",xlim=XL,
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
par(mar=MAR*c(1,1,1,1))
plot.new()
text(0.5,0.7,"Observational Setting",cex=1.4,font=1)
plot.new()
text(0.5,0.7,"Experimental Setting",cex=1.4,font=1)

FACTOR <- 0.8
MAR    <- c(1,0,1,0.5)
par(mar=MAR*c(0,1,0,1))
plot.new()
addImg2(C1, 
       x = -0.02, 
       y = 0.5, 
       height = 0.75*FACTOR)
points(0.46*FACTOR-0.02,0.5,col=COL[[1]],pch=15,cex=4)
points(0.865*FACTOR-0.02,0.5,col=COL[[2]],pch=15,cex=4)

MAR    <- c(1,0,1,0.5)
par(mar=MAR*c(0,1,0,1))
plot.new()
addImg2(C2, 
        x = -0.02, 
        y = 0.5, 
        height = 0.75*FACTOR)
dev.off()