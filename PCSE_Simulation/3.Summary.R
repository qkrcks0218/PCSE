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

N.Grid <- seq(500,2000,by=500)

FL   <- "Result"  
sdU  <- 0
if(sdU==0){
  Figure <- "Simulation_NoU.png"
} else {
  Figure <- "Simulation_U.png"
}

C1 <- readPNG(sprintf("Caption1_%s.png",sdU))
C2 <- readPNG("Caption2.png")

Estimand <- c(2,2)

Eff.Obs <- Estimand
Eff.Exp <- Estimand
AXIS2 <- -7.5

WANTTOSEE <- "Eff"

TT <- 1:4

Est <- SE <- SE.B <- 
  Est.NoU <- SE.NoU <- SE.NoU.B <- list()

N.Text <- c(sprintf("N=%s",N.Grid[TT]))

XL <- c(-0.5,max(TT)+0.5)
YL <- c(-0.55,0.55)
shift <- 0.25
boxwidth <- 0.25
COL <- list()
COL[[1]] <- rgb(0,0,0,0.5)
COL[[2]] <- rgb(0,0,0,0.1)
COL[[3]] <- rgb(0.5,0.5,0.5)
XP1 <- 1.2
XP2 <- 3.0
TW  <- 1.2

png(Figure,height=4,width=9,unit="in",res=500)

MAR    <- c(1,0,0.5,0.5)
par(mar=MAR*c(0,1,0,1))
layout(rbind(c(7,7),
             c(5,6),c(1,3),
             c(8,8),c(2,4)),
       heights=c(1.5, 2, 10, 1.5, 6))

for(sim.type in 1:2){
  
  if(sim.type==1){
    ST <- "obs"
  } else {
    ST <- "exp"
  }
  
  FILE <- c(sprintf("%s/Result_Merge_PMMR_%s_N%s_sdU%d.csv", FL, ST, N.Grid[1], sdU),
            sprintf("%s/Result_Merge_PMMR_%s_N%s_sdU%d.csv", FL, ST, N.Grid[2], sdU),
            sprintf("%s/Result_Merge_PMMR_%s_N%s_sdU%d.csv", FL, ST, N.Grid[3], sdU),
            sprintf("%s/Result_Merge_PMMR_%s_N%s_sdU%d.csv", FL, ST, N.Grid[4], sdU))
  
  FILE.NoU <- c(sprintf("%s/Result_Merge_Ign_%s_N%s_sdU%d.csv", FL, ST, N.Grid[1], sdU),
                sprintf("%s/Result_Merge_Ign_%s_N%s_sdU%d.csv", FL, ST, N.Grid[2], sdU),
                sprintf("%s/Result_Merge_Ign_%s_N%s_sdU%d.csv", FL, ST, N.Grid[3], sdU),
                sprintf("%s/Result_Merge_Ign_%s_N%s_sdU%d.csv", FL, ST, N.Grid[4], sdU))
  
  
  
  for(Est.T in c("Eff")){
    
    for(tt in TT){
      
      Result.temp <- read.csv(FILE[tt])
      Result.temp.NoU <- read.csv(FILE.NoU[tt])
      
      print(nrow(Result.temp))
      print(nrow(Result.temp.NoU))
      
      Result.Est <- Result.temp$Est
      Result.SE <- Result.temp$SE
      Result.SE.B <- Result.temp$SE.B
      Est[[tt]] <- Result.Est
      SE[[tt]]  <- Result.SE
      SE.B[[tt]]  <- Result.SE.B
      
      Est.NoU[[tt]] <- Result.temp.NoU$Est.No
      SE.NoU[[tt]]  <- Result.temp.NoU$SE.No
      SE.NoU.B[[tt]]  <- Result.temp.NoU$SE.B.No
      
      # Result.temp.NoU.OnlyX <- read.csv(FILE.NoU.OnlyX[tt])
      # Est.NoU.OnlyX[[tt]] <- Result.temp.NoU.OnlyX$Est.No
      
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
    
    axis(2,line=AXIS2,at=seq(-0.5,0.5,by=0.25),cex.axis=0.8)
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
    text(0.4,10+YP[1],"Bias (x10)",pos=2,cex=0.85)
    text(0.4,10+YP[2],"ESE (x10)",pos=2,cex=0.85)
    text(0.4,10+YP[3],"ASE (x10)",pos=2,cex=0.85)
    text(0.4,10+YP[4],"BSE (x10)",pos=2,cex=0.85)
    text(0.4,10+YP[5],"Coverage (ASE)",pos=2,cex=0.85)
    text(0.4,10+YP[6],"Coverage (BSE)",pos=2,cex=0.85)
    # text(0,10+YP[7],"Repetition")
    for(tt in TT){
      vvv <- sapply(tt:tt,Summary)
      text(tt-0.25,10+YP[1],sprintf("%0.3f",vvv[1]),cex=0.85)
      text(tt-0.25,10+YP[2],sprintf("%0.3f",vvv[2]),cex=0.85)
      text(tt-0.25,10+YP[3],sprintf("%0.3f",vvv[3]),cex=0.85)
      text(tt-0.25,10+YP[4],sprintf("%0.3f",vvv[4]),cex=0.85)
      text(tt-0.25,10+YP[5],sprintf("%0.3f",vvv[5]),cex=0.85)
      text(tt-0.25,10+YP[6],sprintf("%0.3f",vvv[6]),cex=0.85)
      # text(tt,10+YP[7],sprintf("%d",vvv[7]))
    }
    
    Summary <- function(tt){
      c(mean( Est.NoU[[tt]] - (Eff.Obs) )*10,
        sd( Est.NoU[[tt]] )*10,
        mean( SE.NoU[[tt]] )*10,
        mean( SE.NoU.B[[tt]] )*10,
        mean( abs( Est.NoU[[tt]] - (Eff.Obs) )/SE.NoU[[tt]] <= qnorm(0.975) ),
        mean( abs( Est.NoU[[tt]] - (Eff.Obs) )/SE.NoU.B[[tt]] <= qnorm(0.975) ),
        length(Est.NoU[[tt]])
      )
    } 
    
    for(tt in TT){
      vvv <- sapply(tt:tt,Summary)
      text(tt+0.25,10+YP[1],sprintf("%0.3f",vvv[1]),cex=0.85,col=COL[[3]])
      text(tt+0.25,10+YP[2],sprintf("%0.3f",vvv[2]),cex=0.85,col=COL[[3]])
      text(tt+0.25,10+YP[3],sprintf("%0.3f",vvv[3]),cex=0.85,col=COL[[3]])
      text(tt+0.25,10+YP[4],sprintf("%0.3f",vvv[4]),cex=0.85,col=COL[[3]])
      text(tt+0.25,10+YP[5],sprintf("%0.3f",vvv[5]),cex=0.85,col=COL[[3]])
      text(tt+0.25,10+YP[6],sprintf("%0.3f",vvv[6]),cex=0.85,col=COL[[3]])
      # text(tt,10+YP[7],sprintf("%d",vvv[7]))
    }
    
    
    
  }
  
  
  
}


MAR    <- c(1,0,1,0.5)
par(mar=MAR*c(0,1,0,1))
plot.new()
text(0.5,0.6,"Observational Setting",cex=1.2,font=1)
plot.new()
text(0.5,0.6,"Experimental Setting",cex=1.2,font=1)

FACTOR <- 0.8
MAR    <- c(1,0,1,0.5)
par(mar=MAR*c(0,1,0,1))
plot.new()
addImg2(C1, 
        x = -0.02, 
        y = 0.5, 
        height = 0.8*FACTOR)
# points(0.013+0.3125*FACTOR-0.06,0.5,col=COL[[1]],pch=15,cex=4)
# points(0.014+0.5775*FACTOR-0.092,0.5,col=COL[[2]],pch=15,cex=4)

MAR    <- c(1,0,1,0.5)
par(mar=MAR*c(0,1,0,1))
plot.new()
addImg2(C2, 
        x = -0.02, 
        y = 0.4, 
        height = 0.75*FACTOR)


dev.off()
