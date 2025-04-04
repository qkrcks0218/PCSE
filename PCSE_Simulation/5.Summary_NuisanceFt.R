grid <- (2:16)

ERROR <- matrix(0,15,8)
VAR <- matrix(0,15,8)
XX <- c(8:16)
pp <- XX-1

for(ii in pp){
  
  RRR <- read.csv(sprintf("NF/PSE_Merged_Nuisance_N%0.2d.csv",grid[ii])) 
  ERROR[ii,] <- c(apply( RRR,2,mean ),nrow(RRR))
  print(nrow(RRR))
}

ERROR <- ERROR[pp,]
ERROR <- cbind(ERROR,XX/2)

r.h11 <- coef(lm(I(0.5*log(ERROR[,2]))~ERROR[,ncol(ERROR)])) 
r.h01 <- coef(lm(I(0.5*log(ERROR[,3]))~ERROR[,ncol(ERROR)])) 
r.h00 <- coef(lm(I(0.5*log(ERROR[,4]))~ERROR[,ncol(ERROR)])) 
r.h2  <- coef(lm(I(0.5*log(ERROR[,5]))~ERROR[,ncol(ERROR)])) 
r.q0  <- coef(lm(I(0.5*log(ERROR[,6]))~ERROR[,ncol(ERROR)])) 
r.q11 <- coef(lm(I(0.5*log(ERROR[,7]))~ERROR[,ncol(ERROR)])) 


r.h11[2]+r.q11[2]
r.h11[2]+r.q0 [2]
r.h01[2]+r.q0 [2]
r.h00[2]+r.q0 [2]
r.h2 [2]+r.q0 [2] 

png("ConvergenceRate.png",height=4,width=8,unit="in",res=500)

LL <- matrix(1:6,2,3,byrow=T)
LL <- rbind(LL,c(7,7,7))
LL <- cbind(c(8,8,9),LL)
MAR <- c(2.5,2.5,2,1)

layout(LL,widths=c(0.5,4,4,4),heights=c(4,4,0.5))
par(mar=MAR)
TITLE <- c(expression(h["1"]),
           expression(h["0"]*"(a=1)"),
           expression(h["0"]*"(a=0)"),
           expression(h["2"]),
           expression(q["0"]),
           expression(q["1"]))
YL <- matrix(0,6,2)
YL[1,] <- c(-3.0,-1.4)
YL[2,] <- c(-3.0,-1.7)
YL[3,] <- c(-3.1,-2.0)
YL[4,] <- c(-3.8,-2.2)
YL[5,] <- c(-2.7,-1.2)
YL[6,] <- c(-1.65,-0.4)

COEF <- rbind(r.h11,r.h01,r.h00,r.h2,r.q0,r.q11)

for (tt in 1:6){
  plot(XX/2,0.5*log(ERROR[,tt+1]),
       ylim=YL[tt,],
       cex=1, pch=4, lwd=2,
       xaxt='n')
  axis(1, at=c(2,4,6,8,10,12,14,16)/2, label=c(50*2^(0:7)))
  title(main=TITLE[tt], cex.main=1.5)
  abline(a=COEF[tt,1],b=COEF[tt,2],col=2,lwd=1.25,lty=1)
  
  points(XX/2,0.5*log(ERROR[,tt+1]), 
         cex=1, pch=4, lwd=2,
         xaxt='n')
  
  text(mean(XX)/2+1,YL[tt,1]*0.25+YL[tt,2]*0.75,
       sprintf("Slope = %0.3f",COEF[tt,2]),col=2,pos=3,cex=1.25)
}

par(mar=c(0,1,0,1)*MAR)
plot.new()
text(0.5,0.5,"N",cex=1.5)

par(mar=c(1,0,1,0)*MAR)
plot.new()
text(0.5,0.5,"log ( RMSE )",cex=1.5,srt=90)

dev.off()



