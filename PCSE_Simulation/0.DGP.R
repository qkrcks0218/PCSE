expit <- function(v){1/(1+exp(-v))}

################################################################################
# Parameters
################################################################################

#####################
# DGP
#####################

dimX0 <- 5
dimX1 <- dimX0

betax10.1 <- 0.015
betax1x0.1 <- 0.02
betax1u.1 <- 0.05

betax10.0 <- -0.015
betax1x0.0 <- 0.02
betax1u.0 <- 0.05

betay0.1  <- 4
betayx0.1 <- 0.05
betayx1.1 <- 0.02
betayw.1  <- 0.1
betayu.1  <- 1.5

betay0.0  <- 2
betayx0.0 <- 0.05
betayx1.0 <- 0.02
betayw.0  <- 0.1
betayu.0  <- 1.5

betaz0  <- 0
betazx0 <- 0.1
betazu  <- -1.5
betaza  <- -0.5

betaw0  <- 0
betawx0 <- 0.1
betawu  <- 2

betaa0  <- 0
betaax0 <- 0.02
betaau  <- 0.1

if(Sim.Setup.Obs){
  betad0.0  <- -0.5  # lower -> PrD0|A0 becomes smaller
  betadx0.0 <- 0.01
  betadx1.0 <- 0.02
  betadu.0  <- 0.025
  betadz.0  <- 0
  
  betad0.1  <- -0.5  # lower -> PrD0|A1 becomes smaller
  betadx0.1 <- 0.01
  betadx1.1 <- -0.02
  betadu.1  <- -0.125
  betadz.1  <- 0
} else {
  betad0.0  <- -0.5  # lower -> PrD0|A0 becomes smaller
  betadx0.0 <- 0.01
  betadx1.0 <- 0.02
  betadu.0  <- 0.1
  betadz.0  <- 0.075
  
  betad0.1  <- -0.5  # lower -> PrD0|A1 becomes smaller
  betadx0.1 <- 0.01
  betadx1.1 <- -0.02
  betadu.1  <- -0.2
  betadz.1  <- -0.15
}



deltax10 <- betax10.0 - betax10.1
deltad0  <-  betad0.0  - betad0.1
deltadx0 <-  betadx0.0 - betadx0.1
deltadx1 <-  betadx1.0 - betadx1.1
deltadu  <-  betadu.0  - betadu.1

betaxu  <- 0

PrA1 <- 0.5
PrA0 <- 1-PrA1


#####################################

if(Sim.Setup.Obs){
  
  h11.true <- function(w,x0,x1){
    gamma0  <- (betay0.1 - betayu.1*betaw0/betawu)
    gammax0 <- (betayx0.1-betayu.1*betawx0/betawu)
    gammax1 <- betayx1.1
    gammaw  <- (betayw.1+betayu.1/betawu)
    
    return( 
      gamma0 + 
        gammax0*apply(x0,1,sum) + 
        gammax1*apply(x1,1,sum) + 
        gammaw*w
    )
  }
  
  h01.true <- function(w,x0){
    
    gammaw  <- (betayu.1 + betayw.1*betawu + dimX0*betayx1.1*betax1u.0)/betawu
    deltaw  <- (betadu.0 + dimX0*betadx1.0*betax1u.0)/betawu
    
    gammax0 <- betayx0.1 + betayw.1*betawx0 + betayx1.1*betax1x0.0 - gammaw*betawx0
    deltax0 <- betadx0.0 + betadx1.0*betax1x0.0 - deltaw*betawx0
    
    gamma0  <- betay0.1 + betayw.1*betaw0 + dimX0*betayx1.1*betax10.0 + dimX0*betayx1.1*betadx1.0 - gammaw*deltaw - gammaw*betaw0
    delta0  <- betad0.0 + dimX0*betadx1.0*betax10.0 + 0.5*dimX0*betadx1.0^2 - deltaw*betaw0 - 0.5*deltaw^2
    
    return( 
      (gamma0 + gammax0*apply(x0,1,sum) + gammaw*w)*
        exp(delta0 + deltax0*apply(x0,1,sum) + deltaw*w) ) 
  }
  
  h00.true <- function(w,x0){
    
    gammaw  <- (betayu.0 + betayw.0*betawu + dimX0*betayx1.0*betax1u.0)/betawu
    deltaw  <- (betadu.0 + dimX0*betadx1.0*betax1u.0)/betawu
    
    gammax0 <- betayx0.0 + betayw.0*betawx0 + betayx1.0*betax1x0.0 - gammaw*betawx0
    deltax0 <- betadx0.0 + betadx1.0*betax1x0.0 - deltaw*betawx0
    
    gamma0  <- betay0.0 + betayw.0*betaw0 + dimX0*betayx1.0*betax10.0 + dimX0*betayx1.0*betadx1.0 - gammaw*deltaw - gammaw*betaw0
    delta0  <- betad0.0 + dimX0*betadx1.0*betax10.0 + 0.5*dimX0*betadx1.0^2 - deltaw*betaw0 - 0.5*deltaw^2
    
    return( 
      (gamma0 + gammax0*apply(x0,1,sum) + gammaw*w)*
        exp(delta0 + deltax0*apply(x0,1,sum) + deltaw*w) ) 
  }
  
  h2.true <- function(w,x0){
    deltaw  <- (betadu.0 + dimX0*betadx1.0*betax1u.0)/betawu
    deltax0 <- betadx0.0 + betadx1.0*betax1x0.0 - deltaw*betawx0
    delta0  <- betad0.0 + dimX0*betadx1.0*betax10.0 + 0.5*dimX0*betadx1.0^2 - deltaw*betaw0 - 0.5*deltaw^2
    
    return( 
      exp(delta0 + deltax0*apply(x0,1,sum) + deltaw*w)
    )
  }
  
  q0.true <- function(z,x0){
    gammaz  <- betaau/betazu
    gammax0 <- betaax0 -gammaz*betazx0
    gamma0  <- betaa0 - 0.5*gammaz^2 - gammaz*betaz0
    
    return(
      1 + exp( gamma0 + 
                 gammax0*apply(x0,1,sum) + 
                 gammaz*z )
      
    )
  }
  
  q11.true <- function(z,x0,x1){
    
    eta0 <- 0.5*dimX0*(betax10.1^2-betax10.0^2)
    
    gammaz  <- (deltadu - dimX0*betax1u.0*deltax10 - betaau)/betazu
    gammax1 <- deltadx1 + deltax10
    gammax0 <- deltadx0 - betax1x0.0*deltax10 - betaax0 - gammaz*betazx0
    gamma0  <- eta0 + deltad0  - betaa0  - 0.5*gammaz^2 - gammaz*betaz0 - gammaz*betaza
    
    sgammaz  <- (deltadu - dimX0*betax1u.0*deltax10)/betazu
    sgammax1 <- gammax1
    sgammax0 <- deltadx0 - betax1x0.0*deltax10 - sgammaz*betazx0
    sgamma0  <- eta0 + deltad0  - 0.5*sgammaz^2   - sgammaz*betaz0 - sgammaz*betaza
    
    return(
      exp(gamma0 + 
            gammax1*apply(x1,1,sum) +
            gammax0*apply(x0,1,sum) + 
            gammaz*z) +
        exp(sgamma0 + 
              sgammax1*apply(x1,1,sum) +
              sgammax0*apply(x0,1,sum) + 
              sgammaz*z)
    )
    
  }
  
  
  
}

  