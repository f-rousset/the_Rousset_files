######################################################
##
## Code for Figures 2 and 3
##
######################################################

if ( ! ("RandomFields" %in% installed.packages()) ) {
  install.packages("RandomFields",
                   repos="http://cran.us.r-project.org")
}
library(RandomFields) ## loading Martin Schlather's library for simulation of random field


######################################################
# Simple Mantel test
# X(s) Y(s) zero mean  unit variance 
n <- 50 ## nb of sampling sites
nperm <- 10000 ## nb of permutations
nsim <- 200 ## nb of independent simulations
coord <- matrix(nr=n,nc=2,runif(2*n)) ## spatial coordinates 

SCALE <- c(0,.3,.7) ## scale parameters 
PVAL <- array(dim=c(length(SCALE),nsim,4)) ## arrays to store "p-values" returned by the Mantel tests 

for(iscale in 1:length(SCALE))
    {
      scale.par <- SCALE[iscale]
      for(isim in 1:nsim)
        {
          print(isim)
          if(scale.par==0)
            {
              a <- rnorm(n)
              b <- rnorm(n)}else
          {
            ## simulate data from a Gaussian random field
            ## no trend
            grf <- GaussRF(x=coord[,1],
                           y=coord[,2],
                           grid=FALSE,
                           model="expo",
                           param=c(0,1,0,scale.par), ## unit variance  
                           n=2)
            a <- grf[,1]
            b <- grf[,2]
          }
          A <- as.vector(dist(a))
          B <- as.vector(dist(b))
          rAB <- cor(A,B)
          
          ## permutation
          R <- rep(NA,nperm)
          for(iperm in 1:nperm)
            {
              perm <- sample(n)
              A.perm <- as.vector(dist(a[perm]))
              R[iperm] <- cor(A.perm,B)
            }
          PVAL[iscale,isim,1] <- mean(abs(R)>abs(rAB))
        }
    }



######################################################
# Partial Mantel test
# X(s) Y(s) zero mean  unit variance GRFs
## partial Mantel Dx,Dy | Ds
n <- 50
nperm <- 10000
nsim <- 200
coord <- matrix(nr=n,nc=2,runif(2*n))

SCALE <- c(0,.3,.7)
PVAL <- array(dim=c(length(SCALE),nsim,4))


for(iscale in 1:length(SCALE))
    {
      scale.par <- SCALE[iscale]
      for(isim in 1:nsim)
        {
          print(isim)
          if(scale.par==0)
            {
              a <- rnorm(n)
              b <- rnorm(n)}else
          {
            ## simulate data
            grf <- GaussRF(x=coord[,1],
                           y=coord[,2],
                           grid=FALSE,
                           model="expo",
                           param=c(0,1,0,scale.par),
                           n=2)
            a <- grf[,1]
            b <- grf[,2]
          }
          A <- as.vector(dist(a))
          B <- as.vector(dist(b))
          C <- as.vector(dist(coord))
          rAB <- cor(A,B)
          rAC <- cor(A,C)
          rBC <- cor(B,C)
          rAB.C <- (rAB -  rAC * rBC)/(sqrt(1-rAC^2)*sqrt(1-rBC^2))
          
          ## Permutation: method 1 in Legendre and Fortin, Mol. Ecol.  2010
          A.mat <- matrix(nr=n,nc=n,data=0)
          A.mat[lower.tri(A.mat)] <- A
          A.mat <- A.mat + t(A.mat)
          R <- rep(NA,nperm)
          for(iperm in 1:nperm)
            {
              perm <- sample(n)
              A.mat.perm <- A.mat[perm,perm]
              A.perm <- A.mat.perm[lower.tri(A.mat.perm)]
              r1 <-  cor(A.perm,B)
              r2 <- cor(A.perm,C)
              R[iperm] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
            }
          PVAL[iscale,isim,1] <- mean(abs(R)>abs(rAB.C))
          
          ## ## method 2
          res.A.C <- lm(A~C)$residuals
          res.A.C.mat <- matrix(nr=n,nc=n,data=0)
          res.A.C.mat[lower.tri(res.A.C.mat)] <- res.A.C
          res.A.C.mat <- res.A.C.mat + t(res.A.C.mat)
          R <- rep(NA,nperm)
          for(iperm in 1:nperm)
            {
              perm <- sample(n)
              res.A.C.mat.perm <- res.A.C.mat[perm,perm]
              res.A.C.perm <- res.A.C.mat.perm[lower.tri(res.A.C.mat.perm)]
              r1 <- cor(res.A.C.perm,B)
              r2 <- cor(res.A.C.perm,C)
              R[iperm] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
            }
          PVAL[iscale,isim,2] <- mean(abs(R)>abs(rAB.C))

          ## method 3
          res.A.C <- lm(A~C)$residuals
          res.B.C <- lm(B~C)$residuals
          res.A.C.mat <- matrix(nr=n,nc=n,data=0)
          res.A.C.mat[lower.tri(res.A.C.mat)] <- res.A.C
          res.A.C.mat <- res.A.C.mat + t(res.A.C.mat)
          R <- rep(NA,nperm)
          for(iperm in 1:nperm)
            {
              perm <- sample(n)
              res.A.C.mat.perm <- res.A.C.mat[perm,perm]
              res.A.C.perm <- res.A.C.mat.perm[lower.tri(res.A.C.mat.perm)]
              R[iperm] <- cor(res.A.C.perm, res.B.C)
            }
          PVAL[iscale,isim,3] <- mean(abs(R)>abs(rAB.C))
          
          ## method 4
          res.A.BC <- lm(A~B+C)$residuals
          res.A.BC.mat <- matrix(nr=n,nc=n,data=0)
          res.A.BC.mat[lower.tri(res.A.BC.mat)] <- res.A.BC
          res.A.BC.mat <- res.A.BC.mat + t(res.A.BC.mat)
          R <- rep(NA,nperm)
          for(iperm in 1:nperm)
            {
              perm <- sample(n)
              res.A.BC.mat.perm <- res.A.BC.mat[perm,perm]
              res.A.BC.perm <- res.A.BC.mat.perm[lower.tri(res.A.BC.mat.perm)]
              r1 <- cor(res.A.BC.perm ,B)
              r2 <- cor(res.A.BC.perm ,C)
              R[iperm] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
            }
          PVAL[iscale,isim,4] <- mean(abs(R)>abs(rAB.C))
        }
    }






###################################################################
## Partial Mantel test 
## Model
## X(s)  = linear trend + RF
## Y(s)  = linear trend + RF
## partial Mantel Dx,Dy | Ds



#####
n <- 50
nperm <- 10000 
nsim <- 200 
coord <- matrix(nr=n,nc=2,runif(2*n))

SCALE <- c(0,.3,.7)
PVAL <- array(dim=c(length(SCALE),nsim,4))

for(iscale in 1:length(SCALE))
    {
      scale.par <- SCALE[iscale]
      for(isim in 1:nsim)
        {
          print(isim)
          if(scale.par==0)
            {
              a <- rnorm(n)
              b <- rnorm(n)}else
          {
            ## simulate data
            grf <- GaussRF(x=coord[,1],
                           y=coord[,2],
                           grid=FALSE,
                           model="expo",
                           param=c(0,1,0,scale.par),
                           n=2)
            beta1 <- rnorm(1)
            beta2 <- rnorm(1)
            a <- beta1*coord[,1] + beta2*coord[,2] + grf[,1]
            beta1 <- rnorm(1)
            beta2 <- rnorm(1)
            b <- beta1*coord[,1] + beta2*coord[,2] + grf[,2]
          }
          A <- as.vector(dist(a))
          B <- as.vector(dist(b))
          C <- as.vector(dist(coord))
          rAB <- cor(A,B)
          rAC <- cor(A,C)
          rBC <- cor(B,C)
          rAB.C <- (rAB -  rAC * rBC)/(sqrt(1-rAC^2)*sqrt(1-rBC^2))
          
          ## method 1
          A.mat <- matrix(nr=n,nc=n,data=0)
          A.mat[lower.tri(A.mat)] <- A
          A.mat <- A.mat + t(A.mat)
          R <- rep(NA,nperm)
          for(iperm in 1:nperm)
            {
              perm <- sample(n)
              A.mat.perm <- A.mat[perm,perm]
              A.perm <- A.mat.perm[lower.tri(A.mat.perm)]
              r1 <-  cor(A.perm,B)
              r2 <- cor(A.perm,C)
              R[iperm] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
            }
          PVAL[iscale,isim,1] <- mean(abs(R)>abs(rAB.C))
          
          ## ## method 2
          res.A.C <- lm(A~C)$residuals
          res.A.C.mat <- matrix(nr=n,nc=n,data=0)
          res.A.C.mat[lower.tri(res.A.C.mat)] <- res.A.C
          res.A.C.mat <- res.A.C.mat + t(res.A.C.mat)
          R <- rep(NA,nperm)
          for(iperm in 1:nperm)
            {
              perm <- sample(n)
              res.A.C.mat.perm <- res.A.C.mat[perm,perm]
              res.A.C.perm <- res.A.C.mat.perm[lower.tri(res.A.C.mat.perm)]
              r1 <- cor(res.A.C.perm,B)
              r2 <- cor(res.A.C.perm,C)
              R[iperm] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
            }
          PVAL[iscale,isim,2] <- mean(abs(R)>abs(rAB.C))

          ## method 3
          res.A.C <- lm(A~C)$residuals
          res.B.C <- lm(B~C)$residuals
          res.A.C.mat <- matrix(nr=n,nc=n,data=0)
          res.A.C.mat[lower.tri(res.A.C.mat)] <- res.A.C
          res.A.C.mat <- res.A.C.mat + t(res.A.C.mat)
          R <- rep(NA,nperm)
          for(iperm in 1:nperm)
            {
              perm <- sample(n)
              res.A.C.mat.perm <- res.A.C.mat[perm,perm]
              res.A.C.perm <- res.A.C.mat.perm[lower.tri(res.A.C.mat.perm)]
              R[iperm] <- cor(res.A.C.perm, res.B.C)
            }
          PVAL[iscale,isim,3] <- mean(abs(R)>abs(rAB.C))
          
          ## method 4
          res.A.BC <- lm(A~B+C)$residuals
          res.A.BC.mat <- matrix(nr=n,nc=n,data=0)
          res.A.BC.mat[lower.tri(res.A.BC.mat)] <- res.A.BC
          res.A.BC.mat <- res.A.BC.mat + t(res.A.BC.mat)
          R <- rep(NA,nperm)
          for(iperm in 1:nperm)
            {
              perm <- sample(n)
              res.A.BC.mat.perm <- res.A.BC.mat[perm,perm]
              res.A.BC.perm <- res.A.BC.mat.perm[lower.tri(res.A.BC.mat.perm)]
              r1 <- cor(res.A.BC.perm ,B)
              r2 <- cor(res.A.BC.perm ,C)
              R[iperm] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
            }
          PVAL[iscale,isim,4] <- mean(abs(R)>abs(rAB.C))
        }
    }


####################
## plotting results 
par(mfrow=c(3,3))
par(mai=.4*c(.8,1.5,1.5,.3))
cex <- .8
cex.lab <- 1.
lwd <- .3
for(iscale in 1:length(SCALE))
  {

    xlab <- ifelse(iscale==3,"Quantiles uniform distribution","")

    
    plot((1:nsim)/nsim,sort(PVAL.SIMPLE[iscale,,1]),
         ylab=paste("p-values simulated datasets"),
         xlab=xlab,
         cex.lab=cex.lab,
         main=paste("Scale parameter=",signif(SCALE[iscale],dig=1),"\n",
           "type I error rate=",signif(mean(PVAL.SIMPLE[iscale,,1]<0.05),dig=3)),
         type="n")
    

    points((1:nsim)/nsim,sort(PVAL.SIMPLE[iscale,,1]),col=1,
           pch=1,cex=cex)
    
    abline(0,1,col=1,lty=2)

    plot((1:nsim)/nsim,sort(PVAL.NOTREND[iscale,,1]),
         ylab="",
         xlab=xlab,
         main=paste("Scale parameter=",signif(SCALE[iscale],dig=1),"\n",
           "type I error rate=",signif(mean(PVAL.NOTREND[iscale,,1]<0.05),dig=3)),
         type="n",
           cex.lab=cex.lab)
    
    points((1:nsim)/nsim,sort(PVAL.NOTREND[iscale,,1]),cex=cex,lwd=lwd)
    points((1:nsim)/nsim,sort(PVAL.NOTREND[iscale,,2]),col=2,pch=2,cex=cex,lwd=lwd)
    points((1:nsim)/nsim,sort(PVAL.NOTREND[iscale,,3]),col=3,pch=3,cex=cex,lwd=lwd)
    points((1:nsim)/nsim,sort(PVAL.NOTREND[iscale,,4]),col=4,pch=4,cex=cex,lwd=lwd)
    
    abline(0,1,col=1,lty=2)
    
     plot((1:nsim)/nsim,sort(PVAL.TREND[iscale,,1]),
         ylab="",
         xlab=xlab,
          main=paste("Scale parameter=",signif(SCALE[iscale],dig=1),"\n",
            "type I error rate=",signif(mean(PVAL.TREND[iscale,,1]<0.05),dig=3)),
         type="n" ,
          cex.xlab=cex.lab,)
    
    points((1:nsim)/nsim,sort(PVAL.TREND[iscale,,1]),cex=cex,lwd=lwd)
    points((1:nsim)/nsim,sort(PVAL.TREND[iscale,,2]),col=2,pch=2,cex=cex,lwd=lwd)
    points((1:nsim)/nsim,sort(PVAL.TREND[iscale,,3]),col=3,pch=3,cex=cex,lwd=lwd)
    points((1:nsim)/nsim,sort(PVAL.TREND[iscale,,4]),col=4,pch=4,cex=cex,lwd=lwd)
    
    abline(0,1,col=1,lty=2)

     leg.txt=c("Method 1",
      "Method 2",
      "Method 3",
      "Method 4")
    y.leg <- c(.6,.7,.8,.9)
    pch <- col <- 1:4
    if(iscale==3)
      {
        legend(.005,1,leg.txt,pch=pch,col=1:4,cex=1.1)
      }
  }

######################################################
##
## End code for Figures 2 and 3
##
######################################################


######################################################
##
## Code for Figures 4
##
######################################################


# Partial Mantel test
# X(s) a predictor variable from the original study, 
# Diggle et al., Ann Trop Med. Parasitol, 6, 499-509, 2007
# Y(s) an independent zero mean unit variance GRFs 
#   with correlation function as in original study 
## partial Mantel Dx,Dy | Ds

loa<-read.table("./data/loaloa.txt",sep=" ",skip=1,header=TRUE)
d<-as.matrix(dist(loa[,c("long","lat")]))
corrmat <- exp(-d/0.7) ## 0.7 is the value reported in the original study
cholCorr <- t(chol(corrmat))

n <- nrow(loa) ## i.e., 196
nperm <- 10000
nsim <- 200
coord <- loa[,c("long","lat")] ## geographical coordinates in the same units as the correlation parameter

C<-as.vector(dist(loa[,c("long","lat")]))
SCALE <- NA
PVAL <- array(dim=c(nsim,4))
set.seed(123)

for(isim in 1:nsim)
  {
    print(isim)
    {
      a <- loa$elev ## elevation, the first predictor variable considered in the original study
      ## efficient simulation of multivariate Gaussian with correlation matrix cholCorr %*% t(cholCorr) = corrmat
      b <- cholCorr %*% runif(n) 
    }
    A <- as.vector(dist(a))
    B <- as.vector(dist(b))
    rAB <- cor(A,B)
    rAC <- cor(A,C)
    rBC <- cor(B,C)
    rAB.C <- (rAB -  rAC * rBC)/(sqrt(1-rAC^2)*sqrt(1-rBC^2))
                    
    ## method 1
    A.mat <- matrix(nr=n,nc=n,data=0)
    A.mat[lower.tri(A.mat)] <- A
    A.mat <- A.mat + t(A.mat)
    R <- rep(NA,nperm)
    for(iperm in 1:nperm)
      {
        perm <- sample(n)
        A.mat.perm <- A.mat[perm,perm]
        A.perm <- A.mat.perm[lower.tri(A.mat.perm)]
        r1 <-  cor(A.perm,B)
        r2 <- cor(A.perm,C)
        R[iperm] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
      }
    PVAL[isim,1] <- mean(abs(R)>abs(rAB.C))
    
    ## ## method 2
    res.A.C <- lm(A~C)$residuals
    res.A.C.mat <- matrix(nr=n,nc=n,data=0)
    res.A.C.mat[lower.tri(res.A.C.mat)] <- res.A.C
    res.A.C.mat <- res.A.C.mat + t(res.A.C.mat)
    R <- rep(NA,nperm)
    for(iperm in 1:nperm)
      {
        perm <- sample(n)
        res.A.C.mat.perm <- res.A.C.mat[perm,perm]
        res.A.C.perm <- res.A.C.mat.perm[lower.tri(res.A.C.mat.perm)]
        r1 <- cor(res.A.C.perm,B)
        r2 <- cor(res.A.C.perm,C)
        R[iperm] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
      }
    PVAL[isim,2] <- mean(abs(R)>abs(rAB.C))
    
    ## method 3
    res.A.C <- lm(A~C)$residuals
    res.B.C <- lm(B~C)$residuals
    res.A.C.mat <- matrix(nr=n,nc=n,data=0)
    res.A.C.mat[lower.tri(res.A.C.mat)] <- res.A.C
    res.A.C.mat <- res.A.C.mat + t(res.A.C.mat)
    R <- rep(NA,nperm)
    for(iperm in 1:nperm)
      {
        perm <- sample(n)
        res.A.C.mat.perm <- res.A.C.mat[perm,perm]
        res.A.C.perm <- res.A.C.mat.perm[lower.tri(res.A.C.mat.perm)]
        R[iperm] <- cor(res.A.C.perm, res.B.C)
      }
    PVAL[isim,3] <- mean(abs(R)>abs(rAB.C))
    
    ## method 4
    res.A.BC <- lm(A~B+C)$residuals
    res.A.BC.mat <- matrix(nr=n,nc=n,data=0)
    res.A.BC.mat[lower.tri(res.A.BC.mat)] <- res.A.BC
    res.A.BC.mat <- res.A.BC.mat + t(res.A.BC.mat)
    R <- rep(NA,nperm)
    for(iperm in 1:nperm)
      {
        perm <- sample(n)
        res.A.BC.mat.perm <- res.A.BC.mat[perm,perm]
        res.A.BC.perm <- res.A.BC.mat.perm[lower.tri(res.A.BC.mat.perm)]
        r1 <- cor(res.A.BC.perm ,B)
        r2 <- cor(res.A.BC.perm ,C)
        R[iperm] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
      }
    PVAL[isim,4] <- mean(abs(R)>abs(rAB.C))
  }
}


cex <- 2
lwd <- .2
# pdf("qqplot_loaloa.pdf")
plot((1:nsim)/nsim,sort(PVAL[,1]),
xlab=paste("P values for simulated datasets"),
ylab="Quantiles of a uniform distribution",
type="n")

points((1:nsim)/nsim,sort(PVAL[,1]),cex=cex,lwd=lwd)
points((1:nsim)/nsim,sort(PVAL[,2]),col=2,pch=2,cex=cex,lwd=lwd)
points((1:nsim)/nsim,sort(PVAL[,3]),col=3,pch=3,cex=cex,lwd=lwd)
points((1:nsim)/nsim,sort(PVAL[,4]),col=4,pch=4,cex=cex,lwd=lwd)

abline(0,1,col=1,lty=2)

leg.txt=c("Method 1",
"Method 2",
"Method 3",
"Method 4")
y.leg <- c(.6,.7,.8,.9)
pch <- col <- 1:4
{
  legend(.05,1,leg.txt,pch=2:5,col=1:4,cex=1.1)
}


#dev.off()

######################################################
##
## End code for Figures 4
##
######################################################



######################################################
##
## Code for Figures 5
##
######################################################

if ( ! ("RandomFields" %in% installed.packages()) ) {
  install.packages("RandomFields",
                   repos="http://cran.us.r-project.org")
}
library(RandomFields) ## loading Martin Schlather's library for simulation of random field

######################################################
# Simple Mantel test
# X(s) Y(s) zero mean  unit variance 
n <- 200 ## nb of sampling sites
nsim <- 10000 ## nb of independent simulations
coord <- matrix(nr=n,nc=2,runif(2*n)) ## spatial coordinates 

SCALE <- c(0,.3,.7) ## scale parameters 
RHO.SIMPLE <- array(dim=c(length(SCALE),nsim,8)) ## arrays to store correlation coeff. (true one and the one returned by the Mantel tests )

for(iscale in 1:length(SCALE))
    {
      scale.par <- SCALE[iscale]
      for(isim in 1:nsim)
        {
          print(isim)
          if(scale.par==0)
            {
              a <- rnorm(n)
              b <- rnorm(n)}else
          {
            ## simulate data from a Gaussian random field
            ## no trend
            grf <- GaussRF(x=coord[,1],
                           y=coord[,2],
                           grid=FALSE,
                           model="expo",
                           param=c(0,1,0,scale.par), ## unit variance  
                           n=2)
            a <- grf[,1]
            b <- grf[,2]
          }
          A <- as.vector(dist(a))
          B <- as.vector(dist(b))
          RHO.SIMPLE[iscale,isim,1] <- cor(A,B)
          
          ## permutation
          perm <- sample(n)
          A.perm <- as.vector(dist(a[perm]))
          RHO.SIMPLE[iscale,isim,2] <- cor(A.perm,B)
        }
    }



######################################################
# Partial Mantel test
# X(s) Y(s) zero mean  unit variance GRFs
## partial Mantel Dx,Dy | Ds
RHO.NOTREND <- array(dim=c(length(SCALE),nsim,4))


for(iscale in 1:length(SCALE))
    {
      scale.par <- SCALE[iscale]
      for(isim in 1:nsim)
        {
          print(isim)
          if(scale.par==0)
            {
              a <- rnorm(n)
              b <- rnorm(n)}else
          {
            ## simulate data
            grf <- GaussRF(x=coord[,1],
                           y=coord[,2],
                           grid=FALSE,
                           model="expo",
                           param=c(0,1,0,scale.par),
                           n=2)
            a <- grf[,1]
            b <- grf[,2]
          }
          A <- as.vector(dist(a))
          B <- as.vector(dist(b))
          C <- as.vector(dist(coord))
          rAB <- cor(A,B)
          rAC <- cor(A,C)
          rBC <- cor(B,C)
          RHO.NOTREND[iscale,isim,1] <-  rAB.C <- (rAB -  rAC * rBC)/(sqrt(1-rAC^2)*sqrt(1-rBC^2))
          
          ## Permutation: method 1 in Legendre and Fortin, Mol. Ecol.  2010
          A.mat <- matrix(nr=n,nc=n,data=0)
          A.mat[lower.tri(A.mat)] <- A
          A.mat <- A.mat + t(A.mat)
          perm <- sample(n)
          A.mat.perm <- A.mat[perm,perm]
          A.perm <- A.mat.perm[lower.tri(A.mat.perm)]
          r1 <-  cor(A.perm,B)
          r2 <- cor(A.perm,C)
          RHO.NOTREND[iscale,isim,2] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
        }
    }




###################################################################
## Partial Mantel test 
## Model
## X(s)  = linear trend + RF
## Y(s)  = linear trend + RF
## partial Mantel Dx,Dy | Ds
RHO.TREND <- array(dim=c(length(SCALE),nsim,4))

for(iscale in 1:length(SCALE))
    {
      scale.par <- SCALE[iscale]
      for(isim in 1:nsim)
        {
          print(isim)
          if(scale.par==0)
            {
              a <- rnorm(n)
              b <- rnorm(n)}else
          {
            ## simulate data
            grf <- GaussRF(x=coord[,1],
                           y=coord[,2],
                           grid=FALSE,
                           model="expo",
                           param=c(0,1,0,scale.par),
                           n=2)
            beta1 <- rnorm(1)
            beta2 <- rnorm(1)
            a <- beta1*coord[,1] + beta2*coord[,2] + grf[,1]
            beta1 <- rnorm(1)
            beta2 <- rnorm(1)
            b <- beta1*coord[,1] + beta2*coord[,2] + grf[,2]
          }
          A <- as.vector(dist(a))
          B <- as.vector(dist(b))
          C <- as.vector(dist(coord))
          rAB <- cor(A,B)
          rAC <- cor(A,C)
          rBC <- cor(B,C)
          RHO.TREND[iscale,isim,1] <- rAB.C <- (rAB -  rAC * rBC)/(sqrt(1-rAC^2)*sqrt(1-rBC^2))
          
          ## method 1
          A.mat <- matrix(nr=n,nc=n,data=0)
          A.mat[lower.tri(A.mat)] <- A
          A.mat <- A.mat + t(A.mat)
          perm <- sample(n)
          A.mat.perm <- A.mat[perm,perm]
          A.perm <- A.mat.perm[lower.tri(A.mat.perm)]
          r1 <-  cor(A.perm,B)
          r2 <- cor(A.perm,C)
          RHO.TREND[iscale,isim,2] <- (r1-r2*rBC )/(sqrt(1-r2^2)*sqrt(1-rBC^2))
        }
    }



####################
## plotting results
## pdf("/media/SSD_Gilles/biogillesg/com/articles/2012/Mantel/10/fig/pdf/histo.pdf")
par(mfrow=c(3,3))
par(mai=.4*c(1.3,1.5,1.5,.3))
cex <- .8
cex.lab <- 1.4
lwd <- 1.5
ymax <- 13
for(iscale in 1:length(SCALE))
  {

    main <- paste("Scale parameter=",signif(SCALE[iscale],dig=1))
    xlab <- "Correlation coefficient"#
    plot(density(RHO.SIMPLE[iscale,,1]),lwd=lwd,ylim=c(0,ymax),type="n",
         xlab=xlab,main=main,sub="")
    lines(density(RHO.SIMPLE[iscale,,1]),lwd=lwd,col="cyan4")
    lines(density(RHO.SIMPLE[iscale,,2]),lwd=lwd,lty=2,col="orange")

    plot(density(RHO.NOTREND[iscale,,1]),lwd=lwd,ylim=c(0,ymax),
         xlab=xlab,type="n",main=main,sub="",)
    lines(density(RHO.NOTREND[iscale,,1]),lwd=lwd,col="cyan4")
    lines(density(RHO.NOTREND[iscale,,2]),lwd=lwd,lty=2,col="orange")

    plot(density(RHO.TREND[iscale,,1]),lwd=lwd,ylim=c(0,ymax),
         xlab=xlab,main=main,sub="")
    lines(density(RHO.TREND[iscale,,1]),lwd=lwd,col="cyan4")
    lines(density(RHO.TREND[iscale,,2]),lwd=lwd,lty=2,col="orange")
  }
legend(c("True","Mantel"),
       color=c("cyan4","orange"),lty=c(1,2),lwd=lwd,pos="topleft")

dev.off()



######################################################
##
## End code for Figures 5
##
######################################################
