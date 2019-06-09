## ---- echo=FALSE, message=FALSE------------------------------------------
  library(knitr)
  opts_chunk$set(warning=FALSE)

## ---- echo=FALSE, fig.cap="test caption", fig.width=15-------------------
  library(msm)
  par(mfrow=c(1,2))
  y.seq <- seq(0,4,length.out=100)
  plot(y.seq, ptnorm(y.seq, 1,1,lower=0), col="purple", lwd=3, type="l",
       ylab="cdf", xlab="y")
  segments(x0=0, x1=qtnorm(0.5,1,1,lower=0), y0=0.5, lty=2, lwd=3)
  segments(x0=qtnorm(0.5,1,1,lower=0), y0=0, y1=0.5, lty=2, lwd=3)
  segments(x0=0, x1=qtnorm(0.9,1,1,lower=0), y0=0.9, lty=2, lwd=3)
  segments(x0=qtnorm(0.9,1,1,lower=0), y0=0, y1=0.9, lty=2, lwd=3)
  tau <- seq(0,1,.01)
  plot(tau, qtnorm(tau, 1,1,lower=0), col="red", lwd=3, type="l",
       ylab="quantiles", xlab="tau")
  segments(x0=0, x1=0.5, y0=qtnorm(0.5,1,1,lower=0), lty=2, lwd=3)
  segments(x0=0.5, y0=qtnorm(0.5,1,1,lower=0), y1=0, lty=2, lwd=3)
  segments(x0=0, x1=0.9, y0=qtnorm(0.9,1,1,lower=0), lty=2, lwd=3)
  segments(x0=0.9, y0=qtnorm(0.9,1,1,lower=0), y1=0, lty=2, lwd=3)

## ---- echo=FALSE, fig.width=15-------------------------------------------
par(mfrow=c(1,2))
y.seq <- seq(0,6,length.out=100)
plot(y.seq, ptnorm(y.seq, 1,1,lower=0), col="purple", lwd=3, type="l",
     ylab="cdf", xlab="y")
lines(y.seq, ptnorm(y.seq, 2,2,lower=0), col="green", lwd=3, type="l")
segments(x0=qtnorm(0.5,1,1,lower=0), x1=qtnorm(0.5,2,2,lower=0), y0=0.5, lty=2, lwd=3)
segments(x0=qtnorm(0.9,1,1,lower=0), x1=qtnorm(0.9,2,2,lower=0), y0=0.9, lty=2, lwd=3)

plot(tau, qtnorm(tau, 2,2,lower=0) - qtnorm(tau, 1,1,lower=0), col="blue", lwd=3, type="l", xlab="tau", ylab="qte")

## ------------------------------------------------------------------------
library(qte)

## ------------------------------------------------------------------------
data(lalonde)

## ------------------------------------------------------------------------
class(jt.rand)

## ------------------------------------------------------------------------
summary(jt.rand)

## ------------------------------------------------------------------------
ggqte(jt.rand)

## ------------------------------------------------------------------------
 ggqte(jt.cia)

## ------------------------------------------------------------------------
 ggqte(jt.ciat)

## ------------------------------------------------------------------------
 head(lalonde.psid.panel)[,c("id", "year", "re", "treat",
  "age", "education", "black", "hispanic", "married", "nodegree", "u75")]

## ------------------------------------------------------------------------
 ggqte(jt.pqtet)

## ------------------------------------------------------------------------
 plot(b1) 

