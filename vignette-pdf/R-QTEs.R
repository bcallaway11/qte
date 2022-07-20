## ---- echo=FALSE, message=FALSE-----------------------------------------------
  library(knitr)
  ##opts_chunk$set(warning=FALSE)


## ---- echo=FALSE, fig.cap="test caption", fig.width=15------------------------
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


## ---- echo=FALSE, fig.width=15------------------------------------------------
par(mfrow=c(1,2))
y.seq <- seq(0,6,length.out=100)
plot(y.seq, ptnorm(y.seq, 1,1,lower=0), col="purple", lwd=3, type="l",
     ylab="cdf", xlab="y")
lines(y.seq, ptnorm(y.seq, 2,2,lower=0), col="green", lwd=3, type="l")
segments(x0=qtnorm(0.5,1,1,lower=0), x1=qtnorm(0.5,2,2,lower=0), y0=0.5, lty=2, lwd=3)
segments(x0=qtnorm(0.9,1,1,lower=0), x1=qtnorm(0.9,2,2,lower=0), y0=0.9, lty=2, lwd=3)

plot(tau, qtnorm(tau, 2,2,lower=0) - qtnorm(tau, 1,1,lower=0), col="blue", lwd=3, type="l", xlab="tau", ylab="qte")


## -----------------------------------------------------------------------------
library(qte)


## -----------------------------------------------------------------------------
data(lalonde)


## ----ci_qtet, cache=TRUE------------------------------------------------------
jt.rand <- ci.qtet(re78 ~ treat, data=lalonde.exp, probs=seq(0.05,0.95,0.05), 
 se=T, iters=10)


## -----------------------------------------------------------------------------
class(jt.rand)


## -----------------------------------------------------------------------------
summary(jt.rand)


## -----------------------------------------------------------------------------
ggqte(jt.rand)


## ----ci_qte, cache=TRUE-------------------------------------------------------
 jt.cia <- ci.qte(re78 ~ treat,
   xformla=~age + education + black + hispanic + married + nodegree,
   data=lalonde.psid,
   probs=seq(0.05,0.95,0.05), se=T, iters=10)
 summary(jt.cia)


## -----------------------------------------------------------------------------
 ggqte(jt.cia)


## ---- cache=TRUE--------------------------------------------------------------
 jt.ciat <- ci.qtet(re78 ~ treat,
   xformla=~age + education + black + hispanic + married + nodegree,
   data=lalonde.psid,
   probs=seq(0.05,0.95,0.05), se=T, iters=10)
 summary(jt.ciat)


## -----------------------------------------------------------------------------
 ggqte(jt.ciat)


## -----------------------------------------------------------------------------
 head(lalonde.psid.panel)[,c("id", "year", "re", "treat",
  "age", "education", "black", "hispanic", "married", "nodegree", "u75")]


## ----panel_qtet, cache=TRUE, warning=FALSE------------------------------------
 jt.pqtet <- panel.qtet(re ~ treat, t=1978, tmin1=1975, tmin2=1974,
   tname="year", idname="id",
   xformla=~age + education + black + hispanic + married + nodegree,
   data=lalonde.psid.panel,
   probs=seq(0.05,0.95,0.05), se=T, iters=10)
 summary(jt.pqtet)


## -----------------------------------------------------------------------------
 ggqte(jt.pqtet)


## ----bounds, cache=TRUE-------------------------------------------------------
 b1 <- bounds(re ~ treat, t=1978, tmin1=1975, data=lalonde.psid.panel,
   idname="id", tname="year")
 summary(b1)


## -----------------------------------------------------------------------------
 plot(b1) 


## ----ddid2, cache=TRUE, warning=FALSE-----------------------------------------
  jt.ddid2 <- ddid2(re ~ treat, t=1978, tmin1=1975,
    tname="year", idname="id",
    xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
    data=lalonde.psid.panel,
    probs=seq(0.05,0.95,0.05), se=T, iters=10)
  summary(jt.ddid2)
  ggqte(jt.ddid2)


## ----qdid, cache=TRUE, warning=FALSE------------------------------------------
  jt.qdid <- QDiD(re ~ treat, t=1978, tmin1=1975,
     tname="year", idname="id",
     xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
     data=lalonde.psid.panel,
     probs=seq(0.05,0.95,0.05), se=T, iters=10, panel=T)

summary(jt.qdid)
  ggqte(jt.qdid)


## ----cic, cache=TRUE, warning=FALSE-------------------------------------------
  jt.cic <- CiC(re ~ treat, t=1978, tmin1=1975,
    tname="year", idname="id",
    xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
    data=lalonde.psid.panel,
    probs=seq(0.05,0.95,0.05), se=T, iters=10, panel=T)

summary(jt.cic)
ggqte(jt.cic)

