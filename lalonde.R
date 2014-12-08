require(causalsens)
source("qte/R/qte.r")

data(lalonde.exp)
data(lalonde.psid)


#calculate the average treatment effect on the treated without
#covariates
att = mean(subset(lalonde.exp,treat==1)$re78) - mean(subset(lalonde.exp,treat==0)$re78)

#parameters used below:
#h is the step size passed to various function
h = 0.1
#probs are the values that we compute QTET for
probs = seq(h,1-h,0.02)

#experimental QTE
actual.qte = quantile(subset(lalonde.exp,treat==1)$re78,probs=probs) - quantile(subset(lalonde.exp,treat==0)$re78,probs=probs)
actual.qte.employed = quantile(subset(exp.employed.subset,treat==1)$re78,probs=probs) - quantile(subset(exp.employed.subset,treat==0)$re78, probs=probs)
#plot experimental QTE
plot(probs,actual.qte,type="l")

#set up data for call to fan-yu
#here we make dataset into a panel for calls to functions below
lalonde.psid$id = as.integer(rownames(lalonde.psid))
tempdata3 = cbind(year=1978, lalonde.psid[,c("id","re78","treat","age","education","black",
                                "hispanic","married","nodegree","u74","u75")])
tempdata2 = cbind(year=1975, lalonde.psid[,c("id","re75","treat","age","education","black",
                                "hispanic","married","nodegree","u74","u75")])
tempdata1 = cbind(year=1974, lalonde.psid[,c("id","re74","treat","age","education","black",
                                "hispanic","married","nodegree","u74","u75")])
colnames(tempdata3) = colnames(tempdata2) = colnames(tempdata1) = c("year","id","re","treat","age","education","black",
                                                                    "hispanic","married","nodegree","u74","u75")
lalonde.data = rbind(tempdata3,tempdata2,tempdata1)
lalonde.data$uniqueid = paste(lalonde.data$id,lalonde.data$year,sep="-")

employed.subset = subset(lalonde.data, !(lalonde.data$id %in% lalonde.data[lalonde.data[,"u75"]==1,"id"]))

exp.employed.subset = subset(lalonde.exp, u75==0)

#call firpo (for cross-sectional case) method
lalonde.firpo = firpo(re78 ~ treat, x=c("age","education","black","hispanic",
                      "married","nodegree","u74","u75"),
    data=subset(lalonde.psid, u75==0),
    probs=probs)

#call new panel qte method
lalonde.panel <- panel.qte1(re ~ treat,
                           tname="year",t=1978, tmin1=1975, tmin2=1974,
                           data=employed.subset, idname="id", uniqueid="uniqueid",
                           #x=c("age","education","black","hispanic",
                           #     "married","nodegree","u74","u75"),
                           x=c("age","education"),
                           #x=NULL,
                           #y.seq=seq(0,120000,length.out=20),
                           #dy.seq=seq(-70000,120000,length.out=20), 
                           y.seq=seq(min(lalonde.exp$re78), max(lalonde.exp$re78), length.out=300),#sort(unique(lalonde.exp$re78)),
                           dy.seq=seq(min(lalonde.exp$re78 - lalonde.exp$re75), max(lalonde.exp$re78 - lalonde.exp$re78), length.out=300),#sort(unique(lalonde.exp$re78-lalonde.exp$re75)),
                           probs=probs,
                           dropalwaystreated=FALSE,
                           h=0.37, probevals=500)

#ptm = proc.time()
#Rprof()
#call fan-yu for bounds
lalonde.fy = fan.yu(re ~ treat,
               tname="year",t=1978, tmin1=1975, data=employed.subset,
               idname="id", uniqueid="uniqueid", 
               #x=c("sq"),
               x=NULL,
               y.seq=sort(unique(lalonde.exp$re78)),
               dy.seq=sort(unique(lalonde.exp$re78-lalonde.exp$re75)),
               dropalwaystreated=FALSE,
               probs=probs)
#Rprof(NULL)
#summaryRprof()
#proc.time()-ptm

#ptm = proc.time()
#Rprof()
#call panelDid with 3 periods
lalonde.fy3 = panel.qte(re ~ treat,
                           tname="year",t=1978, tmin1=1975, tmin2=1974,
                           data=employed.subset, idname="id", uniqueid="uniqueid",
                           #x=c("age","education","black","hispanic",
                           #     "married","nodegree","u74","u75"),
                           x=NULL,
                           #y.seq=seq(0,120000,length.out=20),
                           #dy.seq=seq(-70000,120000,length.out=20), 
                           y.seq=seq(min(lalonde.exp$re78), max(lalonde.exp$re78), length.out=300),#sort(unique(lalonde.exp$re78)),
                           dy.seq=seq(min(lalonde.exp$re78 - lalonde.exp$re75), max(lalonde.exp$re78 - lalonde.exp$re78), length.out=300),#sort(unique(lalonde.exp$re78-lalonde.exp$re75)),
                           probs=probs,
                           dropalwaystreated=FALSE,
                           h=0.37, probevals=500)
                           #copula.test=actual.copula)
                           #F.untreated.change.test=actual.F.untreated.change)
                           #F.treated.tmin1.test=actual.F.untreated.initial)
#Rprof(NULL)
#summaryRprof()
#proc.time()-ptm

#ptm = proc.time()
#Rprof()
#call panelDid with covariates
#bw = panel.qte.bw(hvec=c(0.33,0.35,0.37),
lalonde.fy3.cov = panel.qte(re ~ treat, xformla=~age + education + black + 
    hispanic + married + nodegree + u74,
    I(u74*age) + I(u74*education) + I(u74*black) +
    I(u74*hispanic) + I(u74*married) + I(u74*nodegree),
                                tname="year",t=1978, tmin1=1975, tmin2=1974,
                                data=employed.subset, idname="id", uniqueid="uniqueid",
                                #x=c("age","education","black","hispanic",
                                #"married","nodegree","u74"),
                                y.seq=seq(min(lalonde.exp$re78), max(lalonde.exp$re78), length.out=300),#sort(unique(lalonde.exp$re78)),
                           dy.seq=seq(min(lalonde.exp$re78 - lalonde.exp$re75), max(lalonde.exp$re78 - lalonde.exp$re78), length.out=300),#sort(unique(lalonde.exp$re78-lalonde.exp$re75)),
                                probs=probs,
                                dropalwaystreated=FALSE,
                                h=0.37,
    probevals=500)
                                #copula.test=function(u,v) return(u*v))
                                #F.untreated.change.test=actual.F.untreated.change,
                                #F.treated.tmin1.test=actual.F.untreated.initial)
                                
#Rprof(NULL)
#summaryRprof()
#proc.time()-ptm


#ptm = proc.time()
#Rprof()
#call athey-imbens with two periods
lalonde.ai = athey.imbens(re ~ treat, t=1978, tmin1=1975, tname="year", x=NULL,
                     data=lalonde.data, dropalwaystreated=FALSE,
                     idname="id", uniqueid="uniqueid",
                     y.seq=sort(unique(lalonde.exp$re78)),
                     probs=probs)
#Rprof(NULL)
#summaryRprof()
#proc.time()-ptm


###FIGURES###

#1) Experimental Plots
#1.a) Experimental CDF
layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
plot(ecdf(subset(lalonde.exp,treat==1)$re78), main="Experimental Wage Distributions")
lines(ecdf(subset(lalonde.exp,treat==0)$re78),col="blue")
legend("bottomright", c("Treated Wages","Untreated Wages"), 
       col=c("black","blue"), lty=c(1,1),
       lwd=c(3,3))


#1.b) Experimental PDF
plot(density(subset(lalonde.exp,treat==1)$re78),lwd=3,main="Experimental Wage PDFs",
     ylim=c(0,1.5e-04))
lines(density(subset(lalonde.exp,treat==0)$re78),lwd=3,col="blue")
legend("topright", c("Treated Wages","Untreated Wages"), 
       col=c("black","blue"), lty=c(1,1),
       lwd=c(3,3))

#1.c) Experimental QTE
plot(probs,actual.qte,type="l", lwd=3, main="Experimental QTET",
                 xlab="quantile", ylab="QTET")
abline(a=0,b=0,lty=2)
abline(a=att,b=0,lty=2,col="blue")
legend("topleft", c("Experimental QTET", "ATT", "0"), 
       col=c("black","black","blue"), lty=c(1,2,2),
       lwd=c(3,1,1))



#####
#2) No Covariates case
#2.a) Is D-i-D Assumption valid: Go back to first two periods
#for means, these two should be equal
mean.change.treated.tmin1 = mean(subset(lalonde.exp,treat==1)$re75) - mean(subset(lalonde.exp,treat==1)$re74)
mean.change.untreated.tmin1 = mean(subset(lalonde.exp,treat==0)$re75) - mean(subset(lalonde.exp,treat==0)$re74)
png("~/Documents/school/projects/Common App/paper/figures/change-untreated-outcomes-no-cov.png")
plot(density(subset(lalonde.exp,treat==1)$re75 - subset(lalonde.exp,treat==1)$re74), lwd=3, main="Change in Untreated Outcomes", xlim=c(-10000,10000))
lines(density(subset(lalonde.exp,treat==0)$re75 - subset(lalonde.exp,treat==0)$re74), col="blue", lwd=3)
legend("topleft", c("Change Untreated Outcomes Treated Group 1975-1974", "Change Untreated Outcomes Untreated Group 1975-1974"),
       col=c("black","blue"), lty=c(1,1), lwd=c(3,3))
dev.off()

#KS test can test if the two distributions are the same
#ks.test(log(d05.treated.04$applcn),log(d05.untreated.04$applcn))
ks = ks.test((subset(lalonde.exp,treat==1)$re75 - subset(lalonde.exp,treat==1)$re74),
        (subset(lalonde.exp,treat==0)$re75 - subset(lalonde.exp,treat==0)$re74))


#2.b) QTETs
par(mfrow=c(1,1)) #reset plot layout
#png("~/Documents/school/projects/Common App/paper/figures/qtet-estimates-3.png")
plot(probs,actual.qte.employed,type="l", ylim=c(-25000,12000), lwd=3, main="Estimated QTETs",
     xlab="quantile", ylab="QTE")
lines(probs,lalonde.fy3$qte, col="blue", lwd=3)
#lines(probs,lalonde.ai$qte, col="red", lwd=3)
#randomization
lines(probs,
      quantile(subset(lalonde.psid,treat==1)$re78,probs=probs) - 
        quantile(subset(lalonde.psid,treat==0)$re78,probs=probs),
      col="green", lwd=3)
lines(probs, lalonde.firpo$qte, col="purple", lwd=3)
lines(probs, lalonde.fy$ub.qte, lty=2, lwd=3)
lines(probs, lalonde.fy$lb.qte, lty=2, lwd=3)
lines(probs, lalonde.fy3.cov$qte, lwd=3, col="orange")
#thuysbaert
#...
legend("bottomleft", c("Experimental QTE","3 Per.","AI","Random Treatment", "Firpo", "Fan-Yu Bounds", "", "3 Per. w Covariates"), 
       col=c("black","blue","red","green","purple","black","black","orange"), lty=c(1,1,1,1,1,2,2,1),
       lwd=c(3,3,3,3,3,3,3,3))
#dev.off()


#3) covariates case
#3.a) first plot the counterfactual distributions from the covariates case, the no covariates case, and the experimental case
plot(lalonde.fy3$f.treatedcf.t, col="blue", main="counterfactual untreated distribution of wages")
lines(lalonde.fy3.cov$f.treatedcf.t, col="green")
lines(ecdf(subset(lalonde.exp,treat==0)$re78),col="black")
legend("bottomright", c("experimental","3 per. no covariates","3 per. covariates"),
       col=c("black","green","blue"), lty=c(1,1,1),
       lwd=c(3,3,3))


###plot each of the distributions
uvec = seq(min(lalonde.fy3$untreated.change),max(lalonde.fy3$untreated.change),
    length.out=100)
vvec = seq(min(lalonde.fy3$treated.tmin1), max(lalonde.fy3$treated.tmin1),
    length.out=100)
f.joint.t.plotdata = expand.grid(uvec,vvec) #gets every combination
colnames(f.joint.t.plotdata) = c("u","v")
f.joint.t.plotdata$f = lalonde.fy3$f.joint.t(f.joint.t.plotdata$u,
    f.joint.t.plotdata$v)


require(lattice)
ppi <- 1000
#png("~/documents/school/projects/common\ app/paper/figures/figure-%d.png",
#    width=12*ppi, height=12*ppi, res=ppi)
wireframe(f ~ u * v, data=f.joint.t.plotdata,
          xlab="change in outcomes",
          ylab="initial outcome",
          drape=t,
          colorkey=t,
          scales=list(arrows=f))

#levelplot(f ~ u * v, data=f.joint.t.plotdata,
#          col.regions=terrain.colors(100))
dev.off()


###plot the copula functions:
uvec <- seq(0,1,0.05)
vvec <- seq(0,1,0.05)
copula.plotdata <- expand.grid(uvec,vvec)
colnames(copula.plotdata) <- c("u","v")
copula.plotdata$emp.cop <- lalonde.fy3$copula(copula.plotdata[,"u"],
                                              copula.plotdata[,"v"])
copula.plotdata$ind.cop <- copula.plotdata$u * copula.plotdata$v
copula.plotdata$ppd.cop <- apply(copula.plotdata[,c("u","v")], margin=1, fun=min)
copula.plotdata$pnd.cop <- apply(copula.plotdata[,c("u","v")], margin=1,
                                 fun=function(x) max(sum(x) -1, 0))
                                                 

copplot1 <- wireframe(emp.cop ~ u * v, data=copula.plotdata,
          xlab="change in outcomes",
          ylab="initial outcome",
          drape=t,
          colorkey=t,
          col.regions = colorramppalette(c("yellow", "black"))(100),
          scales=list(arrows=f))

copplot2 <- wireframe(ind.cop ~ u * v, data=copula.plotdata,
          xlab="change in outcomes",
          ylab="initial outcome",
          drape=t,
          colorkey=t,
          scales=list(arrows=f))

copplot3 <- wireframe(ppd.cop ~ u * v, data=copula.plotdata,
          xlab="change in outcomes",
          ylab="initial outcome",
          drape=t,
          colorkey=t,
          scales=list(arrows=f))

copplot4 <- wireframe(pnd.cop ~ u * v, data=copula.plotdata,
          xlab="change in outcomes",
          ylab="initial outcome",
          drape=t,
          colorkey=t,
          scales=list(arrows=f))

#png("~/documents/school/projects/common app/paper/figures/figure-copulas.png")
print(copplot1     , split=c(1,1,2,2) , more=true )
print(copplot2    , split=c(2,1,2,2) , more=true )
print(copplot3  , split=c(1,2,2,2) , more=true)
print(copplot4 , split=c(2,2,2,2) )
#dev.off()

#sum(1*(copula.plotdata$emp.cop>copula.plotdata$ppd.cop))
#note that somehow we are getting values for the copula that are above (and below)

#find these problems and get rid of them
tol <- 0.05
problems <- subset(copula.plotdata,
                   (copula.plotdata$emp.cop-tol)>copula.plotdata$ppd.cop |
                   (copula.plotdata$emp.cop+tol)<copula.plotdata$pnd.cop)

#the copula upper (and lower) bound.  this means that there must be some mistake in the calculation of the copula! (not just in the numerical procedure to make
#calculations from it)!.


#plot histogram of propensity score
pscore <- predict(glm(treat ~ age + education + hispanic + married + nodegree +
              u74 + u75, data=lalonde.psid, family=binomial(link="probit")),
                  type="response")
require(ggplot2)
png("~/documents/school/projects/common app/paper/figures/pscore-hist.png")
ggplot(data.frame(pscore=pscore[186:length(pscore)]), aes(x=pscore)) +
    geom_histogram()
dev.off()

sum(1*(pscore[186:length(pscore)]>0.5)) #only 20 observations with pscore>0.5
#in the non-experimental dataset.

#plot the untreated change distributions with and without covariates (these should be different)
png("~/documents/school/projects/common app/paper/figures/change-dist.png")
plot(lalonde.fy3.cov$f.untreated.change, xlim=c(-10000,10000), lwd=3,
     main="change in  outcome between t-1 and t")
lines(lalonde.fy3$f.untreated.change, lwd=3, col="blue")
legend(x="bottomright", legend=c("propensity score reweighted", "unweighted sample"), col=c("black","blue"), lwd=c(3,3))
dev.off()
#lines(lalonde.fy3$f.treated.change.tmin1, col="black")



pscore.dw <- predict(glm(treat ~ age + i(age^2) + education +
                         i(education^2) + nodegree + married + black +
                         hispanic + re74 + i(re74^2) + re75 + i(re75^2) +
                         i(1*(re74==0)) + i(1*(re75==0)) +
                         i(1*(re74==0)*hispanic),
                         data=lalonde.psid,
                         family=binomial(link="probit")), type="response")
require(ggplot2)
#png("~/documents/school/projects/common app/paper/figures/pscore-hist2.png")
ggplot(data.frame(pscore.dw=pscore.dw[186:length(pscore.dw)]),
       aes(x=pscore.dw)) +
    geom_histogram()
#dev.off()


##plot the distribution of outcomes for the treated group in the previous period
#this plot does not depend on covariates whether or not you include them
plot(lalonde.fy3$f.treated.tmin1)



##########plot the estimated outcomes to compare with the actual distn's#####
##plot the estimated copula
print(copplot1, split=c(1,1,2,1), more=t)

#estimate the actual copula
actual.change <- subset(exp.employed.subset,treat==0)$re78 -
    subset(exp.employed.subset,treat==0)$re75
actual.initial <- subset(exp.employed.subset, treat==0)$re75

actual.copula <- getCopula(actual.change, actual.initial)

#this should return a function that takes in u,v and returns
#the value of `thecopula' for that u and v
getCopula <- function(data1, data2, u.seq=seq(0,1,0.02), v.seq=seq(0,1,0.02)) {
    copula.bw <<- npudistbw(~data1 + data2)
    thecopula <<- npcopula(bws=copula.bw,
                          data=data.frame(data1=data1, data2=data2),
                          u=cbind(u.seq,v.seq))
    ###retfun <- function(u,v) {
    #    do.call(what=copulafun, args=u,v,the.npcopula=thecopula)
    #}
    #return(retfun)
    outfun <- function(u,v) copulafun(u,v,the.npcopula=thecopula)
    class(outfun) <- c("copfun", class(outfun))
    outfun
}

  #copula should either
  #(i) take in scalar u and vector v and return vector of length(v)
  #(ii) take in vector u and vector v (of same length) and return
  # vector of length(v)
  copulafun <- function(u,v,the.npcopula) {

      if ( (length(u) != 1) & (length(u) != length(v)) ) {
          stop("u must either be scalar or same length as v")
      }

      #when u is scalar will call this function with vapply
      #call with x=v and z=u
      copula.inner <- function(x,z) {
          which.min((z-the.npcopula$u1)^2 + (x-the.npcopula$u2)^2)
      }

      #when u is vector will call this function with apply
      #call with x=matrix(u,v)
      copula.multiple <- function(x) {
          which.min((x[1]-the.npcopula$u1)^2 + (x[2]-the.npcopula$u2)^2)
      }

      if (length(u) == 1) {
          
          whichys <- vapply(v, FUN=copula.inner, FUN.VALUE=1, z=u)
          
      } else {

          whichys <- apply(cbind(u,v), FUN=copula.multiple, MARGIN=1)

      }
      
       #find the closest point in copula.fun and take its value
      #whichys <- vapply(v,
      #                  fun=function(x) which.min((u-copula.fun$u1)^2 +
      #                      (x-copula.fun$u2)^2),
      #                  fun.value=1)
      #whichval <- (u-copula.fun$u1)^2 + (v-copula.fun$u2)^2
      #whichy <- which.min(whichval)

      actual.copula.fun$copula[whichys]
  }


##Plot copula function
plot.copfun <- function(x, ..., uvec=seq(0,1,0.05), vvec=seq(0,1,0.05)) {
    require(lattice)
    copfun <- x
    plotdata <- expand.grid(uvec,vvec)
    colnames(plotdata) <- c("u","v")
    plotdata$cop <- copfun(plotdata[,"u"], plotdata[,"v"])
    copplot <- wireframe(cop ~ u * v, data=plotdata,...)
          #xlab="Change in outcomes",
          #ylab="Initial Outcome",
          #drape=T,
          #colorkey=T,
          #scales=list(arrows=F))
    copplot
    return(copplot)
}

print.copfun <- function(x,...) {
    cat(x(0.5,0.5))
}

##Estimate the copula for the previous period
actual.change.tmin1 <- subset(exp.employed.subset,treat==0)$re75 -
    subset(exp.employed.subset,treat==0)$re74
actual.initial.tmin1 <- subset(exp.employed.subset, treat==0)$re74

actual.copula.tmin1 <- getCopula(actual.change.tmin1, actual.initial.tmin1)
plot(actual.copula.tmin1)

png("~/Documents/school/projects/Common App/paper/figures/copula-stability-t.png")
plot(actual.copula)
dev.off()

png("~/Documents/school/projects/Common App/paper/figures/copula-stability-tmin1.png")
plot(actual.copula.tmin1)
dev.off()

tcor <- cor(actual.change, actual.initial, method="spearman")
tmin1cor <- cor(actual.change.tmin1, actual.initial.tmin1, method="spearman")


#STILL NEED TO WRITE THIS AS A FUNCTION SO THAT WE CAN CALL IT
####**************FIGURING OUT THE ACTUAL COPULA IS WHERE TO START TOMORROW****##
uvec <- seq(0,1,0.05)
vvec <- seq(0,1,0.05)
actual.copula.plotdata <- expand.grid(uvec,vvec)
actual.colnames(copula.plotdata) <- c("u","v")
actual.copula.plotdata$cop <- actual.copula(actual.copula.plotdata[,"u"],
                                              actual.copula.plotdata[,"v"])
copplot.actual <- wireframe(cop ~ u * v, data=actual.copula.plotdata,
          xlab="Change in outcomes",
          ylab="Initial Outcome",
          drape=T,
          colorkey=T,
          scales=list(arrows=F))

#png("~/Documents/school/projects/Common App/paper/figures/actual-and-estimated-copulas.png")
print(copplot1, split=c(1,1,2,1), more=T, res=1000)
print(copplot.actual, split=c(2,1,2,1))
#dev.off()

png("~/Documents/school/projects/Common App/paper/figures/estimated-copula.png")
copplot1
dev.off()

png("~/Documents/school/projects/Common App/paper/figures/actual-copula.png")
#copplot.actual
plot(actual.copula)
dev.off()

class(lalonde.fy3.cov$copula) <- c("copfun","function")
png("~/Documents/school/projects/Common App/paper/figures/estimated-copula.png")
plot(lalonde.fy3.cov$copula)
dev.off()

#************UNTREATED CHANGE***********#
#png("~/Documents/school/projects/Common App/paper/figures/change-dist.png")
actual.F.untreated.change <- ecdf(actual.change)
plot(actual.F.untreated.change, main="Distribution of change in outcomes",
     ylim=c(0,1.3), xlim=c(-20000,20000))
lines(lalonde.fy3.cov$F.untreated.change, col="blue", lwd=3)
lines(lalonde.fy3$F.untreated.change, col="red", lwd=3)
legend("bottomright", c("Actual", "P-score Reweighted", "Unweighted"),
       col=c("black","blue","red"), lwd=3)
#dev.off()

#untreated initial
png("~/Documents/school/projects/Common App/paper/figures/initial-dist.png")
actual.F.untreated.initial <- ecdf(actual.initial)
plot(actual.F.untreated.initial,
     main="Distribution of initial untreated potential outcomes")
lines(lalonde.fy3.cov$F.treated.tmin1,col="blue") #they are untreated at t-1
#lines(lalonde.fy3$F.treated.tmin1,col="red") #don't need to plot bc they are exactly the same
legend("bottomright", c("Future Experimental Untreated t-1" ,
                        "Future Treated t-1"),
       col=c("black","blue"), lwd=c(3,3))
dev.off()


##make the QTET plot when using the actual copulas
par(mfrow=c(1,1)) #reset plot layout
png("~/Documents/school/projects/Common App/paper/figures/copula-test.png")
plot(probs,actual.qte.employed,type="l", ylim=c(-8000,8000), lwd=3, main="Estimated QTETs",
     xlab="quantile", ylab="QTE")
lines(probs,lalonde.fy3$qte, col="blue", lwd=3)
#lines(probs,lalonde.fy3.copulatest$qte, col="blue", lty=2, lwd=3)
lines(probs,lalonde.fy3.cov$qte, col="orange", lwd=3)
lines(probs,lalonde.fy3.cov.copulatest$qte, col="green", lty=2, lwd=3)
lines(probs,lalonde.fy3.cov.indcopulatest$qte, col="red", lty=2, lwd=3)
#lines(probs, lalonde.fy$ub.qte, lty=2, lwd=3)
#lines(probs, lalonde.fy$lb.qte, lty=2, lwd=3)
#thuysbaert
#...
#legend("bottomleft", c("Experimental QTE","3 Per. Est. Copula","3 Per. Actual Copula","3 Per. w Covariates & Est. Copula", "3 Per. w Covs & Actual Copula", "Fan-Yu Bounds", ""), 
#       col=c("black","blue","blue","orange","orange","black","black"), lty=c(1,1,2,1,2,2,2),
#       lwd=c(3,3,3,3,3,3,3))
legend("bottomleft", c("Experimental QTET", "Unconditional 3 Per.", "Conditional 3 Per.", "Actual Copula Function", "Independence Copula Function"), col=c("black","blue","orange","green","red"), lty=c(1,1,1,2,2), lwd=3)
dev.off()


##make the QTET plot when using the actual change distribution
par(mfrow=c(1,1)) #reset plot layout
png("~/Documents/school/projects/Common App/paper/figures/change-test.png")
plot(probs,actual.qte.employed,type="l", ylim=c(-8000,8000), lwd=3, main="Estimated QTETs",
     xlab="quantile", ylab="QTE")
lines(probs,lalonde.fy3$qte, col="blue", lwd=3)
#lines(probs,lalonde.fy3.changetest$qte, col="blue", lty=2, lwd=3)
lines(probs,lalonde.fy3.cov$qte, col="orange", lwd=3)
lines(probs,lalonde.fy3.changetest$qte, col="green", lty=2, lwd=3)
#lines(probs, lalonde.fy$ub.qte, lty=2, lwd=3)
#lines(probs, lalonde.fy$lb.qte, lty=2, lwd=3)
#thuysbaert
#...
#legend("bottomleft", c("Experimental QTE","3 Per. Observational Change","3 Per. Actual Change Data","3 Per. w Covariates Observational Change", "3 Per. w Covs & Actual Change Data", "Fan-Yu Bounds", ""), 
#       col=c("black","blue","blue","orange","orange","black","black"), lty=c(1,1,2,1,2,2,2),
#       lwd=c(3,3,3,3,3,3,3))
legend("bottomlef", c("Experimental QTET", "Unconditional 3 Per.", "Conditional 3 Per.", "Actual Untreated Change"), col=c("black","blue","orange","green"), lty=c(1,1,1,2), lwd=3)
dev.off()

##make the QTET plot when using the actual initial distribution
par(mfrow=c(1,1)) #reset plot layout
png("~/Documents/school/projects/Common App/paper/figures/initial-test-2.png")
plot(probs,actual.qte.employed,type="l", ylim=c(-25000,12000), lwd=3, main="Estimated QTETs",
     xlab="quantile", ylab="QTE")
lines(probs,lalonde.fy3$qte, col="blue", lwd=3)
lines(probs,lalonde.fy3.initialtest$qte, col="blue", lty=2, lwd=3)
lines(probs,lalonde.fy3.cov$qte, col="orange", lwd=3)
lines(probs,lalonde.fy3.cov.initialtest$qte, col="orange", lty=2, lwd=3)
lines(probs, lalonde.fy$ub.qte, lty=2, lwd=3)
lines(probs, lalonde.fy$lb.qte, lty=2, lwd=3)
#thuysbaert
#...
legend("bottomleft", c("Experimental QTE","3 Per. Observational Initial","3 Per. Actual Initial Data","3 Per. w Covariates Observational Initial", "3 Per. w Covs & Actual Initial Data", "Fan-Yu Bounds", ""), 
       col=c("black","blue","blue","orange","orange","black","black"), lty=c(1,1,2,1,2,2,2),
       lwd=c(3,3,3,3,3,3,3))
dev.off()


##make the QTET plot for the simulation test
par(mfrow=c(1,1)) #reset plot layout
png("~/Documents/school/projects/Common App/paper/figures/sim-test-2.png")
plot(probs,actual.qte.employed,type="l", ylim=c(-25000,12000), lwd=3, main="Estimated QTETs",
     xlab="quantile", ylab="QTE")
lines(probs,lalonde.fy3$qte, col="blue", lwd=3)
lines(probs,lalonde.fy3.simtest$qte, col="blue", lty=2, lwd=3)
lines(probs,lalonde.fy3.cov$qte, col="orange", lwd=3)
lines(probs,lalonde.fy3.cov.simtest$qte, col="orange", lty=2, lwd=3)
lines(probs, lalonde.fy$ub.qte, lty=2, lwd=3)
lines(probs, lalonde.fy$lb.qte, lty=2, lwd=3)
#thuysbaert
#...
legend("bottomleft", c("Experimental QTE","3 Per. Observational Initial","3 Per. Actual Initial Data","3 Per. w Covariates Observational Initial", "3 Per. w Covs & Actual Initial Data", "Fan-Yu Bounds", ""), 
       col=c("black","blue","blue","orange","orange","black","black"), lty=c(1,1,2,1,2,2,2),
       lwd=c(3,3,3,3,3,3,3))
dev.off()


##Finally, plot all the differences (covariates case)
par(mfrow=c(1,1)) #reset plot layout
png("~/Documents/school/projects/Common App/paper/figures/all-test-cov.png")
plot(probs,actual.qte.employed,type="l", ylim=c(-25000,12000), lwd=3, main="Estimated QTETs",
     xlab="quantile", ylab="QTE")
lines(probs,lalonde.fy3.cov$qte, col="orange", lwd=3)
lines(probs,lalonde.fy3.cov.copulatest$qte, col="red", lty=2, lwd=3)
lines(probs,lalonde.fy3.cov.changetest$qte, col="green", lwd=3, lty=2)
lines(probs,lalonde.fy3.cov.initialtest$qte, col="yellow", lty=2, lwd=3)
lines(probs,lalonde.fy3.cov.simtest$qte, col="purple", lty=2, lwd=3)
#lines(probs, lalonde.fy$ub.qte, lty=2, lwd=3)
#lines(probs, lalonde.fy$lb.qte, lty=2, lwd=3)
#thuysbaert
#...
legend("bottomleft", c("Experimental QTE","Estimated QTE","Copula Test","Change Test","Initial Test", "Sim Test"), 
       col=c("black","orange","red","green","yellow","purple"), lty=1,
       lwd=3)
dev.off()


##Plot based on what happens based on different h
par(mfrow=c(1,1)) #reset plot layout
#png("~/Documents/school/projects/Common App/paper/figures/bandwidth-test.png")
plot(probs,actual.qte.employed,type="l", ylim=c(-25000,12000), lwd=3, main="Estimated QTETs",
     xlab="quantile", ylab="QTE")
lines(probs,lalonde.fy3.simtest$qte, col="blue", lwd=3)
lines(probs,lalonde.fy3.simtest.02$qte, col="orange", lwd=3)
lines(probs,lalonde.fy3.simtest.05$qte, col="red", lty=2, lwd=3)
lines(probs,lalonde.fy3.simtest.1$qte, col="green", lwd=3, lty=2)
lines(probs,lalonde.fy3.simtest.2$qte, col="yellow", lty=2, lwd=3)
lines(probs,lalonde.fy3.simtest.3$qte, col="purple", lty=2, lwd=3)
lines(probs,lalonde.fy3.simtest.4$qte, col="gray", lty=2, lwd=3)
lines(probs, lalonde.fy$ub.qte, lty=2, lwd=3)
lines(probs, lalonde.fy$lb.qte, lty=2, lwd=3)
#thuysbaert
#...
legend("bottomleft", c("Experimental QTE","changed v to full range bw=0.3", "bw=0.02","bw=0.05","bw=0.1","bw=0.2", "bw=0.3", "bw=0.4"), 
       col=c("black","blue","orange","red","green","yellow","purple","gray"), lty=1,
       lwd=3)
#dev.off()


##calculate the summary statistics and plot
lalonde.exp$id <- 1
lalonde.exp$cat <- 1
lalonde.psid$cat <- 2*(1-lalonde.psid$treat)
lalonde.all <- rbind(subset(lalonde.exp, treat==0), lalonde.psid)
require(tabular)
require(Hmisc)

tabular( (re78+re75+re74+u75+u74+age+education+black+hispanic+married+nodegree) ~ (n=1) + (Format(digits=2))*(Category=as.factor(cat))*(mean+sd), data=lalonde.all)

#calculate the normalized difference
nd <- function(xvec, yvec) {
    nx <- length(xvec)
    ny <- length(yvec)
    num <- mean(xvec) - mean(yvec)
    denom <- sqrt(var(xvec) + var(yvec))
    num/denom
}

#specialize normalized difference method to work
#for the lalonde dataset called by tabular
lalonde.nd <- function(xvec, yvec)
    #note that we want to calculate everything relative to the treated group
    #so we want to pass only the treated obs on to nd method in yvec
    nd(xvec, yvec[which(lalonde.all$cat==0)])
}

#this creates a table for wages by treated, randomized, observational
latex(tabular( (re78 + re75 + re74) ~
Format(digits=2)*(Treated=(cat==0))*(mean+sd) +
    Format(digits=2)*(Randomized=(cat==1))*(mean+sd+(nd=(Percent(fn=lalonde.nd)))) +
    Format(digits=2)*(Observational=(cat==2))*(mean+sd+(nd=(Percent(fn=lalonde.nd)))),
data=lalonde.all))

#this creates a table for covariates by treated, randomized, observational
latex(tabular( (u75+u74+age+education+black+hispanic+married+nodegree) ~
Format(digits=2)*(Treated=(cat==0))*(mean+sd) +
    Format(digits=2)*(Randomized=(cat==1))*(mean+sd+(nd=(Percent(fn=lalonde.nd)))) +
    Format(digits=2)*(Observational=(cat==2))*(mean+sd+(nd=(Percent(fn=lalonde.nd)))),
data=lalonde.all))
