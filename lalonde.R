require(causalsens)

data(lalonde.exp)
data(lalonde.psid)


#calculate the average treatment effect on the treated without
#covariates
att = mean(subset(lalonde.exp,treat==1)$re78) - mean(subset(lalonde.exp,treat==0)$re78)

#parameters used below:
#h is the step size passed to various function
h = 0.3
#probs are the values that we compute QTET for
probs = seq(h,1-h,0.02)

#experimental QTE
actual.qte = quantile(subset(lalonde.exp,treat==1)$re78,probs=probs) - quantile(subset(lalonde.exp,treat==0)$re78,probs=probs)
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

#call firpo (for cross-sectional case) method
lalonde.firpo = firpo(re ~ treat, x=c("age","education","black","hispanic",
                      "married","nodegree","u74","u75"), data=lalonde.data,
      probs=seq(0.02,0.98,0.02))

#ptm = proc.time()
#Rprof()
#call fan-yu for bounds
lalonde.fy = fan.yu(re ~ treat,
               tname="year",t=1978, tmin1=1975, data=lalonde.data,
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
lalonde.fy3 = threeperiod.fanyu(re ~ treat,
                           tname="year",t=1978, tmin1=1975, tmin2=1974,
                           data=lalonde.data, idname="id", uniqueid="uniqueid",
                           #x=c("age","education","black","hispanic",
                           #     "married","nodegree","u74","u75"),
                           x=NULL,
                           #y.seq=seq(0,120000,length.out=20),
                           #dy.seq=seq(-70000,120000,length.out=20), 
                           y.seq=sort(unique(lalonde.exp$re78)),
                           dy.seq=sort(unique(lalonde.exp$re78-lalonde.exp$re75)),
                           probs=probs,
                           dropalwaystreated=FALSE,
                           h=h, probevals=200)
#Rprof(NULL)
#summaryRprof()
#proc.time()-ptm

#ptm = proc.time()
#Rprof()
#call panelDid with covariates
lalonde.fy3.cov = threeperiod.fanyu(re ~ treat,
                                tname="year",t=1978, tmin1=1975, tmin2=1974,
                                data=lalonde.data, idname="id", uniqueid="uniqueid",
                                x=c("age","education","black","hispanic",
                                "married","nodegree","u74","u75"),
                                #x=NULL,
                                y.seq=sort(unique(lalonde.exp$re78)),
                                dy.seq=sort(unique(lalonde.exp$re78-lalonde.exp$re75)),
                                probs=probs,
                                dropalwaystreated=FALSE,
                                h=0.1, probevals=100)
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
plot(density(subset(lalonde.exp,treat==1)$re75 - subset(lalonde.exp,treat==1)$re74), lwd=3, main="Change in Untreated Outcomes")
lines(density(subset(lalonde.exp,treat==0)$re75 - subset(lalonde.exp,treat==0)$re74), col="blue", lwd=3)
legend("topleft", c("Change Untreated Outcomes Treated Group 1975-1974", "Change Untreated Outcomes Untreated Group 1975-1974"),
       col=c("black","blue"), lty=c(1,1), lwd=c(3,3))

#KS test can test if the two distributions are the same
#ks.test(log(d05.treated.04$applcn),log(d05.untreated.04$applcn))
ks = ks.test((subset(lalonde.exp,treat==1)$re75 - subset(lalonde.exp,treat==1)$re74),
        (subset(lalonde.exp,treat==0)$re75 - subset(lalonde.exp,treat==0)$re74))


#2.b) QTETs
par(mfrow=c(1,1)) #reset plot layout
plot(probs,actual.qte,type="l", ylim=c(-25000,12000), lwd=3, main="Estimated QTETs",
     xlab="quantile", ylab="QTE")
lines(probs,lalonde.fy3$qte, col="blue", lwd=3)
#lines(probs,lalonde.fy3.cov$qte, col="purple", lwd=3)
lines(probs,lalonde.ai$qte, col="red", lwd=3)
#randomization
lines(probs,
      quantile(subset(lalonde.psid,treat==1)$re78,probs=probs) - 
        quantile(subset(lalonde.psid,treat==0)$re78,probs=probs),
      col="green", lwd=3)
lines(probs, lalonde.firpo$qte, col="purple", lwd=3)
lines(probs, lalonde.fy$ub.qte, lty=2, lwd=3)
lines(probs, lalonde.fy$lb.qte, lty=2, lwd=3)
#thuysbaert
#...
legend("bottomleft", c("Experimental QTE","3 Per.","AI","Random Treatment", "Firpo", "Fan-Yu Bounds"), 
       col=c("black","blue","red","green","purple","black","black"), lty=c(1,1,1,1,1,2,2),
       lwd=c(3,3,3,3,3,3,3))


#3) Covariates case
#3.a) First plot the counterfactual distributions from the covariates case, the no covariates case, and the experimental case
plot(lalonde.fy3$F.treatedcf.t, col="blue", main="Counterfactual Untreated Distribution of Wages")
lines(lalonde.fy3.cov$F.treatedcf.t, col="green")
lines(ecdf(subset(lalonde.exp,treat==0)$re78),col="black")
legend("bottomright", c("Experimental","3 Per. No Covariates","3 Per. Covariates"),
       col=c("black","green","blue"), lty=c(1,1,1),
       lwd=c(3,3,3))
