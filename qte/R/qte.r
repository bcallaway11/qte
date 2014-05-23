#TODO: check which of these packages are actually needed
#require(Hmisc)
#require(gmm)
#require(KernSmooth)
require(np)
require(MASS)
require(Rcpp)
require(copula)

#use this line for prototyping; comment for building package
sourceCpp("~/Documents/school/projects/Common App/data/qte/src/qtecpp.cpp")

#####MAIN FUNCTIONS#####

#athey-imbens
#formula should be outcome variable y on treatment variable
# e.g. I(log(applcn)) ~ joined.dummy
# t is the post-treatment time period
# tmin1 is the pre-treatment time period
# this setup only accomodates two time periods; hopefully,
# I will extend this to the full panel case at some point.
#tname is the name of the time variable as a string
#x are covariates - this does not work yet.
#data is a data frame
#dropalwaystreated indicates how to handle observations
# that begin treatment before the pre-treatment period;
# this is not an issue in the original Athey-Imbens setup
# because treatment can only occur in period t; however, 
# in the case of Common App, some colleges are treated
# before the sample period begins.

#TODO: accommodate covariates
#'@title athey.imbens
#' @inheritParams threeperiod.fanyu
#' @export
athey.imbens <- function(formla, t, tmin1, tname, x=NULL,data,y.seq,
                          dropalwaystreated=TRUE, idname, uniqueid, probs) {
  form = as.formula(formla)
  dta = model.frame(terms(form,data=data),data=data) #or model.matrix
  colnames(dta) = c("y","treatment")
  yname="y"
  treat="treatment"
  data=cbind.data.frame(dta,data)
  
  #drop the always treated.  Note that this also relies
  #on no "switchback" or no return to untreated status
  #after joining treatment.
  #first line gets the correct two years of data
  data = subset(data, (data[,tname]==tmin1 | data[,tname]==t))
  data = makeBalancedPanel(data, idname, tname)
  
  #TODO: THIS DOESN'T WORK
  if (dropalwaystreated) {
    #now, drop the observations that are initially treated
    #data$lagged.treat = getLag(xname=treat, idname=idname,
    #                           tname=tname, data=data, 
    #                           withBetweenYears=FALSE)
    #data = subset(data, !(data[,treat]==1 & data[,tname]==tmin1))
    #data = subset(data, !(data[,lagged.treat]==1 & data[,tname]==t))
    #below will be application specific, but I NEED SOMETHING
    # to work.
    data = subset(data, !(d03==1 & d10==1))
  }
  
  #just to make sure the factors are working ok
  data = droplevels(data)
  
  #Setup each of the datasets used below
  #a) get all the treated (in the last period) observations
  treated.t = data[data[,tname]==t & data[,treat]==1,]
  
  #b) set ever.treated to 1 if observation is treated in last period
  data$ever.treated = data$treatment
  data$ever.treated = 1*(data[,idname] %in% treated.t[,idname])  
  ever.treated = "ever.treated"
  
  #Setup each of the datasets used below
  #treated.t = subset(data, data[,treat]==1 & data[,tname]==t)
  #just get the lagged value of y; otherwise keep the same
  #dataset.  Note: this will not work if there are x covariates;
  #well, could follow similar procedure, but as is, would
  #require some modification.
  treated.tmin1 = data[ data[,tname] == tmin1 & 
                          data[,ever.treated] == 1, ]
  #this is right
  untreated.t = data[data[,tname]==t & data[,treat]==0,]
  #get lagged of untreated y
  untreated.tmin1 = data[ data[,tname] == tmin1 &
                            data[,ever.treated] == 0, ]
  
  ###First, get distribution Y_1t | Dt=1
  F.treated.t = ecdf(treated.t[,yname])
  
  ###Now, compute the counterfactual distribution
  ###Y_0t | D_t=1.  There are several steps.
  #1) compute Y_0t | D_t=0
  F.untreated.t = ecdf(untreated.t[,yname])
  
  F.untreated.tmin1 = ecdf(untreated.tmin1[,yname])
  
  #2) compute F^-1_untreated.tmin1
  Finv.untreated.tmin1 <- function(ps) {
    return(quantile(untreated.tmin1[,yname],probs=ps))
  }
  ai.inner = Finv.untreated.tmin1(F.untreated.t(y.seq))
  
  #3) compute distribution Y_0tmin | Dt=1
  F.treated.tmin1 = ecdf(treated.tmin1[,yname])
  
  #3a) use this to compute counterfactual distribution
  F.treatedcf.tval = F.treated.tmin1(ai.inner)
  
  F.treatedcf.t = approxfun(y.seq, F.treatedcf.tval, method="constant",
                                 yleft=0, yright=1, f=0, ties="ordered")
  class(F.treatedcf.t) = c("ecdf", "stepfun", class(F.treatedcf.t))
  assign("nobs", length(ai.inner), envir = environment(F.treatedcf.t))
  
  #4) Plot each distribution
  plot(F.treated.tmin1,col="blue")
  lines(F.untreated.tmin1,col="red")
  lines(F.untreated.t,col="green",lty=2) 
  
  plot(y.seq,F.treated.t(y.seq), type="l")
  #note that F.treated.t is an ecdf and can be used as a function,
  #but F.treatedcf.t is a list of calculated points
  # => There is some fragility in the code in this regard,
  # but it does work as long as you make sure F.treatedcf.t
  #is calculated on the same number of points as you pass
  #to F.treated.t above.
  lines(y.seq,F.treatedcf.t(y.seq), lty=2)
  
  #5) Compute Quantiles
  #a) Quantiles of observed distribution
  q1 = quantile(treated.t[,yname],probs=probs)
  q0 = quantile(F.treatedcf.t,probs=probs)
  
  #6) Plot QTE
  plot(probs, q1-q0, type="l")
  
  #7) Estimate ATT using diff-in-diff
  #straight d-i-d estimate of att
  #make sure that this is right -- esp. what is going on with
  # those that were already treated at t-1
  att = mean(treated.t[,yname]) - (mean(treated.tmin1[,yname])
    + (mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname])))
  #add this to the plot
  abline(a=att, b=0)
  return(list(F.treated.t = F.treated.t, F.treatedcf.t = F.treatedcf.t,
              att=att, qte=(q1-q0),
              treated.t=treated.t, treated.tmin1=treated.tmin1, 
              untreated.t=untreated.t, untreated.tmin1=untreated.tmin1))
}


####Bounds with Fan-yu
##Function to implement the fan-yu procedure
##y.seq is the support of y as a sequence
##dy.seq is the support of the change in y as a sequence
#' @title fan.yu
#' @inheritParams threeperiod.fanyu
#' @export
fan.yu <- function(formla, t, tmin1, tname, x=NULL,data,
                   dropalwaystreated=TRUE, idname, uniqueid,
                   y.seq, dy.seq, probs=seq(0,1,0.1)) {
  form = as.formula(formla)
  dta = model.frame(terms(form,data=data),data=data) #or model.matrix
  colnames(dta) = c("y","treatment")
  yname="y"
  treat="treatment"
  data=cbind.data.frame(dta,data)
  
  #drop the always treated.  Note that this also relies
  #on no "switchback" or no return to untreated status
  #after joining treatment.
  
  data = subset(data, (data[,tname]==tmin1 | data[,tname]==t))
  data = makeBalancedPanel(data, idname, tname)
  
  if (dropalwaystreated) {
    #first line gets the correct two years of data
    #now, drop the observations that are initially treated
    #data$lagged.treat = getLag(xname=treat, idname=idname,
    #                           tname=tname, data=data, 
    #                           withBetweenYears=FALSE)
    #data = subset(data, !(data[,treat]==1 & data[,tname]==tmin1))
    #data = subset(data, !(data[,lagged.treat]==1 & data[,tname]==t))
    #below will be application specific, but I NEED SOMETHING
    # to work.
    data = subset(data, !(d06==1 & d10==1))
  }
  
  #just to make sure the factors are working ok
  data = droplevels(data)
  
  
  #Setup each of the datasets used below
  #treated.t = subset(data, data[,treat]==1 & data[,tname]==t)
  #just get the lagged value of y; otherwise keep the same
  #dataset.  Note: this will not work if there are x covariates;
  #well, could follow similar procedure, but as is, would
  #require some modification.
  #treated.tmin1 = data
  #treated.tmin1$y = getLag(xname=yname, idname=idname,
  #                         tname=tname,data=data,withBetweenYears=FALSE)
  #treated.tmin1 = subset(treated.tmin1, data[,treat]==1 & data[,tname]==t)
  #this is right
  #untreated.t = subset(data, data[,treat]==0 & data[,tname]==t)
  #get lagged of untreated y
  #untreated.tmin1 = data
  #untreated.tmin1$y = getLag(xname=yname, idname=idname,
  #                           tname=tname, data=data,withBetweenYears=FALSE)
  #untreated.tmin1 = subset(untreated.tmin1, data[,treat]==0 & data[,tname]==t)

  #a) get all the treated (in the last period) observations
  treated.t = data[data[,tname]==t & data[,treat]==1,]

  #b) set ever.treated to 1 if observation is treated in last period
  data$ever.treated = data$treatment
  data$ever.treated = 1*(data[,idname] %in% treated.t[,idname])  
  ever.treated = "ever.treated"

  #Setup each of the datasets used below
  #treated.t = subset(data, data[,treat]==1 & data[,tname]==t)
  #just get the lagged value of y; otherwise keep the same
  #dataset.  Note: this will not work if there are x covariates;
  #well, could follow similar procedure, but as is, would
  #require some modification.
  treated.tmin1 = data[ data[,tname] == tmin1 & 
                        data[,ever.treated] == 1, ]
  #this is right
  #untreated.t = subset(data, data[,treat]==0 & data[,tname]==t)
  untreated.t = data[data[,tname]==t & data[,treat]==0,]
  #get lagged of untreated y
  untreated.tmin1 = data[ data[,tname] == tmin1 &
                          data[,ever.treated] == 0, ]
  
  
  ###First, get distribution Y_1t | Dt=1
  F.treated.t = ecdf(treated.t[,yname])
  
  ##2) Compute the lower bound; as currently constructed,
  #this will use very few number of observations
  #2a) Recover the untreated outcomes for the newly treated
  # in the previous period
  #newly.treated.t <<- subset(treated.t, 
   #                          !(treated.t[,idcolnames] %in% treated.tmin1[,idcolnames]))
  #newly.treated.tmin1 <<- subset(data,
   #                              data[,idcolnames] %in% newly.treated.t[,idcolnames] & 
    #                               data[,tname] == tmin1)
  F.treated.tmin1 = ecdf(treated.tmin1[,yname]) #as long as 
  ##dropalwaystreated==TRUE, this distribution is the distribution
  ##of outcomes for the treated group before they were treated.
  
  #2b) Get the distribution of the change in outcomes for the newly treated
  #F.change.newlytreated <<- ecdf(newly.treated.t[,yname] - newly.treated.tmin1[,yname])
  F.treated.change = ecdf(treated.t[,yname]-treated.tmin1[,yname])
  #Actually -- don't think you need that...
  
  #2c) Get the distribution of the change in outcomes for the never treated
  F.untreated.change = ecdf(untreated.t[,yname]-untreated.tmin1[,yname])
  
  #for comparison, compute att first
  att = mean(treated.t[,yname]) - mean(treated.tmin1[,yname]) -
    (mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname]))
  #2c.1) If there are covariates, then above distribution needs to be changed
  if (!(is.null(x))) {
   treated.t$changey = treated.t[,yname] - treated.tmin1[,yname]
   untreated.t$changey = untreated.t[,yname] - untreated.tmin1[,yname]
   pscore.data = rbind(treated.t, untreated.t)
   xmat = pscore.data[,x]
   pscore.reg = glm(pscore.data[,treat] ~ as.matrix(xmat),
                    family=binomial(link="probit"))
   pscore = fitted(pscore.reg)
   pD1 = nrow(treated.t)/nrow(untreated.t)
   p.dy.seq = unique(pscore.data$changey)
   #F.untreated.change = rep(0,length(p.dy.seq))
   distvals = rep(0,length(p.dy.seq))
   for (dy in p.dy.seq) {
     distvals[which(dy==p.dy.seq)] = mean(1*(pscore.data$changey<=dy)*
                                                   (1-pscore.data[,treat])*pscore/((1-pscore)*pD1))
   }
   pscore.data$distvals = distvals
   #instead try to do something like ecdf: cf. ecdf, approxfun, and cumsum
   #documentation
   pscore.data1 = pscore.data[order(pscore.data$changey),]
   #changey = pscore.data1$changey
   #vals = unique(changey)
   #tempfun = approxfun(vals, cumsum(pscore*tabulate(match(changey, vals))/(nrow(pscore.data1))),
             #method="constant", yleft=0, yright=1, f=0, ties="ordered")
   F.untreated.change = approxfun(pscore.data1$changey, pscore.data1$distvals, method="constant",
                       yleft=0, yright=1, f=0, ties="ordered")
   class(F.untreated.change) = c("ecdf", "stepfun", class(F.untreated.change))
   #att using abadie method
   att = mean(((pscore.data$changey)/pD1)*(pscore.data[,treat] - pscore)/(1-pscore))
  }
  
  #2c) Get the distribution of outcomes for the newly treated at (t-1)
  #F.newlytreated.tmin1 <<- ecdf(newly.treated.tmin1[,yname])
  
  #2d) get the lower bound
  #make sure that this is right, but we are taking the smallest
  # over the support (I think) of y
  supy = y.seq #sort(unique(treated.tmin1$y))#this is specific to application
  posvals = dy.seq #these are the values to min over; not sure
   #exactly what they should be, but should cover 0 probably about
   #as wide as the support of y is in each direction
  # and should probably be passed into the function
  # I think that I can pass this as both arguments s, and y below
  # but maybe should separate them esp. if there are issues
  lbs = vapply(supy,FUN=getlb,FUN.VALUE=1, 
               ecdf.change.newlytreated=F.untreated.change,
               ecdf.newlytreated.tmin1=F.treated.tmin1,
               y=posvals)
  ubs = vapply(supy,FUN=getub,FUN.VALUE=1, 
               ecdf.change.newlytreated=F.untreated.change,
               ecdf.newlytreated.tmin1=F.treated.tmin1,
               y=posvals)
  #plot the upper and lower bound distributions
  plot(supy, lbs, type="l", lty=2, ylim=c(0,1), xlab="log(Applications)",
       ylab="F")
  lines(supy, ubs, lty=2)
  
  #get upper bound quantiles
  ub.quantiles = simple.quantile(x=supy, Fx=lbs, probs=probs)
  
  #get lower bound quantiles
  lb.quantiles = simple.quantile(x=supy, Fx=ubs, probs=probs)
  
  #plot bounds on qte
  lb.qte = as.numeric(quantile(treated.t[,yname],probs=probs) - 
                        ub.quantiles)
  ub.qte = as.numeric(quantile(treated.t[,yname],probs=probs) - 
                        lb.quantiles)
  plot(probs, lb.qte, 
       type="l", lty=2, xlab="tau", ylab="QTE",
       ylim=c(-2.5,2.5))
       #ylim=c(min(lb.qte[2:(length(lb.qte)-1)],ub.qte[2:(length(ub.qte)-1)]),max(lb.qte,ub.qte)))
  lines(probs, ub.qte, lty=2)
  abline(a=att, b=0, col="blue")
  print(paste("ATT: ",att))
  
  return(list(lbs=lbs,ubs=ubs, ub.quantiles=ub.quantiles,
              lb.quantiles=lb.quantiles, ub.qte=ub.qte,
              lb.qte = lb.qte, att=att))
}


#####Three Period Fan-Yu#####
#Idea here is that we can use information from a third period
#to point identify counterfactual distribution of outcomes
#for the treated group
#Change the name of this function: "panelDiD" e.g.
#main function should take in x,y,d
#call plot function, summary function, formula function, etc. later
#add functionality to pass in pscore
#' @title Three Period Fan Yu
#'
#' @description
#' \code{threeperiod.fanyu} uses third period of data, combined with Distributional Difference in Differences assumption (Fan and Yu, 2012) to point identify QTET.
#' 
#' @param formla outcome variable on treatment
#' @param t last time period
#' @param tmin1 middle time period
#' @param tmin2 initial time period
#' @param x additional covariates if using propensity score reweighting technique
#' @param dropalwaystreated boolean indicating whether in true panel data context (when treatment can occur in any period) whether or not previously treated observations should be dropped from the sample.  This functionality is not currently implemented
#' @param idname an id that identifies individual units over time
#' @param uniqueid an id that identifies unique 'observations'; i.e. each individual-time observation should have its own uniqueid
#' @param zname not used / needs to be deleted
#' @param y.seq the support of y (need to implement functionality so that this is not required)
#' @param dy.seq the support of the change in y (need to implement functionality so that this is not required)
#' @param probs the values at which the quantile treatment effects should be computed
#' @param h the step size used in the simulated integral calculation
#' @param bootstrap.iter boolean passed that is passed in when this method is used to compute standard errors
#' @param copBool boolean that indicates whether parametric copula should be used.  Default FALSE indicates that nonparametric copula should be used.
#'
#' @return threeperiod object
#'
#' @examples
#' \dontrun{
#' threeperiod.fanyu(y~d)
#' }
#' @export
threeperiod.fanyu <- function(formla, t, tmin1, tmin2,
                              tname, x=NULL,data, 
                              dropalwaystreated=TRUE, idname, uniqueid, zname,
                              y.seq=NULL, dy.seq=NULL, probs=seq(0,1,0.1),
                              h, probevals, bootstrap.iter=FALSE, copBool=FALSE) {
  form = as.formula(formla)
  dta = model.frame(terms(form,data=data),data=data) #or model.matrix
  colnames(dta) = c("y","treatment")
  yname="y"
  treat="treatment"
  data=cbind.data.frame(dta,data)
  
  #drop the always treated.  Note that this also relies
  #on no "switchback" or no return to untreated status
  #after joining treatment.
  #TODO: this doesn't work, can move some duplicate functionality
  #across methods to outside function.
  if (dropalwaystreated & !bootstrap.iter) {
    #first line gets the correct two years of data
    data = subset(data, (data[,tname]==tmin1 | data[,tname]==t | 
                           data[,tname]==tmin2))
    data = makeBalancedPanel(data, idname, tname)
    #now, drop the observations that are initially treated
    #data$lagged.treat = getLag(xname=treat, idname=idname,
    #                           tname=tname, data=data, 
    #                           withBetweenYears=FALSE)
    #data = subset(data, !(data[,treat]==1 & data[,tname]==tmin1))
    #data = subset(data, !(data[,lagged.treat]==1 & data[,tname]==t))
    #below will be application specific, but I NEED SOMETHING
    # to work.
    data = subset(data, !(d03==1 & d10==1))
    data = subset(data, !(d06==1 & d10==1))
    #probably want to make above line drop those who are treated
    #in the middle period, but leave for now.
  }
  
  #just to make sure the factors are working ok
  data = droplevels(data)
  
  #1) set up a dummy variable indicating whether the individual is 
  #treated in the last period.
  
   #a) get all the treated (in the last period) observations
   treated.t = data[data[,tname]==t & data[,treat]==1,]
   
   #b) set ever.treated to 1 if observation is treated in last period
   data$ever.treated = data$treatment
   data$ever.treated = 1*(data[,idname] %in% treated.t[,idname])  
   ever.treated = "ever.treated"
   treated.t$ever.treated = 1
  
  #Generate subsets of the panel data based on time period and
  #whether or not the observation is treated.  These will be used
  #to calculate distributions below.  
  treated.tmin1 = data[ data[,tname] == tmin1 & 
                           data[,ever.treated] == 1, ]
  #treated at t-2
  treated.tmin2 = data[ data[,tname] == tmin2 &
                           data[,ever.treated] == 1, ]
  
  #untreated at t
  untreated.t = data[data[,tname]==t & data[,treat]==0,]
  
  #untreated at t-1 & t-2
  untreated.tmin1 = data[ data[,tname] == tmin1 &
                             data[,ever.treated] == 0, ]
  untreated.tmin2 = data[ data[,tname] == tmin2 &
                             data[,ever.treated] == 0, ]
 
  #3) First, generate the copula function using the (t-2) time period.
  # The copula function will capture the dependence structure between
  # the untreated level of outcomes in period (t-2) and the change in
  # untreated outcomes between (t-1) and (t-2) for the eventually treated
  # observations.
  
   #a) Get distribution of y0.tmin2 | Dt=1
   F.treated.tmin2 = ecdf(treated.tmin2[,yname])
  
   #b) Get distribution of y0.tmin1 - y0tmin2 | Dt=1
   F.treated.change.tmin1 = ecdf(treated.tmin1[,yname] -
                                   treated.tmin2[,yname])
  
   #b.1) Get correlations between level at (t-1) and change between
   # (t-1) and (t-2).  These are not used, but returned by the function.
   pearsons.cor = cor(treated.tmin1[,yname] - treated.tmin2[,yname],
                      treated.tmin2[,yname])
   spearmans.rho = cor(treated.tmin1[,yname] - treated.tmin2[,yname],
                       treated.tmin2[,yname],method="spearman")
   kendalls.tau = cor(treated.tmin1[,yname] - treated.tmin2[,yname],
                      treated.tmin2[,yname],method="kendall")
   #OLD: below is up for deletion possibly
   #spearmans.rho2=cor(F.treated.change.tmin1(seq(-2,2,length.out=100)),
   #                   F.treated.tmin2(seq(5,11,length.out=100)))
   #see what happens to the above if you use the actual points in the
   #distribution rather than feeding in extra point; e.g. a lot of 
   #simulataneous zeros would tend to artificially increase the
   #estimated correlation coefficient.
  
   #c) Get joint distribution of y0.tmin2 | Dt=1 and
    # y0.tim1 - y0tmin2 | Dt=1
   #this will return the joint distribution as a function
   #the joint distribution between the change and the initial
   # is the key piece of informationt that we need to recover the copula
   F.joint.tmin1 = getJoint(as.numeric(treated.tmin1[,yname] - 
                              treated.tmin2[,yname]),
                            as.numeric(treated.tmin2[,yname]))
   
   
  #d) Get the copula function
   #internal function returns value of copula for two values (u,v)
  #between 0 & 1
  #account for dependence outside function
  #@title copula
  #@description this function uses the joint distribution computed just above
  #and is just and application of Sklar's Theorem.
  #@param u NumericVector between 0 and 1
  #@param v NumericVector between 0 and 1, same length as u (no don't have to be)
  #@return The value of the copula function (between 0 & 1)
  copula1 = function(u,v) {

    #if (length(u) != length(v)) {
    #    stop("copula: u and v must be of the same length")
    #}

    if ( any(u>1) | any(v>1) | any(u<0) | any(v<0) ) {
        stop("all arguments must be between 0 and 1")
    }
    #essentially F.joint.tmin1 is fast (approx 1 second used)
    #but quantile is slow (about 15 seconds)
    #temp1 = quantile(F.treated.change.tmin1,u)
    #temp2 = quantile(F.treated.tmin2,v)
    
    #c++ method is much faster; it takes about 2 seconds compared to 15
    #Note: These probably do some interpolating, but interpolation
    # should not cause boundary problem.
    change.quant = quantileCPP(treated.tmin1[,yname]-treated.tmin2[,yname],u)
    initial.quant  = quantileCPP(treated.tmin2[,yname],v)
    #F.joint.tmin1(quantile(F.treated.change.tmin1,u),
    #              quantile(F.treated.tmin2,v))

    F.joint.tmin1(change.quant,initial.quant) #old, this is the straightforward
    # way but it was causing some problems
    #try just putting this between the bounds...
    #Fval <- F.joint.tmin1(change.quant, initial.quant)
    #ub <- apply(data.frame(u=u,v=v), FUN=min, MARGIN=1)
    #lb <- apply(data.frame(u=u,v=v), FUN=function(x) max(x[1] + x[2] - 1, 0),
    #            MARGIN=1)
    #lb*(Fval < lb) + ub*(Fval > ub) + Fval*(Fval>=lb & Fval<=ub)
  }

  #alternatively, get the copula function using the np package npcopula function
  change <- treated.tmin1[,yname]-treated.tmin2[,yname]
  initial <- treated.tmin2[,yname]
  copula.bws <- npudistbw(~ change + initial)
  u.seq = v.seq = seq(0,1,0.02)
  copula.fun <- npcopula(bws=copula.bws,
                         data=data.frame(change=change, initial=initial),
                         u=cbind(u.seq,v.seq))

  #copula should either
  #(i) take in scalar u and vector v and return vector of length(v)
  #(ii) take in vector u and vector v (of same length) and return
  # vector of length(v)
  copula <- function(u,v) {

      if ( (length(u) != 1) & (length(u) != length(v)) ) {
          stop("u must either be scalar or same length as v")
      }

      #when u is scalar will call this function with vapply
      #call with x=v and z=u
      copula.inner <- function(x,z) {
          which.min((z-copula.fun$u1)^2 + (x-copula.fun$u2)^2)
      }

      #when u is vector will call this function with apply
      #call with x=matrix(u,v)
      copula.multiple <- function(x) {
          which.min((x[1]-copula.fun$u1)^2 + (x[2]-copula.fun$u2)^2)
      }

      if (length(u) == 1) {
          
          whichys <- vapply(v, FUN=copula.inner, FUN.VALUE=1, z=u)
          
      } else {

          whichys <- apply(cbind(u,v), FUN=copula.multiple, MARGIN=1)

      }
      
       #find the closest point in copula.fun and take its value
      #whichys <- vapply(v,
      #                  FUN=function(x) which.min((u-copula.fun$u1)^2 +
      #                      (x-copula.fun$u2)^2),
      #                  FUN.VALUE=1)
      #whichval <- (u-copula.fun$u1)^2 + (v-copula.fun$u2)^2
      #whichy <- which.min(whichval)

      copula.fun$copula[whichys]
  }
    
  #4) Get known distributions for period t that we will
  # apply the copula function that we just computed to. 
  
  #a) change for the untreated between t-1 and t
  F.untreated.change = ecdf(untreated.t[,yname]-untreated.tmin1[,yname])
  
  #for comparison, compute att: first with no covariates
  att = mean(treated.t[,yname]) - mean(treated.tmin1[,yname]) -
    (mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname]))
  
  #a.1) If there are covariates need to satisfy the Distributional D-i-D
  #then we will need to modify the distribution of the changes in outcomes
  #using the method presented in the paper.
  #This section does this.  For most flexibility, the user should
  #be able to pass in the propensity score using estimated using any
  #method that he chooses.  In the case where there are covariates,
  #but no propensity score passed in, then this section estimates
  #a propensity score using a simple probit on each of the variables
  #entered additively.
  #TODO: add functionality to allow propensity score values to be
  #passed in.
  if (!(is.null(x))) {
    treated.t$changey = treated.t[,yname] - treated.tmin1[,yname]
    untreated.t$changey = untreated.t[,yname] - untreated.tmin1[,yname]
    pscore.data = rbind(treated.t, untreated.t)
    xmat = pscore.data[,x]
    pscore.reg = glm(pscore.data[,treat] ~ as.matrix(xmat),
                     family=binomial(link="probit"))
    pscore = fitted(pscore.reg)
    pD1 = nrow(treated.t)/nrow(untreated.t)
    p.dy.seq = pscore.data$changey #unique(pscore.data$changey)
    #F.untreated.change = rep(0,length(p.dy.seq))
    #TODO: What is this?  Need to come up with better name for this variable
    distvals = rep(0,length(p.dy.seq))
    for (dy in p.dy.seq) {
      distvals[which(dy==p.dy.seq)] = mean(1*(pscore.data$changey<=dy)*
                                             (1-pscore.data[,treat])*pscore/((1-pscore)*pD1))
    }
    pscore.data$distvals = distvals
    
    pscore.data1 = pscore.data[order(pscore.data$changey),]
    
    #Note: this is also likely to have some smoothing
    F.untreated.change = approxfun(pscore.data1$changey, pscore.data1$distvals, method="constant",
                                   yleft=0, yright=1, f=0, ties="ordered")
    class(F.untreated.change) = c("ecdf", "stepfun", class(F.untreated.change))
    assign("nobs", length(p.dy.seq), envir = environment(F.untreated.change))
    
    #att using abadie method
    att = mean(((pscore.data$changey)/pD1)*(pscore.data[,treat] - pscore)/(1-pscore))
  }
  
  #b) distribution of outcomes for treated in previous period
  F.treated.tmin1 = ecdf(treated.tmin1[,yname])
  #first compute the correct bandwidth
  #F.treated.tmin1.bw = npudistbw(~treated.tmin1[,yname])
  #F.treated.tmin1 = function(y, bw=F.treated.tmin1.bw) {
  #  return(fitted(npudist(edat=y, bws=bw)))
  #}

  #c) use known copula function to get joint distribution
  F.joint.t = function(x,y) {
    copula(F.untreated.change(x),F.treated.tmin1(y))
  }
  
  #5) Recover the counterfactual distribution:  This involves simulating
  #from the joint distribution of the change and the initial computed in step 4
  #F.joint.t.  To do this we follow the method given in Nelson (2005) which
  #involves taking the partial derivative of the copula function with respect
  #to its first argument, and then simulating.

  #a)#partial1.copula estimates the partial derivative of the
  #copula with respect to the first
  #argument numerically; for a given step size h
  #this function is needed to calculate the function of the joint distribution
  #noted in the paper.
  partial1.copula = function(u,v,h) {
    ##try commenting out below to see if has any effect
    if (u > 1-h) {
      return(rep(1,length(v)))
    } else if (u < h) {
      return(rep(0,length(v)))
      ##
    } else {
      #first two lines trim values above 1
      #out = (copula(u+h,v)-copula(u-h,v))/(2*h)
      #out = 1*(out>1) + out*(out<=1)
      out = partial1CopulaCPP(u,v,h,copula)
      return(out)
    }
  }

  #getPartialQuant is a helper function to simulate from the
  #copula function.  This follows Nelson (2005)'s method for simulating
  #from a joint distribution.
  #@param u vector of random uniforms of length = #probevals
  #@param t vector of random uniforms of length = #probevals; this is only
  #used for randomly choosing quantiles of the distribution below.
  #TODO: this function needs to be cleaned up
  getPartialQuant = function(u,t) {
    #v.seq = seq(0.02,0.98,0.02)
    v.seq = seq(h,1-h,0.02) #TODO (def shouldn't have step size h); maybe pass
    #this as parameter?
    
    #this new way all u to be passed in as a vector
    #funcval is a list of returns from partial1.copula
    #with same length as length(u)
    #funcval = lapply(u,FUN=partial1.copula,v=v.seq,h=0.1)
    #this call takes about 3 seconds; so look to improve eventually
    #funcval is a list containing the value of copula function
    #for each different value of u vector; holding v.seq 
    #fixed for each u.
    funcval = getListPartialQuantCPP(u,v.seq,h,copula)
    
    #this applies approxfun1 to each above list
    #and produces length(u) number of function approximations
    #only about 2.5 seconds; okish
    ###func = lapply(funcval, FUN=approxfun1, x=v.seq, method="constant",
    ###              yleft=0, yright=1, f=0, ties="ordered")
    #func should contain a list of function return values (a vector of
    #length length(v.seq)).
    #below is old, but I know that it works (just slow)
    #funcval = partial1.copula(u,v.seq,h)
    #func = approxfun1(y=funcval, x=v.seq,
    #                 method="constant",
    #                 yleft=0, yright=1, f=0, ties="ordered")
    #func = approxfun(v.seq, funcval,
    #                          method="constant",
    #                          yleft=0, yright=1, f=0, ties="ordered")
    #class(func) = c("ecdf", "stepfun", class(func))
    #assign("nobs", length(v.seq), envir = environment(func))
    #return(quantile(func,probs=t))
    
    #this line changes the class, so that we can call quantile
    ###func = lapply(func, FUN=modClass, 
    ###              classNames=c("ecdf","stepfun","function"),
    ###              num=length(v.seq))
    
    #at this point, we have a list of ecdfs, we need to call
    #quantile on each of them for a variety of t-quantiles
    #each element of list contains the t-quantiles specified
    #for each passed in value of u
    #actually, this computes way too many quantiles, we only need
    #to compute one for each u, but here we compute them for each t
    #that we pass in (this is pretty fast half a second cost only)
    #list.quanties = lapply(func, FUN=quantile, probs=t)
    #3 seconds; this can be improved.
    ###list.quanties.old=list()
    ###list.quanties=list()
    ###for (i in 1:length(t)) {
      ###old, this is probably right?  TODO: HAVE TO MAKE SURE WHICH OF THESE
      ###IS RIGHT; EVERYTHING IN THIS FUNCTION WITH ### PRECEDING IS 
      ###WHAT I WAS DOING BEFORE AND THINKING WAS WORKING.
      ###list.quanties.old[[i]] = quantile(func[[i]],t[i])
      ###list.quanties[[i]] = quantile(funcval[[i]], t[i])
    ###}
    list.quanties = getListQuantilesCPP(funcval, t)
    return(list.quanties)
  }

  ##some helper functions for the nonparametric cdf estimation

  #@description the integrated kernel used for esimating
  #a univariate cdf
  int.kernel <- function(u) {
      return(pnorm(u))
  }

  #@return the kernel density cdf estimate at a particular point y
  F.smoothed <- function(y,obs,bandwidth=NULL) {
      #n <- length(obs)
      #(1/n)*sum(int.kernel((y-obs)/bandwidth))
  }

  #@title leavoneoutF
  leaveoneoutF <- function(y, obs, bw) {
      n <- length(obs)
      looF <- rep(0,n)
      for (i in 1:n) {
          looF[i] <- F.smoothed(y, obs[-i], bw)
      }
      return(looF)
  }
  
  #@title cross.validation
  cv <- function(bw, y, obs) {
      n <- length(obs)
      (1/n)*sum((1*(obs<=y) - leaveoneoutF(y, obs, bw))^2)
  }

  cv.switchargs <- function(y, bw, obs) {
      cv(bw, y, obs)
  }

  cv.global <- function(bw, obs) {
      sum(vapply(y.seq, FUN=cv.switchargs, FUN.VALUE=1, bw=bw, obs=obs))
  }

  #@return the correct bandwidth using cross validation
  bwminfun <- function(obs) {
      return(optimize(f=cv.global, lower=0, upper=max(obs),
               tol=(max(obs)/length(obs)), obs=obs))
  }
  
  
  #internal function that computes the value of the counterfactual 
  #distribution at a particular value of y.
  F.treatedcf.tfun = function(y, bw=NULL) {
    #simchange = quantile(F.untreated.change, probs=randu)
    #siminit = quantile(F.treated.tmin1, probs=randv)
    #distval = simchange + siminit
    #return(sum(1*(simchange+siminit<y))/probevals)
    #instead get a smoothed estimate of this!
    #return(vapply(y, FUN=F.smoothed, FUN.VALUE=1,
    #              obs=simchange+siminit, bandwidth=10000))
    #npudistbw(~distval)
    return(fitted(npudist(edat=y, bws=bw)))
  }

  #This section uses all the above functions to generate some random
  #numbers simulated from the joint distribution.
  
  #seed used for random number generation.  It needs to change
  #based on every bootstrap iteration.
  seed = round(1000*runif(1)) 
  set.seed(seed)
  #do the random draws need to incorporate the copula?
  #not sure that I am doing this right...
  #TODO: this is where I left off
  r.unifu = runif(probevals)
  r.unifu = h*(r.unifu<h) + r.unifu*(r.unifu>=h)
  r.unifu = (1-h)*(r.unifu>(1-h)) + r.unifu*(r.unifu<=(1-h))
  r.unift = runif(probevals)
  #TODO: still stuck here; also not sure if getPartialQuant
  #is working correctly either...
  randv = rep(0,probevals)
  #for (i in 1:probevals) {
  #  randv[i] = getPartialQuant(r.unifu[i],r.unift[i])
  #}

  #this is the call to getPartialQuant that does all of the work
  #in the nonparametric case.
  numtrims=0 #to count number of trims
  if (!copBool) {
      randv = unlist(getPartialQuant(r.unifu, r.unift))
      #for some reason, we get values greater than 1 sometimes;
      #do some trimming here, but record how many times that it occurs
      numtrims = sum(1*(randv>1))
      randv = vapply(randv, FUN=min, FUN.VALUE=1, 1)
  }
  
  #get partial quant takes up about 10 out of 11 seconds
  randu = r.unifu
  
  ###Alternatively, use copula function; comment this out if
  #want to do things nonparametrically.  TODO: pass type of copula
  #as a parameter
  if(copBool) {
   disdat = cbind(treated.tmin1[,yname]-treated.tmin2[,yname],
                  treated.tmin2[,yname])
   normdat = pobs(disdat)
   #estimate copula
   clayton.cop = claytonCopula(3,dim=2)
   fit.clayton = fitCopula(clayton.cop, normdat, method="ml")
   clayton.cop$parameters = fit.clayton$estimate
   clayton.cop@parameters = fit.clayton@estimate
   rdraws = rCopula(probevals,clayton.cop)
   randu = rdraws[,1]
   randv = rdraws[,2]
   ###
  }
  
  #helper function to change class
  modClass = function(obj, classNames, num) {
    class(obj) = classNames
    assign("nobs", num, envir=environment(obj))
    return(obj)
  }
  
  #helper function to change order of arguments
  #to approxfun
  approxfun1 = function(y,x,...) {
    return(approxfun(x=x,y=y,...))
  }

  #write the counterfactual distribution as a function that we can call.
  sim.cf <- quantile(F.untreated.change, probs=randu) +
    quantile(F.treated.tmin1, probs=randv)
  F.treatedcf.t.bws <- npudistbw(~sim.cf)
  #note that this is "artificially making the smallest possible value
  #for the counterfactual = 0 because we never pass it anything smaller
  F.treatedcf.t = approxfun(y.seq, #old puts y.seq here
                           vapply(y.seq,FUN=F.treatedcf.tfun,FUN.VALUE=1,
                                  bw=F.treatedcf.t.bws),
                           method="constant",
                           yleft=0, yright=1, f=0, ties="ordered")
  class(F.treatedcf.t) = c("ecdf", "stepfun", class(F.treatedcf.t))
  assign("nobs", length(y.seq), envir = environment(F.treatedcf.t))

  
  #this is the actual distribution of outcomes for the treated
  F.treated.t = ecdf(treated.t[,yname])
  
  q1 = quantile(treated.t[,yname], probs=probs)
  #q0 = simple.quantile(F.treatedcf.t$xrange, F.treatedcf.t$dist,
  #                     probs=seq(0,0.9,0.1))
  q0 = quantile(F.treatedcf.t, probs=probs)
  
  #TODO: make the plot outside of the function
  plot(probs,(q1-q0),type="l",xlab="tau", ylab="qte")
  abline(a=0,b=0,lty=2)
  abline(a=att, b=0, col="blue")
  
  #TODO: make the summary outside of the function
  print("QTE:")
  print(q1-q0)
  print(paste("ATT: ",att))

  retObj = list(F.treatedcf.t = F.treatedcf.t, F.joint.t=F.joint.t,
                F.joint.tmin1 = F.joint.tmin1, F.treated.t=F.treated.t,
      F.treated.tmin1=F.treated.tmin1, F.treated.tmin2=F.treated.tmin2,
      F.treated.change.tmin1=F.treated.change.tmin1,
      F.untreated.change=F.untreated.change,
      treated.t=treated.t[,yname],
      untreated.t=untreated.t[,yname],
      treated.tmin1=treated.tmin1[,yname],
      treated.tmin2=treated.tmin2[,yname],
      untreated.tmin1=untreated.tmin1[,yname],
      untreated.tmin2=untreated.tmin2[,yname],
      untreated.change=untreated.t[,yname]-untreated.tmin1[,yname],
                copula=copula, pearsons.cor=pearsons.cor, qte=(q1-q0),
                spearmans.rho=spearmans.rho, kendalls.tau=kendalls.tau,
                spearmans.rho2=0, att=att, randu=randu, randv=randv,
                randt=r.unift, h=h, numtrims=numtrims)
  class(retObj) = "panelDiD"
  
  return(retObj)
}


#'@title bootstrap.threeperiod.fanyu
#' @inheritParams threeperiod.fanyu
#' @export
bootstrap.threeperiod.fanyu <- function(formla, t, tmin1, tmin2,
                                        tname, x=NULL,data, 
                                        dropalwaystreated=TRUE, idname, uniqueid, zname,
                                        y.seq=NULL, dy.seq=NULL, probs=seq(0,1,0.1),
                                        iters, h, alp=0.05, probevals) {
  form = as.formula(formla)
  dta = model.frame(terms(form,data=data),data=data) #or model.matrix
  colnames(dta) = c("y","treatment")
  yname="y"
  treat="treatment"
  data=cbind.data.frame(dta,data)
  
  ##This is pasted from Athey-Imbens method!
  #drop the always treated.  Note that this also relies
  #on no "switchback" or no return to untreated status
  #after joining treatment.
  if (dropalwaystreated) {
    #first line gets the correct two years of data
    data = subset(data, (data[,tname]==tmin1 | data[,tname]==t) |
                    data[,tname]==tmin2)
    data = makeBalancedPanel(data, idname, tname)
    #now, drop the observations that are initially treated
    #data$lagged.treat = getLag(xname=treat, idname=idname,
    #                           tname=tname, data=data, 
    #                           withBetweenYears=FALSE)
    #data = subset(data, !(data[,treat]==1 & data[,tname]==tmin1))
    #data = subset(data, !(data[,lagged.treat]==1 & data[,tname]==t))
    #below will be application specific, but I NEED SOMETHING
    # to work.
    data = subset(data, !(d03==1 & d10==1))
    data = subset(data, !(d06==1 & d10==1))
  }
  
  #just to make sure the factors are working ok
  data = droplevels(data)
  ##
  #first calculate the actual estimate
  fy3 = threeperiod.fanyu(formla=formla, t=t, tmin1=tmin1, tmin2=tmin2, 
                          tname=tname, x=x,data=data,
                          dropalwaystreated=dropalwaystreated,
                          idname=idname, uniqueid=uniqueid, zname=zname,
                          y.seq=y.seq, dy.seq=dy.seq, probs=probs, h=h,
                          probevals=probevals)
  #now calculate the bootstrap confidence interval
  eachIter = list()
  #Need to build dataset by sampling individuals, and then
  #taking all of their time periods
  all.ids = unique(data[,idname])
  n = length(all.ids)
  out.bootdatalist <<- list()
  for (i in 1:iters) {
    #reset boot.data
    out.iter <<- i
    boot.data = data[0,]
    randy = sample(1:n, n, replace=T)
    #there has to be a way to do this faster, but go with the loop
    #for now
    for (j in all.ids[randy]) {
      boot.data = rbind(boot.data, data[(data[,idname]==j),])
    }
    boot.data = process.bootdata(boot.data, idname, uniqueid)
    out.bootdatalist[[i]] <<- boot.data
    thisIter = threeperiod.fanyu(formla=formla, t=t, tmin1=tmin1, tmin2=tmin2, 
                                 tname=tname, x=x,data=boot.data,
                                 dropalwaystreated=dropalwaystreated,
                                 idname="id", uniqueid=uniqueid, zname=zname,
                                 y.seq=y.seq, dy.seq=dy.seq, probs=probs, h=h,
                                 probevals=probevals, bootstrap.iter=TRUE)
    eachIter[[i]] = list(att = thisIter$att, qte=thisIter$qte)
  }
  
  #For now, just plot the qte and att with standard errors
  #90% confidence interval
  #att
  #helper function to get the first element out of a list
  getElement <- function(Lst, elemNum) {
    return(as.numeric(unlist((Lst[elemNum])))) #as.numeric is a trick to 
    #get numerical value of qte
  }
  all.att = unlist(sapply(eachIter, FUN=getElement,elemNum=1))
  #reorder asc
  all.att = all.att[order(all.att)]
  att.upper = all.att[min(iters,round((1-alp/2)*iters))]
  att.lower = all.att[max(1,round((alp/2)*iters))]
  #set up a blank plot
  plot(1,type="n",ylim=c(-2,2)) #these numbers come from having a good idea
  #what the result is, and then subject to some modifications
  abline(a=fy3$att, b=0)
  abline(a=att.upper, b=0, lty=2)
  abline(a=att.lower, b=0, lty=2)
  
  #now get CI for qte:
  all.qte = lapply(eachIter, FUN=getElement, elemNum=2)
  qte.mat = do.call(rbind,lapply(all.qte, FUN=as.numeric, ncol=length(all.qte[[1]]), byrow=TRUE))
  #order each column
  sorted.qtemat = apply(qte.mat, 2, sort)
  qte.upper = sorted.qtemat[round((1-alp/2)*iters),]
  qte.lower = sorted.qtemat[max(1,round((alp/2)*iters)),]
  
  plot(seq(0,1,length.out=ncol(qte.mat)), fy3$qte, type="l", ylim=c(-2,1))
  lines(seq(0,1,length.out=ncol(qte.mat)), qte.upper, lty=2)
  lines(seq(0,1,length.out=ncol(qte.mat)), qte.lower, lty=2)
  abline(a=0,b=0)  
  
  return(list(qte=fy3$qte, qte.upper=qte.upper, qte.lower=qte.lower,
         att=fy3$att, att.upper=att.upper, att.lower=att.lower))
}


####FIRPO########
#'@title firpo
#' @inheritParams threeperiod.fanyu
#' @export
firpo = function(formla, x=NULL, data, 
                 y.seq=NULL, dy.seq=NULL, probs=seq(0,1,0.1)) {
  form = as.formula(formla)
  dta = model.frame(terms(form,data=data),data=data) #or model.matrix
  colnames(dta) = c("y","treatment")
  yname="y"
  treat="treatment"
  data=cbind.data.frame(dta,data)
  
  #setup the data
  treated = data[data[,treat]==1,]
  untreated = data[data[,treat]==0,]
  
  n = nrow(data)
  #estimate the propensity score
  pscore=fitted(glm(data[,treat] ~ as.matrix(data[,x]),
                    family=binomial(link="probit")))
  p = rep(nrow(treated)/(nrow(treated) + nrow(untreated)), n)
  #there are alternatives for how to compute the quantiles of 
  #treated outcomes for the treated group:
  #1) compute quantiles directly
  treated.quantiles = quantile(treated[,yname], probs=probs)
  #2) use firpo method
  #checkfun will be called by the various functions to be minimized
  #in this routine
  checkfun = function(a, tau) {
    return(a*(tau - (1*(a<=0))))
  }
  treated.weights = data[,treat] / sum(data[,treat])
  minfun.inner = function(q, tau, weights) {
    retval = sum(weights*checkfun(data[,yname]-q, tau))
    return(retval)
  }

  get.firpo.quantiles = function(tau, weights) {
    return(optimize(minfun.inner, 
                    lower=min(data[,yname]),
                    upper=max(data[,yname]),
                    tau=tau,weights=weights)$minimum)
  }
  treated.firpo.quantiles = vapply(probs, FUN=get.firpo.quantiles,
                                   FUN.VALUE=1, weights=treated.weights)
  
  untreated.weights = (pscore/(1-pscore))*((1-data[,treat])/sum(data[,treat]))
  untreated.firpo.quantiles = vapply(probs, FUN=get.firpo.quantiles,
                                     FUN.VALUE=1, weights=untreated.weights)
  
  firpo.qte = treated.firpo.quantiles - untreated.firpo.quantiles

  #just for completeness, calculate the distribution of each potential outcomes
  #using moment conditions.
  F.treated <- ecdf(treated[,yname])
  F.treatedcf.fun <- function(y) {
      pterm <- pscore/((1-pscore)*p)
      Dterm <- 1 - data[,treat]
      yterm <- 1*(data[,yname] < y)
      mean(pterm*Dterm*yterm)
  } #something appears to be off here for 0 wages, otherwise, everything good!
  y.seq <- seq(min(data[,yname]), max(data[,yname]), length.out=500)
  F.treatedcf = approxfun(y.seq,
      vapply(y.seq, FUN=F.treatedcf.fun, FUN.VALUE=1)
      , method="constant", yleft=0, yright=1, f=0, ties="ordered")
  class(F.treatedcf) = c("ecdf", "stepfun", class(F.treatedcf.fun))
  assign("nobs", nrow(treated), envir = environment(F.treatedcf))
  
  
  return(list(qte=firpo.qte, F.treated=F.treated, F.treatedcf=F.treatedcf,
              pscore=pscore, p=p, F.treatedcf.fun=F.treatedcf.fun))
  
}


######GENERAL HELPER FUNCTIONS#######

##change function takes a panel data set, and takes the time 
##difference of a particular variable, returning NA
##in the first year
change <- function(xname, idname, tname, data, withBetweenYears=TRUE) {
  mint = min(data[,tname])
  maxt = max(data[,tname])
  data$outvec = NA #we will change this value, but only
  #for the years 2-end
   
  #debug
  out.xname <<-xname
  out.data <<- data
  out.withBetweenYears <<- withBetweenYears
  out.idname <<- idname
  out.tname <<-tname
  
  if (withBetweenYears) {
    for(t in (mint+1):maxt) {
      for(id in unique(data[,idname])) {
        data[which(t == data[,tname] & id == data[,idname]),]$outvec =
          data[which(t == data[,tname] & id == data[,idname]),xname] -
          data[which((t-1) == data[,tname] & id == data[,idname]),xname]
      }
    }
  } else {
    #the else case is if you only want to get the lag between two years
    for (id in unique(data[,idname])) {
      ## Debug
      out.id <<- id
      out.data <<- data
      out.outvec <<- out.data$outvec
      out.breakmatch <<- data[which(mint == data[,tname] & id == data[,idname]),xname]
      ##
      data[which(maxt == data[,tname] & id == data[,idname]),]$outvec =
        data[which(maxt == data[,tname] & id == data[,idname]),xname] -
        data[which(mint == data[,tname] & id == data[,idname]),xname]
    }
  }
  return(data$outvec)
} 
  
#getLag is the exact same as change function
#except just recover the lagged value instead
#of taking difference with something you already
#have.
#Changes: add parameter withBetweenYears
# when this is set to false, the function will only
# compute the lag over the two time periods in the data
# [probably should just adjust this so that it can 
# work with whatever data is passed, but this will
# work for now.]
getLag <- function(xname, idname, tname, data, withBetweenYears=TRUE) {
  mint = min(data[,tname])
  maxt = max(data[,tname])
  data$outvec = NA #we will change this value, but only
  #for the years 2-end
  
  #debug
  out.xname <<-xname
  out.data <<- data
  out.withBetweenYears <<- withBetweenYears
  out.idname <<- idname
  out.tname <<-tname
  
  if (withBetweenYears) {
    for(t in (mint+1):maxt) {
      for(id in unique(data[,idname])) {
        data[which(t == data[,tname] & id == data[,idname]),]$outvec =
          data[which((t-1) == data[,tname] & id == data[,idname]),xname]
      }
    }
  } else {
    #the else case is if you only want to get the lag between two years
    for (id in unique(data[,idname])) {
      ## Debug
      out.id <<- id
      out.data <<- data
      out.outvec <<- out.data$outvec
      out.breakmatch <<- data[which(mint == data[,tname] & id == data[,idname]),xname]
      ##
      data[which(maxt == data[,tname] & id == data[,idname]),]$outvec =
        data[which(mint == data[,tname] & id == data[,idname]),xname]
      #test.data <<- data ...Not sure what this does
    }
  }
  return(data$outvec)
}

###makeBalancedPanel is a function to take a dataset
# and make sure that all years are available for 
# all observations.  If some years are not available,
# then that observation is dropped.
#'@title makeBalancedPanel
#' @description This function drops observations from data.frame
#' that are not part of balanced panel data set.
#' @param data data.frame used in function
#' @param idname unique id
#' @param tname time period name
#' @export
makeBalancedPanel <- function(data, idname, tname) {
  data=droplevels(data)
  allt = unique(data[,tname])
  allid = unique(data[,idname])
  
  #loop over each id in the dataset
  for (id in allid) {
    #get the number of time periods for that id
    this.allt = unique(data[data[,idname]==id,tname])
    
    #check if its equal to the largest number of time
    #periods in the dataset
    if (!(length(this.allt) == length(allt))) {
      #if it is fewer, then drop all observations
      #from that id from the dataset
      data = data[!(data[,idname] == id),]
    }
  }
  return(data)
}


#pass in two sets of observations and the support of each
#getJoint returns the joint distribution as a matrix
getJoint <- function(x1, y1, xrange=NULL, yrange=NULL) {
  #num = length(x) #this works as long as x and y have the same length
  #outmat = matrix(0, nrow=length(yrange), ncol=length(xrange))
  #rownames(outmat) = yrange
  #colnames(outmat) = xrange
  #for (i in xrange) {
  #  for (j in yrange) {
  #    outmat[which(yrange==j),which(xrange==i)] = (sum(1*(x <= i & y<=j)))/num
  #  }
  #}
  #return(list(outmat=outmat, xrange=xrange, yrange=yrange))
  #changey = as.numeric((treated.tmin1[,yname] - treated.tmin2[,yname]))
  #tmin2y=as.numeric(treated.tmin2[,yname])
  #F.joint.tmin1 = npudist(~changey+tmin2y,ckernel="uniform")
  
  #this works, but is slow:
  #joint = npudist(~x1+y1, ckernel="uniform")
  return(function(u,v){
    #this is really slow
    #predict(joint, newdata=data.frame(x1=u,y1=v))
    #out = rep(0,length(v))
    #for (i in 1:length(v)) {
    #  out[i] = sum(1*(x1<u & y1<v[i]))/length(y1)
    #}
    #for empirical CDF, C++ code is faster
    return(getJointUVCPP(u,v,x1,y1))
    #return(sum(1*(x1<u & y1<v))/length(y1))
  })
}


#getDistVal takes in a distribution as a List and
# a particular point and returns the value
#of the distribution at that particular point
getDistVal <- function(distObj, point) {
  disty = distObj[[1]]
  xvals = distObj[[2]]
  retval = disty[which(xvals >= point)[1]]
  return(retval)
}

#getJointDistVal takes in a joint distribution
#and two points and returns the value of the
#joint distribution at those two points
getJointDistVal <- function(jointDistObj, pointx, pointy) {
  disty = jointDistObj[[1]]
  xvals = jointDistObj$xrange
  yvals = jointDistObj$yrange
  retval = disty[ which(yvals >= pointy)[1], 
                  which(xvals >= pointx)[1]]
  return(retval)
}

#getCondDist takes in a joint distribution
#and returns a conditional distribution for a
#particular value of the 2nd r.v. in the joint
#distribution
#****This doesn't work correctly
getCondDist <- function(jointDistObj, pointy) {
  disty = jointDistObj[[1]]
  yvals =  jointDistObj$yrange
  retDisty = disty[which(yvals >= pointy)[1],]
  return(list(dist=retDisty, xrange=jointDistObj$xrange, yval=pointy))
}

##getSumDist takes in a joint distribution G(x,y) and a point z
#and calculates the probability P(X+Y<z)
getSumDist <- function(jointPdfObj, z) {
  jointpdf = jointPdfObj[[1]]
  xvals = jointPdfObj$xrange
  yvals = jointPdfObj$yrange
  sumout = 0
  for (y in yvals) {
    for (x in xvals[xvals<(z-y)]) {
      sumout = sumout + jointpdf[which(yvals==y), which(xvals==x)]
    }
  }
  return(sumout)
}

##jointCDF.to.jointPDF takes a discrete CDF and turns
# it into a pdf
jointCDF.to.jointPDF <- function(jointDistObj) {
  disty = jointDistObj[[1]]
  xvals = jointDistObj$xrange
  yvals = jointDistObj$yrange
  jointpdf = matrix(0,nrow=(length(yvals)-1),
                    ncol=(length(xvals)-1))
  colnames(jointpdf) = xvals[2:length(xvals)]
  rownames(jointpdf) = yvals[2:length(yvals)]
  for (x in xvals[2:length(xvals)]) {
    for (y in yvals[2:length(yvals)]) {
      jointpdf[which(yvals==y)[1]-1, which(xvals==x)[1]-1] =
        #old way, subtract diagonal (above and left).
        #I think it overcounts
        #disty[which(yvals==y)[1], which(xvals==x)[1]] -
        #disty[which(yvals==y)[1]-1, which(xvals==x)[1]-1]
        
        #2nd try, subtract one above, and one to the left
        #individually, and take diff in diff, I think will work
        ((disty[which(yvals==y)[1], which(xvals==x)[1]] -
              disty[which(yvals==y)[1], which(xvals==x)[1]-1]) -
        (disty[which(yvals==y)[1]-1, which(xvals==x)[1]] -
           disty[which(yvals==y)[1]-1, which(xvals==x)[1]-1]))
    }
  }
  return(list(jointpdf=jointpdf, xrange=xvals[2:length(xvals)], 
              yrange=yvals[2:length(yvals)]))
}

####HELPER FUNCTIONS FOR ATHEY-IMBENS####

#bootstrap.athey.imbens is a function that computes bootstrap
#standard errors for quantile treatment effects
#'@title bootstrap.athey.imbens
#' @inheritParams threeperiod.fanyu
#' @export
bootstrap.athey.imbens <- function(formla, t, tmin1, tname, x=NULL,data,
                                  dropalwaystreated=TRUE, idname, 
                                   uniqueid, iters=20) {
  
  ##This is pasted from Athey-Imbens method!
  #drop the always treated.  Note that this also relies
  #on no "switchback" or no return to untreated status
  #after joining treatment.
  if (dropalwaystreated) {
    #first line gets the correct two years of data
    data = subset(data, (data[,tname]==tmin1 | data[,tname]==t))
    data = makeBalancedPanel(data, idname, tname)
    #now, drop the observations that are initially treated
    #data$lagged.treat = getLag(xname=treat, idname=idname,
    #                           tname=tname, data=data, 
    #                           withBetweenYears=FALSE)
    #data = subset(data, !(data[,treat]==1 & data[,tname]==tmin1))
    #data = subset(data, !(data[,lagged.treat]==1 & data[,tname]==t))
    #below will be application specific, but I NEED SOMETHING
    # to work.
    data = subset(data, !(d03==1 & d10==1))
  }
  
  #just to make sure the factors are working ok
  data = droplevels(data)
  ##
  
  #first calculate the actual estimate
  ai = athey.imbens(formla, t, tmin1, tname, x, data,
                    dropalwaystreated, idname, uniqueid)
  #now calculate the bootstrap confidence interval
  eachIter = list()
  #Need to build dataset by sampling individuals, and then
  #taking all of their time periods
  all.ids = unique(data[,idname])
  n = length(all.ids)
  out.bootdatalist <<- list()
  for (i in 1:iters) {
    #reset boot.data
    out.iter <<- i
    boot.data = data[0,]
    randy = sample(1:n, n, replace=T)
    #there has to be a way to do this faster, but go with the loop
    #for now
    for (j in all.ids[randy]) {
      boot.data = rbind(boot.data, data[(data[,idname]==j),])
    }
    boot.data = process.bootdata(boot.data, idname, uniqueid)
    out.bootdatalist[[i]] <<- boot.data
    thisIter = athey.imbens(formla, t, tmin1, tname, x, boot.data,
                            dropalwaystreated, idname, uniqueid)
    eachIter[[i]] = list(att = thisIter$att, qte=thisIter$qte)
  }
  
  #For now, just plot the qte and att with standard errors
  #90% confidence interval
  #att
  #helper function to get the first element out of a list
  getElement <- function(Lst, elemNum) {
    return(as.numeric(unlist((Lst[elemNum])))) #as.numeric is a trick to 
     #get numerical value of qte
  }
  all.att = unlist(sapply(eachIter, FUN=getElement,elemNum=1))
  #reorder asc
  all.att = all.att[order(all.att)]
  att.upper = all.att[round(0.95*iters)]
  att.lower = all.att[round(0.05*iters)]
  #set up a blank plot
  plot(1,type="n",ylim=c(0.3,0.55)) #these numbers come from having a good idea
   #what the result is, and then subject to some modifications
  abline(a=ai$att, b=0)
  abline(a=att.upper, b=0, lty=2)
  abline(a=att.lower, b=0, lty=2)
  
  #now get CI for qte:
  all.qte = lapply(eachIter, FUN=getElement, elemNum=2)
  qte.mat = do.call(rbind,lapply(all.qte, FUN=as.numeric, ncol=length(all.qte[[1]]), byrow=TRUE))
  #order each column
  sorted.qtemat = apply(qte.mat, 2, sort)
  qte.upper = sorted.qtemat[round(0.95*iters),]
  qte.lower = sorted.qtemat[round(0.05*iters),]
  
  plot(seq(0,1,length.out=ncol(qte.mat)), ai$qte, type="l", ylim=c(-1,3))
  lines(seq(0,1,length.out=ncol(qte.mat)), qte.upper, lty=2)
  lines(seq(0,1,length.out=ncol(qte.mat)), qte.lower, lty=2)
  abline(a=0,b=0)
  
}

#this function adjusts bootdata idnames
#so that there are no more duplicates
process.bootdata <- function(data, idname, uniqueid) {
  #duplicated function works perfectly for this, I think
  this.data = data
  add.to.end = 1
  while (any(duplicated(this.data[,uniqueid]))) {
   this.data[duplicated(this.data[,uniqueid]), idname] = 
     sapply(this.data[duplicated(this.data[,uniqueid]), idname],
            FUN=paste, add.to.end, sep="-")
   #this doesn't follow exact naming pattern as previous
   #unique ids, but I think it will be ok because we don't
   #use it anywhere but here; and it is easy to implement
   this.data[duplicated(this.data[,uniqueid]), uniqueid] = 
     sapply(this.data[duplicated(this.data[,uniqueid]), uniqueid],
            FUN=paste, add.to.end, sep="-")
   add.to.end = add.to.end + 1
   if (add.to.end > 20) {
     #do this just in case we start spinning
     
   }
  }
  return(this.data)
}



#####HELPER FUNCTIONS FOR FAN-YU######

#getlb calculates the lower bound
#s a particular value of distribution for 
# which to calculate the bound; usually
# will call this function with vapply
# so that it computes lots of bounds all at once
# ecdf.change.newlytreated an ecdf object (so that)
# we can use it as a function of the distribution
# of the change in outcomes for the newly treated
# ecdf.newlytreated.tmin1 is the distribution
# of outcomes for the newly treated in the previous
# period (when they were part of the no treatment regime)
# What is y?
getlb <- function(s, ecdf.change.newlytreated, ecdf.newlytreated.tmin1, y) {
  #OLD - this worked before iv stuff
  #return(max(F.changey.treated(y) + F.Y0tmin1.Dt1(s-y)-1,0))
  
  #new
  return(max(ecdf.change.newlytreated(y) + ecdf.newlytreated.tmin1(s-y) - 1,0))
}

#Get the Fan-Yu upper bound
#same as the lower bound, but with slightly different formula
getub <- function(s, ecdf.change.newlytreated, ecdf.newlytreated.tmin1, y) {
  #OLD - this worked before iv stuff
  #return(max(F.changey.treated(y) + F.Y0tmin1.Dt1(s-y)-1,0))
  
  #new
  return(1 + min(ecdf.change.newlytreated(y) + ecdf.newlytreated.tmin1(s-y) - 1,0))
}

#this is probably simpler than what R's quantile function does
#but I think that it works; otherwise, couldn't figure out
#how to invert cdf that I was generating.
simple.quantile <- function(x,Fx,probs=c(0,0.25,0.5,0.75,1)) {
  inner.quantfun <- function(prob,x,Fx) {
    ind = which.max(Fx>=prob) #this works by picking out smallest (because
    #of the way ties are handled) F(x) satisfying condition
    return(x[ind])
  }
  out = vapply(probs,FUN=inner.quantfun,FUN.VALUE=1,x=x,Fx=Fx)
  names(out) = probs
  return(out)
}
