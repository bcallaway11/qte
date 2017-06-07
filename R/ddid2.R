utils::globalVariables(c("yname", "treat", "panel", "treated.t", "treated.tmin1", "untreated.t", "untreated.tmin1", "F.treated.t", "F.treated.tmin1", "F.untreated.t", "F.untreated.tmin1", "xformla", "data", "x", "wname", "probs", "method", "eachIter"))


#####Two Period DDID#####
##Idea here is that we can use information from a third period
##to point identify counterfactual distribution of outcomes
##for the treated group
##call plot function, summary function, formula function, etc. later
##add functionality to pass in pscore
#' @title compute.ddid2
#'
#' @description
#' \code{compute.ddid2} uses two periods of data (repeated cross sections
#' or panel) to estimate the Quantile Treatment Effect on the Treated (QTET)
#'
#' @param qp QTEparams object containing the parameters passed to ddid2
#'
#' @importFrom quantreg rq
#' @import BMisc
#'
#' @return QTE object
#'
#' @keywords internal
compute.ddid2 <- function(qp) {

    setupData(qp)
    bootstrapiter <- qp$bootstrapiter

    ##b)
    if (panel) {
        untreated.change.t <- untreated.t[,yname] - untreated.tmin1[,yname]

    } else {
        untreated.change.t <- cs2panel(untreated.t, untreated.tmin1, yname)
    }

    F.untreated.change.t <- stats::ecdf(untreated.change.t)
    ##calculate the distribution of the change for the treated group;
    ## this will be changed if there are covariates
    F.treated.change.t <- F.untreated.change.t
    

    ##calculate the att; this will be changed if there are covariates
    att = mean(treated.t[,yname]) - mean(treated.tmin1[,yname]) -
        (mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname]))

    ##now compute the average over the treated observations
    ##quantys1 <- quantile(F.treated.change.t,
    ##                     probs=F.untreated.change.t(untreated.t[,yname] -
    ##                         untreated.tmin1[,yname]), type=1)


    quantys1 <- untreated.change.t

    quantys2 <- stats::quantile(F.treated.tmin1,
                         probs=F.untreated.tmin1(untreated.tmin1[,yname]),
                         type=1)

    y.seq <- unique((quantys1+quantys2)[order(quantys1 + quantys2)])

    F.treated.t.cf.val <- vapply(y.seq,
                                 FUN=function(y) { mean(1*(quantys1 + quantys2 <=
                                                           y)) }, FUN.VALUE=1)

    F.treated.t.cf <- makeDist(y.seq, F.treated.t.cf.val)

    pscore.reg <- NULL
    if (!(is.null(x))) {

        stop("method not implemented with covariates")

        ## this.formla <- y ~ x ##just set up dummy formula first
        ## stats::lhs(this.formla) <- as.name(yname)
        ## stats::rhs(this.formla) <- stats::rhs(xformla)

        ## tauu <- seq(0,1, length.out=nrow(untreated.tmin1))

        ## condQ.untreated.tmin1.qr <- rq(this.formla, tau=tauu,
        ##                                data=untreated.tmin1)
        
        ## condF.untreated.tmin1 <- stats::predict(condQ.untreated.tmin1.qr,
        ##                                  newdata=untreated.tmin1,
        ##                                  stepfun=T, type="Fhat")

        ## qs <- c()
        ## for (i in 1:nrow(untreated.tmin1)) {
        ##     qs[i] <- condF.untreated.tmin1[[i]](untreated.tmin1[,yname][i])
        ## }

        ## tauu1 <- seq(0,1,length.out=nrow(untreated.t))
        ## condQ.untreated.t.qr <- rq(this.formla, tau=tauu1, data=untreated.t)

        ## condQ.untreated.t <- stats::predict(condQ.untreated.t.qr,
        ##                                newdata=untreated.tmin1,
        ##                                stepfun=T, type="Qhat")
        ## quantys1 <- c()
        ## for (i in 1:nrow(untreated.tmin1)) {
        ##     quantys1[i] <- condQ.untreated.t[[i]](qs[i])
        ## }

        ## taut <- seq(0,1,length.out=nrow(treated.tmin1))
        ## condQ.treated.tmin1.qr <- rq(this.formla, tau=taut, data=treated.tmin1)

        ## condQ.treated.tmin1 <- stats::predict(condQ.treated.tmin1.qr,
        ##                                newdata=untreated.tmin1,
        ##                                stepfun=T, type="Qhat")
        ## quantys2 <- c()
        ## for (i in 1:nrow(untreated.tmin1)) {
        ##     quantys2[i] <- condQ.treated.tmin1[[i]](qs[i])
        ## }

        ## innerdf <- untreated.tmin1
        ## innerdf$y <- quantys1+quantys2-untreated.tmin1[,yname]
                              
        ## condQ.treated.t.cf.qr <- rq(this.formla, tau=tauu, data=innerdf)

        ## taut1 <- seq(0,1,length.out=nrow(treated.t))
        ## condQ.treated.t.qr <- rq(this.formla, tau=taut1, data=treated.t)



        ## return(list(condQ.treated.t.qr=condQ.treated.t.qr,
        ##             condQ.treated.t.cf.qr=condQ.treated.t.cf.qr))
                                                 
                                          


        ## ## temp <- dr.predict(untreated.tmin1[,yname], data.frame(sex="Female"), condF.untreated.tmin1)

        ## ## temp2 <- dr.predict(untreated.tmin1[,yname], data.frame(sex="Male"), condF.untreated.tmin1)


        ## ## graphics::plot(makeDist(untreated.tmin1[,yname], temp))
        ## ## plot(predict(condQ.untreated.tmin1, newdata=data.frame(sex="Female"),
        ## ##              stepfun=T, type="Fhat"), add=T, col="blue")

        
        ## ## plot(makeDist(untreated.tmin1[,yname], temp2), col="blue", add=T)


        








        
        
        ##  ##Step 2: Distribution regression / quantile regression        
        ## xmat0 <- as.matrix(untreated.tmin1[,x]) #drop intercept
        ## xmat1 <- as.matrix(treated.tmin1[,x]) #drop intercept

        ## yvals <- unique(stats::quantile(F.untreated.tmin1, probs=seq(.01,.99,.01)))
        ## dr.list <- list()
        ## formly1 <- "I(1*(y <= yvals[i])) ~ "
        ## formly2 <- paste(x[-1], collapse=" + ")
        ## formly <- paste(formly1, formly2, sep="")
        ## formly <- as.formula(formly)
        ## for (i in 1:length(yvals)) {
        ##     dr.list[[i]] <- glm(formly, data=untreated.tmin1, family=binomial(link=logit))
        ## }

        
       
        ## tauvals <- seq(.01,.99,.01)
        ## formlq1 <- "y ~"
        ## formlq2 <- paste(x[-1], collapse=" + ")
        ## formlq <- paste(formlq1, formlq2, sep="")
        ## formlq <- as.formula(formlq)
        ## qr <- rq(formlq, data=treated.tmin1, tau=tauvals)
        ## round.ytmin1 <- round(untreated.tmin1[,yname],2)

        ## ##function to take in tau and x and return F^-1(tau|x)
        ## quant.treated.tmin1.x <- function(tau, x) {
        ##     tau.idx <- which.min(abs(tauvals-tau))
        ##     return(predict(qr, newdata=x)[tau.idx]) ##this step can probably be sped up because can predict all x's at the same time
        ##     ##also might be good to replace this with distribution regression
        ##     ## and then inverting because of the 0s.
        ## }

        ## ##function to get everything for single values of delta, y_t-1, and x
        ## itfun <- function(dta,y) {
        ##     waits <- dta["waits"]
        ##     change <- dta["changey"]
        ##     ytmin1 <- dta["ytmin1"]
        ##     xx <- data.frame(t(dta[x[-1]]))
        ##     waits*1*(change + quant.treated.tmin1.x(F.untreated.tmin1.x(ytmin1, xx), xx) <= y)
        ## }

        ## itfun.y <- function(y) {
        ##     mean(apply(untreated.t, 1, itfun, y)) / mean(untreated.t$waits) ##the last division is to scale for estimated propensity score -- asymptotically it is equal to 1.
        ## }

        ## temp <- vapply(yvals, itfun.y, 1.0)

        ## F.treated.t.cf = stats::approxfun(yvals,
        ##     temp, method="constant",
        ##     yleft=0, yright=1, f=0, ties="ordered")
        ## class(F.treated.t.cf) = c("ecdf", "stepfun",
        ##          class(F.treated.t.cf))
        ## assign("nobs", length(yvals), envir = environment(F.treated.t.cf))

        ## ##compute the att using abadie-2005
        ## ##TODO: this will only work for panel
        ## att <- mean(((pscore.data$changey)/pD1)*(pscore.data[,treat] - pscore) /
        ##             (1-pscore))
        
    }





    ##compare this to the actual outcomes
    F.treated.t <- stats::ecdf(treated.t[,yname])

    qte <- stats::quantile(F.treated.t, probs=probs) -
        stats::quantile(F.treated.t.cf, probs=probs)

    out <- QTE(F.treated.t=F.treated.t,
               F.treated.tmin1=F.treated.tmin1,
               F.untreated.change.t=F.untreated.change.t,
               F.untreated.t=F.untreated.t,
               F.untreated.tmin1=F.untreated.tmin1,
               F.treated.t.cf=F.treated.t.cf,
               qte=qte, pscore.reg=pscore.reg,  ate=att, probs=probs)
    class(out) <- "QTE"
    return(out)
}




#' @title ddid2
#'
#' @description \code{ddid2} computes the Quantile Treatment Effect
#' on the Treated (QTET) using the method of Callaway, Li, and Oka (2015).  
#'
#' @inheritParams ci.qte
#' @param formla The formula y ~ d where y is the outcome and d is the
#'  treatment indicator (d should be binary)
#' @param xformla A optional one sided formula for additional covariates that
#'  will be adjusted for.  E.g ~ age + education.  Additional covariates can
#'  also be passed by name using the x paramater.
#' @param t The 3rd time period in the sample (this is the name of the column)
#' @param tmin1 The 2nd time period in the sample (this is the name of the
#'  column)
#' @param tname The name of the column containing the time periods
#' @param data The name of the data.frame that contains the data
#' @param panel Boolean indicating whether the data is panel or repeated cross
#'  sections
#' @param dropalwaystreated How to handle always treated observations
#'  in panel data case (not currently used)
#' @param idname The individual (cross-sectional unit) id name
#' @param probs A vector of values between 0 and 1 to compute the QTET at
#' @param iters The number of iterations to compute bootstrap standard errors.
#'  This is only used if se=TRUE
#' @param alp The significance level used for constructing bootstrap
#'  confidence intervals
#' @param method The method for estimating the propensity score when covariates
#'  are included
#' @param se Boolean whether or not to compute standard errors
#' @param retEachIter Boolean whether or not to return list of results
#'  from each iteration of the bootstrap procedure
#' @param seedvec Optional value to set random seed; can possibly be used
#'  in conjunction with bootstrapping standard errors.
#'
#' @examples
#' ##load the data
#' data(lalonde)
#'
#' ## Run the panel.qtet method on the experimental data with no covariates
#' pq1 <- ddid2(re ~ treat, t=1978, tmin1=1975, tname="year",
#'  x=NULL, data=lalonde.psid.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(pq1)
#'
#' ## Run the panel.qtet method on the observational data with no covariates
#' 
#'
#' @references
#' Callaway, Brantly, Tong Li, and Tatsushi Oka.  ``Quantile Treatment Effects
#' in Difference in Differences Models under Dependence Restrictions and with
#' Only Two Time Periods.'' Working Paper, 2015.
#'
#' @return \code{QTE} object
#' 
#' @export
ddid2 <- function(formla, xformla=NULL, t, tmin1,
                      tname, data, panel=FALSE,
                      dropalwaystreated=TRUE, idname=NULL, probs=seq(0.05,0.95,0.05),
                      iters=100, alp=0.05, method="logit", se=TRUE,
                      retEachIter=FALSE, seedvec=NULL, pl=FALSE, cores=NULL) {


    qp <- QTEparams(formla=formla, xformla=xformla, t=t, tmin1=tmin1,
                    tname=tname, data=data, panel=panel,
                    idname=idname, probs=probs,
                    iters=iters, alp=alp, method=method,
                    se=se, retEachIter=retEachIter, seedvec=seedvec,
                    pl=pl, cores=cores)
                    
    
    
    ## maybe, move this to setupData
    ##treated.t = data[data[,tname]==t & data[,treat]==1,]
    ##treated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 1, ]
    ##untreated.t = data[data[,tname]==t & data[,treat]==0,]
    ##untreated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 0, ]

    ##first calculate the actual estimate
    pqte = compute.ddid2(qp)

    if (se) {

        qp$bootstrapiter <- TRUE

        ##bootstrap the standard errors
        SEobj <- bootstrap(qp, pqte, compute.ddid2)



        ##could return each bootstrap iteration w/ eachIter
        ##but not currently doing that
        out <- QTE(qte=pqte$qte, qte.upper=SEobj$qte.upper,
                    qte.lower=SEobj$qte.lower, ate=pqte$ate,
                    ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                    qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                    F.treated.t=pqte$F.treated.t,
                    F.untreated.t=pqte$F.untreated.t,
                    F.treated.t.cf=pqte$F.treated.t.cf,
                    F.treated.tmin1=pqte$F.treated.tmin1,
                    F.treated.tmin2=pqte$F.treated.tmin2,
                    F.treated.change.tmin1=pqte$F.treated.change.tmin1,
                    F.untreated.change.t=pqte$F.untreated.change.t,
                    F.untreated.change.tmin1=pqte$F.untreated.change.tmin1,
                    F.untreated.tmin1=pqte$F.untreated.tmin1,
                    F.untreated.tmin2=pqte$F.untreated.tmin2,
                    pscore.reg=pqte$pscore.reg,
                    eachIterList=eachIter,
                    probs=probs)
        return(out)
    } else {
        return(pqte)
    }
}
