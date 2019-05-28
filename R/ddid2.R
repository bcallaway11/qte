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
#' @export
compute.ddid2 <- function(qp) {

    setupData(qp)
    bootstrapiter <- qp$bootstrapiter

    ##calculate the att; this will be changed if there are covariates
    att = mean(treated.t[,yname]) - mean(treated.tmin1[,yname]) -
        (mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname]))

    ##now compute the average over the treated observations
    ##quantys1 <- quantile(F.treated.change.t,
    ##                     probs=F.untreated.change.t(untreated.t[,yname] -
    ##                         untreated.tmin1[,yname]), type=1)


    quantys1 <- untreated.change.t

    ## will update this term if there are covariates
    quantys2 <- stats::quantile(F.treated.tmin1,
                         probs=F.untreated.tmin1(untreated.tmin1[,yname]),
                         type=1)


    pscore.reg <- NULL
    qr0.reg <- NULL
    u <- seq(.01,.99,.01)
    if (!(is.null(x))) {

        n0 <- nrow(untreated.t)
        n1 <- nrow(treated.t)

        yformla <- toformula("y",BMisc::rhs.vars(xformla))
        QR1tmin1 <- rq(yformla, data=treated.tmin1, tau=u)
        QR0tmin1 <- rq(yformla, data=untreated.tmin1, tau=u)
        rank0tmin1 <- predict(QR0tmin1, type="Fhat", stepfun=TRUE) ## rank untreated at tmin1 (inner step)
        rank0tmin1 <- sapply(1:n0, function(i) rank0tmin1[[i]](untreated.tmin1$y[i]))
        ytmin1 <- predict(QR1tmin1, newdata=untreated.tmin1, type="Qhat", stepfun=TRUE) ## transformation from callaway-li-oka
        ytmin1 <- sapply(1:n0, function(i) ytmin1[[i]](rank0tmin1[i])) ## actual value
        quantys2 <- ytmin1

        ## uncomment this if want qr results for untreated potential outcomes
        ## untreated.tmin1$ddy <- quantys1 + quantys2
        ## ddyformla <- toformula("ddy",BMisc::rhs.vars(xformla))
        ## qr0.reg <- rq(ddyformla, data=untreated.tmin1, probs=tau)
    }

    ## build counterfactual distribution
    y.seq <- unique((quantys1+quantys2)[order(quantys1 + quantys2)])

    F.treated.t.cf.val <- vapply(y.seq,
                                 FUN=function(y) { mean(1*(quantys1 + quantys2 <=
                                                           y)) }, FUN.VALUE=1)

    F.treated.t.cf <- makeDist(y.seq, F.treated.t.cf.val)

    ##compare this to the actual outcomes
    F.treated.t <- stats::ecdf(treated.t[,yname])

    qte <- stats::quantile(F.treated.t, probs=probs) -
        stats::quantile(F.treated.t.cf, probs=probs)

    if (!is.null(x)) {
        att <- mean(treated.t$y) - sum(quantile( F.treated.t.cf, probs=u, type=1 ))/length(u)
    }
    
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
#' ## Run the ddid2 method on the observational data with no covariates
#' d1 <- ddid2(re ~ treat, t=1978, tmin1=1975, tname="year",
#'  data=lalonde.psid.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(d1)
#'
#' ## Run the ddid2 method on the observational data with covariates
#' d2 <- ddid2(re ~ treat, t=1978, tmin1=1975, tname="year",
#'  data=lalonde.psid.panel, idname="id", se=FALSE,
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(d2)
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
                      tname, data, panel=TRUE,
                      dropalwaystreated=TRUE, idname=NULL, probs=seq(0.05,0.95,0.05),
                      iters=100, alp=0.05, method="logit", se=TRUE,
                      retEachIter=FALSE, seedvec=NULL, pl=FALSE, cores=NULL) {

    if (!panel) {
        stop("method not implemented with repeated cross sections data...\n  In this case, try change in changes method...")
    }

    data <- panelize.data(data, idname, tname, t, tmin1)   

    
    qp <- QTEparams(formla=formla, xformla=xformla, t=t, tmin1=tmin1,
                    tname=tname, data=data, panel=panel,
                    idname=idname, probs=probs,
                    iters=iters, alp=alp, method=method,
                    se=se, retEachIter=retEachIter, seedvec=seedvec,
                    pl=pl, cores=cores)
                    

    panel.checks(qp)
    
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
                   c=SEobj$c,
                   F.treated.t=pqte$F.treated.t,
                   F.untreated.t=pqte$F.untreated.t,
                   F.treated.t.cf=pqte$F.treated.t.cf,
                   F.treated.tmin1=pqte$F.treated.tmin1,
                   ##F.treated.tmin2=pqte$F.treated.tmin2,
                   ##F.treated.change.tmin1=pqte$F.treated.change.tmin1,
                   F.untreated.change.t=pqte$F.untreated.change.t,
                   ##F.untreated.change.tmin1=pqte$F.untreated.change.tmin1,
                   F.untreated.tmin1=pqte$F.untreated.tmin1,
                   ##F.untreated.tmin2=pqte$F.untreated.tmin2,
                   pscore.reg=pqte$pscore.reg,
                   eachIterList=eachIter,
                   probs=probs)
        return(out)
    } else {
        return(pqte)
    }
}
