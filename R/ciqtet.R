utils::globalVariables(c("yname", "treat", "treated", "x", "xformla", "data", "wname", "probs", "method", "treated", "untreated", "eachIter"))
####Cross-sectional QTET method using Firpo (2007)########
#' @title compute.ci.qtet
#'
#' @description \code{compute.ci.qtet} computes the QTT under
#' selection on observables
#' 
#' @inheritParams compute.ci.qte
#'
#' @import Hmisc
#' 
#' @keywords internal
#' 
#' @return QTE object
compute.ci.qtet = function(qp) {

    setupData(qp)
    bootstrapiter <- qp$bootstrapiter


    att <- wtd.mean(treated[,yname], treated[,wname]) -
        wtd.mean(untreated[,yname], untreated[,wname])

    treated.firpo.quantiles <- wtd.quantile(treated[,yname],
                                            treated[,wname],
                                            probs=probs, type="i/n")
    untreated.firpo.quantiles <- wtd.quantile(untreated[,yname],
                                              untreated[,wname],
                                              probs=probs, type="i/n")


    qte <- treated.firpo.quantiles - untreated.firpo.quantiles

    
    n <- nrow(data)

    ##set these up to access later
    pscore.reg <- NULL

    if (!is.null(x)) {

        ##estimate the propensity score
        D <- data[,treat]
        y <- data[,yname]
        w <- data[,wname]
        this.formla <- y ~ x
        formula.tools::lhs(this.formla) <- as.name(treat)
        formula.tools::rhs(this.formla) <- formula.tools::rhs(xformla)
        pscore.reg <- glm(this.formla, data=data,
                          family=binomial(link=method))
        pscore <- fitted(pscore.reg)
        pscore <- fitted(pscore.reg)
        p = rep(nrow(treated)/(nrow(treated) + nrow(untreated)), n)
        ##there are alternatives for how to compute the quantiles of 
        ##treated outcomes for the treated group:
        ##1) compute quantiles directly
        treated.quantiles <- treated.firpo.quantiles
        ##2) use firpo method
        ##checkfun will be called by the various functions to be minimized
        ##in this routine
        untreated.weights = w * (pscore/(1-pscore))*((1-D)/sum(D))

        untreated.firpo.quantiles = getWeightedQuantiles(probs, y,
                                                         untreated.weights)
        
        qte <- treated.firpo.quantiles - untreated.firpo.quantiles

        att <- wtd.mean(y, w) - wtd.mean(y, w*(1-D)*pscore/(p*(1-pscore)))
        

        ##Alternative method for calculating the distribution of each
        ##potential outcome using moment conditions.
        ##comment this out (unused)
        ##F.treated <- ecdf(treated[,yname])
        ##F.treatedcf.fun <- function(y) {
        ##    pterm <- pscore/((1-pscore)*p)
        ##    Dterm <- 1 - data[,treat]
        ##    yterm <- 1*(data[,yname] < y)
        ##    mean(pterm*Dterm*yterm)
        ##} #something appears to be off here for 0 wages, otherwise, everything good!
        ##y.seq <- seq(min(data[,yname]), max(data[,yname]), length.out=500)
        ##F.treatedcf = approxfun(y.seq,
        ##    vapply(y.seq, FUN=F.treatedcf.fun, FUN.VALUE=1)
        ##    , method="constant", yleft=0, yright=1, f=0, ties="ordered")
        ##class(F.treatedcf) = c("ecdf", "stepfun", class(F.treatedcf.fun))
        ##assign("nobs", nrow(treated), envir = environment(F.treatedcf))
    }

    F.treated.t <- ecdf(treated[,yname])
    F.treated.t.cf <- ecdf(untreated.firpo.quantiles)
    

    out <- QTE(F.treated.t=F.treated.t,
               F.treated.t.cf=F.treated.t.cf,
               qte=qte, pscore.reg=pscore.reg,  ate=att, probs=probs)
    return(out)
    
}

#' @title ci.qtet
#'
#' @description The \code{ci.qtet} method implements estimates the Quantile
#' Treatment Effect on the Treated (QTET) under a Conditional Independence
#' Assumption (sometimes this is called Selection on Observables) developed
#' in Firpo (2007).  This method using propensity score re-weighting
#' and minimizes a check function to compute the QTET.  Standard errors
#' (if requested) are computed using the bootstrap.
#' 
#' @inheritParams panel.qtet
#' @inheritParams ci.qte
#' @param method Method to compute propensity score.  Default is logit; other
#'  option is probit.
#' @param indsample Binary variable for whether to treat the samples as
#'  independent or dependent.  This affects bootstrap standard errors.  In
#'  the job training example, the samples are independent because they
#'  are two samples collected independently and then merged.  If the data is
#'  from the same source, usually should set this option to be FALSE.
#' @param printIter For debugging only; should leave at default FALSE unless
#'  you want to see a lot of output
#'
#' @references
#' Firpo, Sergio.   ``Efficient Semiparametric Estimation of Quantile Treatment
#'  Effects.'' Econometrica 75.1, pp. 259-276, 2015.
#' 
#' @examples
#' ## Load the data
#' data(lalonde)
#'
#' ##Estimate the QTET of participating in the job training program;
#' ##This is the no covariate case.  Note: Because individuals that participate
#' ## in the job training program are likely to be much different than
#' ## individuals that do not (e.g. less experience and less education), this
#' ## method is likely to perform poorly at estimating the true QTET
#' q1 <- ci.qtet(re78 ~ treat, x=NULL, data=lalonde.psid, se=FALSE,
#'  probs=seq(0.05,0.95,0.05))
#' summary(q1)
#' 
#' ##This estimation controls for all the available background characteristics.
#' q2 <- ci.qtet(re78 ~ treat, 
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  data=lalonde.psid, se=FALSE, probs=seq(0.05, 0.95, 0.05))
#' summary(q2)
#'
#' @return QTE object
#' @export
ci.qtet <- function(formla, xformla=NULL, w=NULL, data,
                    probs=seq(0.05,0.95,0.05), se=TRUE,
                 iters=100, alp=0.05, method="logit",
                 retEachIter=FALSE, indsample=TRUE,
                 printIter=FALSE, pl=FALSE, cores=2) {

    qp <- QTEparams(formla, xformla, t=NULL, tmin1=NULL, tmin2=NULL, tname=NULL, data=data, idname=NULL, probs=probs, iters=iters, alp=alp, method=method, se=se, retEachIter=retEachIter, bootstrapiter=FALSE, pl=pl, cores=cores)
    setupData(qp) ##may be able to get rid of this too

    
    ##first calculate the actual estimate
    firpo.qtet <- compute.ci.qtet(qp)

    if (se) {

        qp$bootstrapiter <- TRUE

        ##bootstrap the standard errors
        SEobj <- bootstrap(qp, firpo.qtet, compute.ci.qtet)

        out <- QTE(qte=firpo.qtet$qte, qte.upper=SEobj$qte.upper,
                   qte.lower=SEobj$qte.lower, ate=firpo.qtet$ate,
                   ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                   qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                   c=SEobj$c,
                   pscore.reg=firpo.qtet$pscore.reg,
                   F.treated.t=firpo.qtet$F.treated.t,
                   F.treated.t.cf=firpo.qtet$F.treated.t.cf,
                   eachIterList=eachIter,
                   probs=probs)
        return(out)
    } else {
        return(firpo.qtet)
    }

}
