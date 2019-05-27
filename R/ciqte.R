utils::globalVariables(c("yname", "treat", "treated", "x", "xformla", "data", "wname", "probs", "method", "treated", "untreated", "eachIter"))
####Cross-sectional QTE method using Firpo (2007)########
#' @title compute.ci.qte
#'
#' @description \code{compute.ci.qte} computes the QTE under selection
#' on observables
#' 
#' @param qp QTEparams object containing the parameters passed to ciqte
#'
#' @import Hmisc
#'
#' @keywords internal
#' 
#' @return QTE object
compute.ci.qte <- function(qp) {

    setupData(qp)
    bootstrapiter <- qp$bootstrapiter

    ##don't use weights now
    ##no covariate att - will update if there are covariates
    #ate <- wtd.mean(treated[,yname], treated[,wname]) -
    #    wtd.mean(untreated[,yname], untreated[,wname])


    #treated.firpo.quantiles <- wtd.quantile(treated[,yname],
    #                                        treated[,wname],
    #                                        probs=probs,
    #                                        normwt=TRUE)
    #untreated.firpo.quantiles <- wtd.quantile(untreated[,yname],
    #                                          untreated[,wname],
    #                                          probs=probs,
    #                                          normwt=TRUE)

    ate <- getWeightedMean(treated[,yname], treated[,wname]) -
        getWeightedMean(untreated[,yname], untreated[,wname])

    treated.firpo.quantiles <- getWeightedQuantiles(probs, treated[,yname],
                                                   treated[,wname])

    untreated.firpo.quantiles <- getWeightedQuantiles(probs, untreated[,yname],
                                                     untreated[,wname])

    qte <- treated.firpo.quantiles - untreated.firpo.quantiles

    
    n = nrow(data)

    ##set these up to access later
    pscore.reg <- NULL
    if (!is.null(x)) {
        p <- rep(nrow(treated)/(nrow(treated) + nrow(untreated)), n)
        D <- data[,treat]
        y <- data[,yname]
        w <- data[,wname]
        ##estimate the propensity score
        this.formla <- y ~ x
        formula.tools::lhs(this.formla) <- as.name(treat)
        formula.tools::rhs(this.formla) <- formula.tools::rhs(xformla)
        pscore.reg <- glm(this.formla, data=data,
                          family=binomial(link=method))
        pscore <- fitted(pscore.reg)
        ##there are alternatives for how to compute the quantiles of 
        ##treated outcomes for the treated group:
        ##1) compute quantiles directly
        ##treated.quantiles = quantile(treated[,yname], probs=probs)
        ##2) use firpo method
        ##checkfun will be called by the various functions to be minimized
        ##in this routine

        treated.weights = w * D / pscore

        treated.firpo.quantiles <- getWeightedQuantiles(probs, y,
                                                        treated.weights, norm=TRUE)
        
        untreated.weights = w * (1-D) /  (1-pscore)

        untreated.firpo.quantiles <- getWeightedQuantiles(probs, y,
                                                          untreated.weights, norm=TRUE)
        
        qte <- treated.firpo.quantiles - untreated.firpo.quantiles

        ate <- getWeightedMean(y, treated.weights) -
            getWeightedMean(y, untreated.weights)
    }
    
    ##set up the distributions too if they are needed
    if (is.null(untreated.firpo.quantiles)) {
        F.treated.t.cf <- ecdf(untreated[,yname])
        F.treated.t <- ecdf(treated[,yname]) ##wrong name
    } else {
        F.treated.t.cf <- ecdf(untreated.firpo.quantiles)
        F.treated.t <- ecdf(treated.firpo.quantiles)  ##wrong name, but just need something in object
    }

    if (bootstrapiter) { ## do this to decrease size of each iteration
        out <- QTE(qte=qte, ate=ate, probs=probs)
    } else {
        out <- QTE(F.treated.t=F.treated.t,
                   F.treated.t.cf=F.treated.t.cf,
                   qte=qte, pscore.reg=pscore.reg,  ate=ate, probs=probs)
    }
    return(out)
}

#' @title ci.qte
#'
#' @description The \code{ci.qtet} method implements estimates the Quantile
#' Treatment Effect (QTE) under a Conditional Independence
#' Assumption (sometimes this is called Selection on Observables) developed
#' in Firpo (2007).  This method using propensity score re-weighting
#' and minimizes a check function to compute the QTET.  Standard errors
#' (if requested) are computed using the bootstrap.
#' 
#' @inheritParams panel.qtet
#' @param x Vector of covariates.  Default is no covariates
#' @param method Method to compute propensity score.  Default is logit; other
#'  option is probit.
#' @param w an additional vector of sampling weights
#' @param pl boolean for whether or not to compute bootstrap error in parallel.
#'  Note that computing standard errors in parallel is a new feature and may
#'  not work at all on Windows.
#' @param cores the number of cores to use if bootstrap standard errors are
#'  computed in parallel
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
#' q1 <- ci.qte(re78 ~ treat, x=NULL, data=lalonde.psid, se=FALSE,
#'  probs=seq(0.05,0.95,0.05))
#' summary(q1)
#' 
#' ##This estimation controls for all the available background characteristics.
#' q2 <- ci.qte(re78 ~ treat,
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  data=lalonde.psid, se=FALSE, probs=seq(0.05, 0.95, 0.05))
#' summary(q2)
#'
#' @return QTE object
#' @export
ci.qte <- function(formla, xformla=NULL, x=NULL, data, w=NULL,
                   probs=seq(0.05,0.95,0.05), se=TRUE,
                   iters=100, alp=0.05, method="logit",
                   retEachIter=FALSE, 
                   printIter=FALSE, pl=FALSE, cores=2) {

    qp <- QTEparams(formla, xformla, t=NULL, tmin1=NULL, tmin2=NULL, tname=NULL, data=data, w=w, idname=NULL, probs=probs, iters=iters, alp=alp, method=method, se=se, retEachIter=retEachIter, bootstrapiter=FALSE, pl=pl, cores=cores)
    ##setupData(qp) ##may be able to get rid of this too

    
    ##first calculate the actual estimate
    firpo.qte <- compute.ci.qte(qp)

    if (se) {

        qp$bootstrapiter <- TRUE

        ##bootstrap the standard errors
        SEobj <- bootstrap(qp, firpo.qte, compute.ci.qte)

        out <- QTE(qte=firpo.qte$qte, qte.upper=SEobj$qte.upper,
                   qte.lower=SEobj$qte.lower, ate=firpo.qte$ate,
                   ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                   qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                   c=SEobj$c,
                   pscore.reg=firpo.qte$pscore.reg,
                   F.treated.t=firpo.qte$F.treated.t,
                   F.treated.t.cf=firpo.qte$F.treated.t.cf,
                   eachIterList=eachIter,
                   probs=probs)
        return(out)
    } else {
        return(firpo.qte)
    }
}
