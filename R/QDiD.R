utils::globalVariables(c("yname", "treat", "treated.t", "treated.tmin1", "untreated.t", "untreated.tmin1", "x", "wname", "probs", "method", "eachIter", "xformla"))


###Quantile Difference-in-Differences
##Note that you need to pass in data where treated status is noted in
##every period.  Data is form of (year-individual-outcome-x-evertreated)
#' @title Quantile Difference in Differences
#' @description
#' \code{compute.QDiD} computes the Quantile Difference in
#' Differences estimator
#' 
#' @param qp QTEparams object containing the parameters passed to QDiD
#' @return QTE object
#' @keywords internal
compute.QDiD <- function(qp) {

    setupData(qp)
    bootstrapiter <- qp$bootstrapiter

    ## this is the old way of adjusting for covariates, by just residualizing and then using unconditional method       ##adjust for covariates
    ##after adjustment then everything should proceed as before
    ## if (!(is.null(x))) {
    ##     cov.data <- data
    ##     cov.data$group <- as.factor(paste(cov.data[,treat],
    ##                                       cov.data[,tname],sep="-"))
    ##     group <- "group"
    ##     xmat = cov.data[,x]
    ##     first.stage <- lm(cov.data[,yname] ~ -1 + cov.data[,group] +
    ##                       as.matrix(xmat))
    ##     ##get residuals not including group dummies
    ##     bet <- coef(first.stage)[5:length(coef(first.stage))]
    ##     yfit <- cov.data[,yname] - as.matrix(xmat)%*%bet
    ##     data[,yname] <- yfit
    ## } 
    
    
    ##5) Compute Quantiles
    ##a) Quantiles of observed distribution
    q1 = stats::quantile(treated.t[,yname],probs=probs)
    q0 = stats::quantile(treated.tmin1[,yname] ,probs=probs) + stats::quantile(untreated.t[,yname] ,probs=probs) - stats::quantile(untreated.tmin1[,yname] ,probs=probs)
    
   
    ##7) Estimate ATT using A-I
     att = mean(treated.t[,yname]) - ( mean(treated.tmin1[,yname]) +
         mean(stats::quantile(untreated.t[,yname],
                  probs=stats::ecdf(treated.tmin1[,yname])(treated.tmin1[,yname]))) -
         mean(stats::quantile(untreated.tmin1[,yname],
                       probs=stats::ecdf(treated.tmin1[,yname])(treated.tmin1[,yname]))) )
   

    ##now with covariates
    pscore.reg <- NULL
    condQ.treated.t.qr <- NULL
    condQ.treated.t.cf.qr <- NULL
    if(!is.null(x)) {
                
        this.formla <- y ~ x ##just set up dummy formula first
        formula.tools::lhs(this.formla) <- as.name(yname)
        formula.tools::rhs(this.formla) <- formula.tools::rhs(xformla)

        ##condF.untreated.tmin1 <- dr(this.formla, untreated.tmin1,
        ##                            unique(untreated.tmin1[,yname]))

        tau <- probs

        taut <- seq(0,1, length.out=nrow(treated.t))
        taut <- seq(.01,.99,.01)

        uncondQTT <- TRUE

        ##need to run many more quantile regressions if you want unconditional
        ## quantiles
        if (uncondQTT==TRUE) {
            tau <- taut ##seq(0,1, length.out=nrow(untreated.tmin1))
        }
        condQ.untreated.tmin1.qr <- rq(this.formla, tau=tau,
                                       data=untreated.tmin1)

        condQ.untreated.t.qr <- rq(this.formla, tau=tau,
                                   data=untreated.t)

        condQ.treated.tmin1.qr <- rq(this.formla, tau=tau,
                                     data=treated.tmin1)

        condQ.treated.t.qr <- rq(this.formla, tau=tau,
                                 data=treated.t)


        
        condQ.treated.t.cf.qr <- list()
        condQ.treated.t.cf.qr$tau <- tau ##this should be the same as for treated.t and it is the last thing that is set, so ok
        condQ.treated.t.cf.qr$coefficients <- stats::coef(condQ.treated.tmin1.qr) + stats::coef(condQ.untreated.t.qr) -
            stats::coef(condQ.untreated.tmin1.qr)
        condQ.treated.t.cf.qr$terms <- condQ.treated.tmin1.qr$terms
        class(condQ.treated.t.cf.qr) <- "rqs"

        ## average the conditional distribution
        condF.dist <- stats::predict(condQ.treated.t.cf.qr, newdata=treated.t,
                              type="Fhat", stepfun=TRUE)

        yvals <- unique(treated.t[,yname])
        yvals <- yvals[order(yvals)]

        ##lf should contain 
        Ef <- function(y) {
            outvec <- c()
            for (i in 1:nrow(treated.t)) {
                outvec[i] <- condF.dist[[i]](y)
            }
            return(outvec)
        }

        uncFvals.all <- lapply(yvals, Ef)
        uncFvals <- unlist(lapply(uncFvals.all, mean))

        ##TODO: not sure if ok to reorder here or do it elsewhere...
        uncFvals <- uncFvals[order(uncFvals)]
        uncF <- makeDist(yvals, uncFvals, TRUE)

        q0 <- stats::quantile(uncF, probs, type=1)
                
        
        condQTT <- list()
        condQTT$tau <- probs
        condQTT$coefficients <- stats::coef(condQ.treated.t.qr) - stats::coef(condQ.treated.t.cf.qr)
        class(condQTT) <- "rqs"
    }

    if (bootstrapiter) {
        out <- QTE(ate=att, qte=q1-q0, probs=probs)
    } else{
        out <- QTE(ate=att, qte=(q1-q0), condQ.treated.t=condQ.treated.t.qr,
                   condQ.treated.t.cf=condQ.treated.t.cf.qr, probs=probs)
    }
    return(out)
}

##QDiD is a function that computes bootstrap
##standard errors for quantile treatment effects
#' @title Quantile Difference in Differences
#' @description \code{QDiD} is a Difference in Differences type method for
#' computing the QTET.
#'
#' The method can accommodate conditioning on covariates though it does so
#' in a restrictive way:  It specifies a linear model for outcomes conditional
#' on group-time dummies and covariates.  Then, after residualizing (see details
#' in Athey and Imbens (2006)), it computes the Change in Changes model
#' based on these quasi-residuals.
#' @inheritParams panel.qtet
#' @inheritParams CiC
#' @inheritParams ci.qte
#' 
#' @references
#' Athey, Susan and Guido Imbens.  ``Identification and Inference in Nonlinear
#'  Difference-in-Differences Models.'' Econometrica 74.2, pp. 431-497,
#'  2006.
#' 
#' @examples
#' ## load the data
#' data(lalonde)
#'
#' ## Run the Quantile Difference in Differences method conditioning on
#' ## age, education, black, hispanic, married, and nodegree
#' qd1 <- QDiD(re ~ treat, t=1978, tmin1=1975, tname="year",
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  data=lalonde.psid.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(qd1)
#'
#' @return QTE Object
#' 
#' @export
QDiD <- function(formla, xformla=NULL, t, tmin1, tname, x=NULL,data,
                 dropalwaystreated=TRUE, panel=FALSE, se=TRUE,
                 idname=NULL, method=NULL,
                 uniqueid=NULL, alp=0.05, probs=seq(0.05,0.95,0.05), iters=100,
                 retEachIter=FALSE, seedvec=NULL, printIter=F,
                 pl=FALSE, cores=NULL) {

    qp <- QTEparams(formla=formla, xformla=xformla, t=t, tmin1=tmin1,
                    tname=tname, data=data, panel=panel,
                    idname=idname, probs=probs,
                    iters=iters, bootstrapiter=FALSE, alp=alp, method=method,
                    se=se, retEachIter=retEachIter, seedvec=seedvec,
                    pl=pl, cores=cores)
    
    
    
    ## maybe, move this to setupData
    ##treated.t = data[data[,tname]==t & data[,treat]==1,]
    ##treated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 1, ]
    ##untreated.t = data[data[,tname]==t & data[,treat]==0,]
    ##untreated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 0, ]

    ##first calculate the actual estimate
    qdid = compute.QDiD(qp)


    if (se) {

        qp$bootstrapiter <- TRUE


        ##bootstrap the standard errors
        SEobj <- bootstrap(qp, qdid, compute.QDiD)

        ##could return each bootstrap iteration w/ eachIter
        ##but not currently doing that
        out <- QTE(qte=qdid$qte, qte.upper=SEobj$qte.upper,
                    qte.lower=SEobj$qte.lower, ate=qdid$ate,
                    ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                    qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                    F.treated.t=qdid$F.treated.t,
                    F.untreated.t=qdid$F.untreated.t,
                    F.treated.t.cf=qdid$F.treated.t.cf,
                    F.treated.tmin1=qdid$F.treated.tmin1,
                    F.treated.tmin2=qdid$F.treated.tmin2,
                    F.treated.change.tmin1=qdid$F.treated.change.tmin1,
                    F.untreated.change.t=qdid$F.untreated.change.t,
                    F.untreated.change.tmin1=qdid$F.untreated.change.tmin1,
                    F.untreated.tmin1=qdid$F.untreated.tmin1,
                   F.untreated.tmin2=qdid$F.untreated.tmin2,
                   condQ.treated.t=qdid$condQ.treated.t,
                   condQ.treated.t.cf=qdid$condQ.treated.t.cf,
                    pscore.reg=qdid$pscore.reg,
                    eachIterList=eachIter,
                    probs=probs)
        return(out)
    } else {
        return(qdid)
    }
}
