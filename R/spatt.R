utils::globalVariables(c("yname", "treat", "treated", "x", "wname", "probs", "method", "treated", "untreated", "eachIter"))


#####Semiparametric Difference in Differences#####
#' @title compute.spatt
#'
#' @description
#' \code{compute.spatt} implements the method of Abadie (2005); this is
#'  computed automatically in several other methods in the qte package
#'  but this function provides a standalone result when quantiles are not
#'  wanted/needed.
#' @description
#' \code{compute.ddid2} uses two periods of data (repeated cross sections
#' or panel) to estimate the Quantile Treatment Effect on the Treated (QTET)
#'
#' @param qp QTEparams object containing the parameters passed to ciqte
#'
#' @importFrom quantreg rq
#'
#' @return QTE object
#'
#' @keywords internal
compute.spatt <- function(qp) {

    setupData(qp)
    bootstrapiter <- qp$bootstrapiter

    ##calculate the att; this will be changed if there are covariates
    att = getWeightedMean(treated.t[,yname], treated.t[,wname]) -
        getWeightedMean(treated.tmin1[,yname], treated.tmin1[,wname]) -
                                             (getWeightedMean(untreated.t[,yname], untreated.t[,wname]) -
                                              getWeightedMean(untreated.tmin1[,yname], untreated.tmin1[,wname]))

    if(panel) {

        dta <- panel2cs(data, yname, qp$idname, tname)
        n <- nrow(dta)
        D <- dta[,treat]
        p <- sum(D)/n
        dy <- dta$dy
        ##estimate the propensity score
        this.formla <- y ~ x
        lhs(this.formla) <- as.name(treat)
        rhs(this.formla) <- rhs(xformla)
        pscore.reg <- glm(this.formla, data=dta,
                          family=binomial(link=method))
        pscore <- fitted(pscore.reg) ## TODO: does this make sense for repeated cross sections;
        ## above, I am just pooling both periods.

        waits <- (D-pscore)/(p*(1-pscore))


        ##TODO: notice that we are not accounting for sampling weight

        att <- getWeightedMean(dy, waits, norm=FALSE)

    } else {

   
    ##functionality with covariates is not yet implemented
        pscore.reg <- NULL #do this in case no covariates as we return this value
        if (!(is.null(x))) {
            ntt <- nrow(treated.t)
            nttmin1 <- nrow(treated.tmin1)
            nut <- nrow(untreated.t)
            nutmin1 <- nrow(untreated.tmin1)
            p <- (ntt+nttmin1)/(ntt+nttmin1+nut+nutmin1)
            D <- data[,treat]
            T <- 1*(data[,tname]==t)
            y <- data[,yname]
            w <- data[,wname]
            ##estimate the propensity score
            this.formla <- y ~ x
            lhs(this.formla) <- as.name(treat)
            rhs(this.formla) <- rhs(xformla)
            pscore.reg <- glm(this.formla, data=data,
                              family=binomial(link=method))
            pscore <- fitted(pscore.reg) ## TODO: does this make sense for repeated cross sections;
            ## above, I am just pooling both periods.
            
            lam <- (ntt+nut)/(ntt+nttmin1+nut+nutmin1) ## the fraction of observations in the last period
            
            waits1 <- (T-lam)/(lam*(1-lam))
            waits2 <- (D-pscore)/(p*(1-pscore))
            waits <- waits1 * waits2
            
            ##TODO: notice that we are not accounting for sampling weight
            
            att <- getWeightedMean(y, waits)
            
            ##att <- getWeightedMean(y, treated.weights) -
            ##    getWeightedMean(y, untreated.weights)
        }

    }

    out <- QTE(qte=NULL, pscore.reg=pscore.reg,  ate=att, probs=NULL)
    class(out) <- "QTE"
    return(out)
}


#' @title spatt
#'
#' @description \code{spatt} computes the Average Treatment Effect on the
#'  Treated (ATT) using the method of Abadie (2005)
#'
#' @inheritparams ciqte
#' @param formla The formula y ~ d where y is the outcome and d is the
#'  treatment indicator (d should be binary)
#' @param xformla A optional one sided formula for additional covariates that
#'  will be adjusted for.  E.g ~ age + education.  Additional covariates can
#'  also be passed by name using the x paramater.
#' @param t The 3rd time period in the sample (this is the name of the column)
#' @param tmin1 The 2nd time period in the sample (this is the name of the
#'  column)
#' @param tname The name of the column containing the time periods
#' @param x An optional vector of covariates (the name of the columns).
#'  Covariates can also be passed in formulat notation using the
#'  xformla paramter.
#' @param data The name of the data.frame that contains the data
#' @param panel Boolean indicating whether the data is panel or repeated cross
#'  sections
#' @param dropalwaystreated How to handle always treated observations
#'  in panel data case (not currently used)
#' @param idname The individual (cross-sectional unit) id name
#' @param iters The number of iterations to compute bootstrap standard errors.
#'  This is only used if se=TRUE
#' @param alp The significance level used for constructing bootstrap
#'  confidence intervals
#' @param method The method for estimating the propensity score when covariates
#'  are included
#' @param plot Boolean whether or not the estimated QTET should be plotted
#' @param se Boolean whether or not to compute standard errors
#' @param retEachIter Boolean whether or not to return list of results
#'  from each iteration of the bootstrap procedure
#' @param seedvec Optional value to set random seed; can possibly be used
#'  in conjunction with bootstrapping standard errors.
#' @param weights additional (usually sampling) weights to be applied (not yet
#'  implemented)
#'
#' @examples
#' ##load the data
#' data(lalonde)
#'
#' ## Run the panel.qtet method on the experimental data with no covariates
#' att1 <- spatt(re ~ treat, t=1978, tmin1=1975, tname="year",
#'  x=NULL, data=lalonde.psid.panel, idname="id", se=FALSE)
#' summary(att1)
#'
#' ## Run the panel.qtet method on the observational data with no covariates
#' 
#'
#' @references
#' Abadie (2005)
#'
#' @return \code{QTE} object
#' 
#' @export
spatt <- function(formla, xformla=NULL, t, tmin1,
                      tname, data, w=NULL, panel=FALSE,
                      idname=NULL, 
                      iters=100, alp=0.05, method="logit", plot=FALSE, se=TRUE,
                      retEachIter=FALSE, seedvec=NULL, pl=FALSE, cores=2) {
    
    qp <- QTEparams(formla=formla, xformla=xformla, t=t, tmin1=tmin1,
                    tname=tname, data=data, panel=panel,
                    idname=idname, probs=NULL,
                    iters=iters, alp=alp, method=method,
                    se=se, retEachIter=retEachIter, seedvec=seedvec,
                    pl=pl, cores=cores)
                    
    
    satt = compute.spatt(qp)

    if (se) {

        qp$bootstrapiter <- TRUE

        ##bootstrap the standard errors
        SEobj <- bootstrap(qp, satt, compute.spatt)

       
        ##could return each bootstrap iteration w/ eachIter
        ##but not currently doing that
        out <- QTE(qte=NULL, pscore.reg=satt$pscore.reg, ate=satt$ate,
                   ate.se=SEobj$ate.se, probs=NULL)
                   
        return(out)
    } else {
        return(satt)
    }

}
