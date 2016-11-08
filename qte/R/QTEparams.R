#####SETUP CLASSES################
#' @title QTEparams
#'
#' @description \code{QTEparams} is an object that contains all the
#'  parameters passed to QTE methods
#'
#' @param formla Should be some y on treatment variable
#' @param xformla a formula for the other covariates such as ~ x1 + x2
#' @param t The last period (not always used)
#' @param tmin1 The last pre-treatment period (not always used)
#' @param tmin2 The 2nd to last pre-treatment period (not always used)
#' @param tname The name of the column containing time periods (not always used)
#' @param data The name of the data frame
#' @param idname The name of the id column used with panel data (not always used)
#' @param probs Which quantiles to produce quantile treatment effects for
#' @param iters The number of bootstrap iterations to use to compute standard errors
#' @param alp The significance level
#' @param method The method to compute the propensity score
#' @param plot boolean for whether or not to plot qtes
#' @param se boolean whether or not to compute standard errors
#' @param retEachIter boolean whether or not to return the full results from
#'  each bootstrap iteration
#' @param seedvec A seed to compute the same bootstrap standard errors each
#'  time the method is called (not always used)
#' 
#'
#' @export
QTEparams <- function(formla, xformla=NULL, t=NULL, tmin1=NULL, tmin2=NULL,
                      tname=NULL, data, 
                      idname=NULL, probs,
                      alp=NULL, method=NULL, plot=NULL, se=NULL, iters=NULL,
                      retEachIter=NULL, bootstrapiter=NULL, seedvec=NULL,
                      pl=NULL, cores=NULL) {
    out <- list(formla=formla, xformla=xformla, t=t, tmin1=tmin1,
                tmin2=tmin2, tname=tname, data=data, weights=weights,
                idname=idname, probs=probs, method=method,
                plot=plot, se=se, iters=iters, retEachIter=retEachIter,
                bootstrapiter=bootstrapiter, seedvec=seedvec, pl=pl, cores=cores)
    class(out) <- "QTEparams"
    return(out)
}
