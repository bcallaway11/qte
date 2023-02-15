###Change in Changes (Athey-Imbens-2006)
##Note that you need to pass in data where treated status is noted in
##every period.  Data is form of (year-individual-outcome-x-evertreated)
#'@title athey.imbens
#'
#' @description
#' \code{compute.CiC} does the computational
#' work for the Change in Changes model
#' of Athey and Imbens, 2006.
#' 
#' @inheritParams panel.qtet
#'
#' @keywords internal
#' @export
compute.CiC <- function(qp) {

  setupData(qp)
  
  ##just to make sure the factors are working ok
  data = droplevels(data)

  ## will update this if there are covariates...
  F.treatedcf.t <- ecdf(quantile(untreated.t[,yname],
                                 probs=F.untreated.tmin1(treated.tmin1[,yname]), type=1))
  att = mean(treated.t[,yname]) -
    mean(quantile(untreated.t[,yname],
                  probs=F.untreated.tmin1(treated.tmin1[,yname]),type=1)) #See A-I p.441 Eq. 16

  QR0tQ <- NULL
  QR1t <- NULL
  ##adjust for covariates
  ##after adjustment then everything should proceed as before
  if (!(is.null(xformla))) {

    u <- seq(.01,.99,.01)
    n1t <- nrow(treated.t)
    n1tmin1 <- nrow(treated.tmin1)
    n0t <- nrow(untreated.t)
    n0tmin1 <- nrow(untreated.tmin1)

    yformla <- toformula("y", rhs.vars(xformla))
    QR0t <- rq(yformla, data=untreated.t, tau=u)
    QR0tmin1 <- rq(yformla, data=untreated.tmin1, tau=u)
    QR1t <- rq(yformla, data=treated.t,tau=u)

    ## in athey and imbens, think about k(cic) transformation
    QR0tmin1F <- predict(QR0tmin1, newdata=treated.tmin1, type="Fhat", stepfun=TRUE)
    F0tmin1 <- sapply(1:n1tmin1, function(i) QR0tmin1F[[i]](treated.tmin1$y[i]))
    QR0tQ <- predict(QR0t, newdata=treated.tmin1, type="Qhat", stepfun=TRUE)
    y0t <- sapply(1:n1tmin1, function(i) QR0tQ[[i]](F0tmin1[i]))## these are pseudo counterfactual outcomes (in the sense that they share the same distribution as Y_t(0) but are not necessarily equal)

    F.treatedcf.t <- ecdf(y0t)
    
    att <- mean(treated.t[,yname]) - mean(y0t)
    ## old regression-based approach
    ## cov.data <- data
    ## cov.data$group <- as.factor(paste(cov.data[,treat],
    ##                                   cov.data[,tname],sep="-"))
    ## group <- "group"
    ## xmat = cov.data[,x]
    ## first.stage <- lm(cov.data[,yname] ~ -1 + cov.data[,group] +
    ##                   as.matrix(xmat))
    ## ##get residuals not including group dummies
    ## bet <- coef(first.stage)[5:length(coef(first.stage))]
    ## yfit <- cov.data[,yname] - as.matrix(xmat)%*%bet
    ## data[,yname] <- yfit
  }    

  
  ##5) Compute Quantiles
  ##a) Quantiles of observed distribution
  q1 = quantile(treated.t[,yname],probs=probs,type=1)
  q0 = quantile(F.treatedcf.t,probs=probs,type=1)
  
  
  
  out <- QTE(F.treated.t = F.treated.t, F.treated.t.cf = F.treatedcf.t,
             F.treated.tmin1=F.treated.tmin1,
             F.untreated.t=F.untreated.t,
             F.untreated.tmin1=F.untreated.tmin1,
             condQ.treated.t.cf=QR0tQ,
             condQ.treated.t=QR1t,
             ate=att, qte=(q1-q0), probs=probs)
  class(out) <- "QTE"
  
  return(out)
}


##CiC is a function that computes bootstrap
##standard errors for quantile treatment effects
#' @title Change in Changes
#'
#' @description \code{CiC} computes the Quantile Treatment Effect on the
#'  Treated (QTET) using the method of Athey and Imbens (2006).  \code{CiC}
#'  is a Difference in Differences type method.  It requires
#'  having two periods of data that can be either  repeated cross sections
#'  or panel data.
#'
#' The method can accommodate conditioning on covariates though it does so
#' in a restrictive way:  It specifies a linear model for outcomes conditional
#' on group-time dummies and covariates.  Then, after residualizing (see details
#' in Athey and Imbens (2006)), it computes the Change in Changes model
#' based on these quasi-residuals.
#'
#' @inheritParams panel.qtet
#' @param panel Binary variable indicating whether or not the dataset is
#'  panel.  This is used for computing bootstrap standard errors correctly.
#'
#' @examples
#' ## load the data
#' data(lalonde)
#' ## Run the Change in Changes model conditioning on age, education,
#' ## black, hispanic, married, and nodegree
#' c1 <- CiC(re ~ treat, t=1978, tmin1=1975, tname="year",
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  data=lalonde.psid.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(c1)
#' 
#'
#' @return QTE Object
#'
#' @references
#' Athey, Susan and Guido Imbens.  ``Identification and Inference in Nonlinear
#'  Difference-in-Differences Models.'' Econometrica 74.2, pp. 431-497,
#'  2006.
#' 
#' @export
CiC <- function(formla, xformla=NULL, t, tmin1, tname, data,
                panel=FALSE,
                se=TRUE, idname=NULL, 
                alp=0.05, probs=seq(0.05,0.95,0.05), iters=100,
                pl=FALSE, cores=2,
                retEachIter=FALSE) {
  
  if (panel) {
    data <- panelize.data(data, idname, tname, t, tmin1)
  } else { ## repeated cross sections case
    data <- subset(data, (data[,tname]==tmin1 | data[,tname]==t))
  }
  
  qp <- QTEparams(formla=formla, xformla=xformla, t=t, tmin1=tmin1,
                  tname=tname, data=data, panel=panel,
                  idname=idname, probs=probs,
                  iters=iters, bootstrapiter=FALSE, alp=alp, 
                  se=se, retEachIter=retEachIter, 
                  pl=pl, cores=cores)

  if (panel) {
    panel.checks(qp)
  }
  
  ##first calculate the actual estimate
  cic = compute.CiC(qp)

  if (se) {

    qp$bootstrapiter <- TRUE


    ##bootstrap the standard errors
    SEobj <- bootstrap(qp, cic, compute.CiC)


    out <- QTE(qte=cic$qte, qte.upper=SEobj$qte.upper,
               F.treated.t=cic$F.treated.t,
               F.untreated.t=cic$F.untreated.t,
               F.treated.t.cf=cic$F.treated.t.cf,
               F.treated.tmin1=cic$F.treated.tmin1,
               F.untreated.tmin1=cic$F.untreated.tmin1,
               condQ.treated.t=cic$condQ.treated.t,
               condQ.treated.t.cf=cic$condQ.treated.t.cf,
               qte.lower=SEobj$qte.lower, ate=cic$ate,
               ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
               qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
               c=SEobj$c, alp=alp, 
               eachIterList=eachIter,
               probs=probs)
    return(out)
  } else {
    return(cic)
  }
}

#' @title cic_attgt
#'
#' @inheritParams pte::did_attgt
#'
#' @return pte::attgt_noif object.  `cic_attgt` computes attgt using
#'  the CIC approach.  It also returns distributions of observed outcomes
#'  for the treated group (F1), the counterfactual distribution
#'  of untreated potential potential outcomes for the treated group (F0),
#'  and the distribution of the treatment effect under the assumption
#'  of rank invariance over time (Fte) all through the extra_gt_returns
#'  argument to `pte::attgt_noif` object.
#'
#' @export
cic_attgt <- function(gt_data, xformla=~1, ...) {

  #-----------------------------------------------------------------------------
  # handle covariates
  #-----------------------------------------------------------------------------
  # for outcome regression, get pre-treatment values
  Xpre <- model.frame(xformla, data=subset(gt_data,name=="pre"))

  # convert two period panel into one period
  gt_data_outcomes <- tidyr::pivot_wider(gt_data[,c("D","id","period","name","Y")], id_cols=c(id, D),
                                           names_from=c(name),
                                           values_from=c(Y))

  # merge outcome and covariate data
  gt_dataX <- cbind.data.frame(gt_data_outcomes, Xpre)

  # treatment dummy variable
  D <- gt_dataX$D

  # pre- and post-treatment outcomes
  Y_post <- gt_dataX$post
  Y_pre <- gt_dataX$pre

  # drop missing levels to be safe
  gt_dataX <- droplevels(gt_dataX)

  #-----------------------------------------------------------------------------
  # make computations, this code is copied / slightly modified
  # from compute.CiC function
  #-----------------------------------------------------------------------------

  
  # will update this if there are covariates...
  kcic <- quantile(Y_post[D==0], probs=ecdf(Y_pre[D==0])(Y_pre[D==1]), type=1)
  att <- mean(Y_post[D==1]) - mean(kcic)

  F0 <- ecdf(kcic)
  F1 <- ecdf(Y_post[D==1])

  # distribution of the treatment effect under rank invariance
  Fte <- ecdf(Y_post[D==1] - kcic)

  # adjust for covariates
  if (length(rhs.vars(xformla)) > 0) {

    u <- seq(.01,.99,.01)
    n1 <- sum(D)
    n0 <- sum(1-D)

    post_formula <- BMisc::toformula("post", rhs.vars(xformla))
    pre_formula <- BMisc::toformula("pre", rhs.vars(xformla))
    QR0t <- rq(post_formula, data=gt_dataX[D==0,], tau=u)
    QR0tmin1 <- rq(pre_formula, data=gt_dataX[D==0,], tau=u)
    QR1t <- rq(post_formula, data=gt_dataX[D==1,],tau=u)

    # compute counterfactual outcomes
    QR0tmin1F <- predict(QR0tmin1, newdata=gt_dataX[D==1,], type="Fhat", stepfun=TRUE)
    F0tmin1 <- sapply(1:n1, function(i) QR0tmin1F[[i]](gt_dataX[D==1,]$pre[i]))
    
    # check for violations of support conditions
    if ( mean( F0tmin1 >=.99 ) + mean(F0tmin1 <= .01) > 0.1 ) warning("lots of very high/low \"ranks\" for treated units => CIC support conditions are likely violated...")

    QR0tQ <- predict(QR0t, newdata=gt_dataX[D==1,], type="Qhat", stepfun=TRUE)
    y0t <- sapply(1:n1, function(i) QR0tQ[[i]](F0tmin1[i]))## these are pseudo counterfactual outcomes (in the sense that they share the same distribution as Y_t(0) but are not necessarily equal)

    F0 <- ecdf(y0t)
    
    att <- mean(Y_post[D==1]) - mean(y0t)

    Fte <- ecdf(Y_post[D==1] - y0t)

  }    

  # return attgt
  attgt_noif(attgt=att, extra_gt_returns=list(F0=F0, F1=F1, Fte=Fte))
}


#' @title cic2
#'
#' @description This is a multi-period implementation of the change-in-changes
#'  approach from Athey and Imbens (2006, Econometrica).  This function
#'  is in a beta release and users should use caution when using this function
#'  in emprical work.
#'
#'  The function builds on the `pte` package and will return an overall
#'  treatment effect parameter as well as an event study.  See, in particular,
#'  the argument `ret_quantile` below.
#'
#' @inheritParams pte::pte
#' @param ret_quantile This parameter determines which quantile will be reported
#'  by the cic2 function.  By default `ret_quantile=NULL`; in this case, the
#'  function will return an estimate of the overall ATT and an event study for
#'  the ATT.  Other choices should be between 0 and 1.  For example, if the
#'  user specifies `ret_quantile=0.9`, then the function will return overall
#'  and event study parameters for the QTT(0.9).  These ...would be better to return the overall distribution and then to average and invert in later steps...
#' @param ret_dist If set to be true, the function returns the observed
#'  distribution of outcomes and counterfactual distribution of outcomes
#'  for each (g,t) through the `extra_gt_returns` element of `group_time_att`
#'  object.
#'
#' @export
cic2 <- function(yname,
                 gname,
                 tname,
                 idname,
                 data,
                 xformla=~1,
                 ret_quantile=NULL,
                 gt_type="att",
                 anticipation=0,
                 cband=TRUE,
                 alp=0.05,
                 boot_type="empirical",
                 biters=100,
                 cl=1) {

  if (boot_type != "empirical") {
    stop("only empirical bootstrap currently implemented")
  }
  
  res <- pte(yname=yname,
             gname=gname,
             tname=tname,
             idname=idname,
             data=data,
             setup_pte_fun=setup_pte,
             subset_fun=two_by_two_subset,
             attgt_fun=cic_attgt,
             xformla=xformla,
             anticipation=anticipation,
             cband=cband,
             alp=alp,
             boot_type=boot_type,
             biters=biters,
             cl=cl,
             ret_quantile=ret_quantile,
             gt_type=gt_type)

  res
}

