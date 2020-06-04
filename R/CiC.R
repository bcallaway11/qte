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
    ##now calculate the bootstrap confidence interval
    ## eachIter = list()
    ## ##Need to build dataset by sampling individuals, and then
    ## ##taking all of their time periods
    ## ##when it's a panel make draws by individual
    ## if (panel) {
    ##     ##all.ids = unique(data[,idname])
    ##     ##here we rely on having a balanced panel to get the right obs.
    ##     treated.t <- treated.t[order(treated.t[,idname]),]
    ##     treated.tmin1 <- treated.tmin1[order(treated.tmin1[,idname]),]
    ##     untreated.t <- untreated.t[order(untreated.t[,idname]),]
    ##     untreated.tmin1 <- untreated.tmin1[order(untreated.tmin1[,idname]),]
    ##     nt <- nrow(treated.t)
    ##     nu <- nrow(untreated.t)
    ##     ##out.bootdatalist <<- list()
    ##     for (i in 1:iters) {
    ##         ##reset boot.data
    ##         ##boot.data = data[0,]
    ##         if(!is.null(seedvec)) {
    ##             set.seed(seedvec[i])
    ##         }
    ##         randy.t = sample(1:nt, nt, replace=T)
    ##         randy.u <- sample(1:nu, nu, replace=T)
    ##         ##there has to be a way to do this faster, but go with the loop
    ##         ##for now
    ##         ##for (j in all.ids[randy]) {
    ##         ##    boot.data = rbind(boot.data, data[(data[,idname]==j),])
    ##         ##}
    ##         ##these.ids <- data[,idname][randy]
    ##         boot.data.treated.t <- treated.t[randy.t, ]
    ##         boot.data.treated.tmin1 <- treated.tmin1[randy.t, ]
    ##         boot.data.untreated.t <- untreated.t[randy.u, ]
    ##         boot.data.untreated.tmin1 <- untreated.tmin1[randy.u, ]
    ##         boot.data <- rbind(boot.data.treated.t, boot.data.untreated.t,
    ##                            boot.data.treated.tmin1,
    ##                            boot.data.untreated.tmin1)
    ##         ##boot.data = process.bootdata(boot.data, idname, uniqueid)
    ##         ##out.bootdatalist[[i]] <<- boot.data
    ##         thisIter = compute.CiC(qp)
    ##         ##already have a balanced panel so can increase speed by calling
    ##         ##with panel option set to F.
    ##         eachIter[[i]] = QTE(ate = thisIter$ate, qte=thisIter$qte,
    ##                     probs=probs)

    ##         if (printIter==T) {
    ##             print(i)
    ##         }
    ##     }
    ## } else { #make draws within each sample
    ##     treated.t = data[data[,tname]==t & data[,treat]==1,]
    ##     treated.tmin1 = data[data[,tname]==tmin1 & data[,treat]==1,]
    ##     untreated.t = data[data[,tname]==t & data[,treat]==0,]

    ##     untreated.tmin1 = data[data[,tname]==tmin1 & data[,treat]==0,]

    ##     for (i in 1:iters) {
    ##         if(!is.null(seedvec)) {
    ##             set.seed(seedvec[i])
    ##         }
    ##         n <- nrow(treated.t)
    ##         ran <- sample(1:n, n, replace=T)
    ##         boot.treated.t <- treated.t[ran,]

    ##         n <- nrow(treated.tmin1)
    ##         ran <- sample(1:n, n, replace=T)
    ##         boot.treated.tmin1 <- treated.tmin1[ran,]

    ##         n <- nrow(untreated.t)
    ##         ran <- sample(1:n, n, replace=T)
    ##         boot.untreated.t <- untreated.t[ran,]

    ##         n <- nrow(untreated.tmin1)
    ##         ran <- sample(1:n, n, replace=T)
    ##         boot.untreated.tmin1 <- untreated.tmin1[ran,]

    ##         boot.data <- rbind(boot.treated.t, boot.untreated.t,
    ##                            boot.treated.tmin1, boot.untreated.tmin1)
    ##         thisIter = compute.CiC(formla, xformla, t, tmin1, tname,
    ##             x, boot.data, 
    ##             dropalwaystreated, panel, idname, uniqueid, probs)
    ##         eachIter[[i]] = QTE(ate = thisIter$ate, qte=thisIter$qte,
    ##                     probs=probs)

    ##         if (printIter==T) {
    ##             print(i)
    ##         }
    ##     }
    
    ## }

    ## SEobj <- computeSE(eachIter, cic, alp=alp)

    ## if(!retEachIter) {
    ##     eachIter=NULL
    ## }

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
               c=SEobj$c,
               eachIterList=eachIter,
               probs=probs)
    return(out)
  } else {
    return(cic)
  }
}
