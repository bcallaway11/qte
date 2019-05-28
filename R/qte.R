#####MAIN FUNCTIONS#####

#####Panel QTET#####
##Idea here is that we can use information from a third period
##to point identify counterfactual distribution of outcomes
##for the treated group
##call plot function, summary function, formula function, etc. later
##add functionality to pass in pscore
#' @title compute.panel.qtet
#'
#' @description
#' \code{compute.panel.qtet} uses third period of data,
#' combined with Distributional
#' Difference in Differences assumption (Fan and Yu, 2012)
#' to point identify QTET.
#' 
#' @inheritParams compute.ci.qte
#' 
#' @return QTE object
#'
#' @export 
compute.panel.qtet <- function(qp) {
    
    setupData(qp)


    ##calculate the distribution of the change for the treated group;
    ## this will be changed if there are covariates
    F.treated.change.t <- F.untreated.change.t

    
    ##calculate the att; this will be changed if there are covariates
    att = mean(treated.t[,yname]) - mean(treated.tmin1[,yname]) -
        (mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname]))


    ## build counterfactual distribution depending on which case we are in
    
    ##a.1) If there are covariates need to satisfy the Distributional D-i-D
    ##then we will need to modify the distribution of the changes in outcomes
    ##using the method presented in the paper.
    
    pscore.reg <- NULL #do this in case no covariates as we return this value
    qr.reg <- NULL

    ## setup the data some
    if (!(is.null(x))) {    
        treated.t$dy = treated.t[,yname] - treated.tmin1[,yname]
        treated.tmin1$dy <- treated.tmin1[,yname] - treated.tmin2[,yname]
        untreated.t$dy = untreated.t[,yname] - untreated.tmin1[,yname]
        untreated.tmin1$dy <- untreated.tmin1[,yname] -
            untreated.tmin2[,yname]
    }

    if (is.null(x)) method <- "pscore" ## just trick to reuse some code in pscore and no covariates case

    if (method == "pscore") {
    
        if (!(is.null(x))) {

            ## estimate pscore
            this.formla <- y ~ x
            formula.tools::rhs(this.formla) <- formula.tools::rhs(xformla)
            formula.tools::lhs(this.formla) <- as.name(treat)
            pscore.data  <-  rbind(treated.tmin2, untreated.tmin2)
            pscore.reg <- glm(this.formla, data=pscore.data,
                              family=binomial(link="logit"))
            pscore <- fitted(pscore.reg)

            

            ## calculate ATT and QTT
            dtat <- rbind.data.frame(treated.t, untreated.t)
            dtat$pscore <- pscore
            pD1 = nrow(treated.t)/nrow(untreated.t)

            ##this contains the support of the change in y
            p.dy.seq <- unique(dtat$dy)
            p.dy.seq <- p.dy.seq[order(p.dy.seq)]
            distvals <- rep(0, length(p.dy.seq))
            for (i in 1:length(p.dy.seq)) {
                distvals[i] <- getWeightedMean( 1*(dtat$dy <= p.dy.seq[i]), (1-dtat[,treat])*pscore/((1-pscore)*pD1), norm=TRUE)
                ## mean(1*(pscore.data$changey <= p.dy.seq[i])*
                ##             (1-pscore.data[,treat])*pscore/((1-pscore)*pD1)) /
                ## mean( (1-pscore.data[,treat])*pscore/((1-pscore)*pD1) )
            }
            
            F.untreated.change.t  <- makeDist(p.dy.seq, distvals)        

            ##after we have the propensity score 
            ##use it to estimate the att using abadie's method.
            att <- getWeightedMean(y=dtat$dy,
                                   weights=(dtat[,treat]-pscore)/((1-pscore)*pD1), norm=TRUE)## mean(((pscore.data$changey)/pval)*(pscore.data[,treat] - pscore) /
            ##      (1-pscore))

            ## update the lag of the untreated change so that we can
            ## do pre-testing if desired
            pscore.data.tmin1 <- rbind(treated.tmin1, untreated.tmin1)
            posvals.seq <- pscore.data.tmin1$dy
            distvals.tmin1 <- rep(0, length(posvals.seq))
            for (dy in posvals.seq) {
                distvals.tmin1[which(dy==posvals.seq)] =
                    mean(1*(pscore.data.tmin1$dy<=dy)*
                         (1-pscore.data.tmin1[,treat])*pscore/((1-pscore)*pD1))
            }
            pscore.data.tmin1$distvals <- distvals.tmin1
            pscore.data1.tmin1 <- pscore.data.tmin1[order(pscore.data.tmin1$dy),]
            F.untreated.change.tmin1 <- makeDist(pscore.data1.tmin1$dy,
                                                 pscore.data1.tmin1$distvals)

            
            
        }
        
        
        ##compute counterfactual distribution
        quantys1 <- quantile(F.treated.tmin1,
                             probs=F.treated.tmin2(treated.tmin2[,yname]))

        quantys2 <- quantile(F.untreated.change.t,
                             probs=F.treated.change.tmin1(treated.tmin1[,yname] -
                                                          treated.tmin2[,yname]))

        y.seq <- (quantys1+quantys2)[order(quantys1 + quantys2)]

        F.treated.t.cf.val <- vapply(y.seq,
                                     FUN=function(y) { mean(1*(quantys2 <=
                                                               y - quantys1)) }, FUN.VALUE=1)

        F.treated.t.cf <- makeDist(y.seq, F.treated.t.cf.val)

    } else if (method == "qr") {
        u <- seq(.01,.99,.01) ## hard-coded for now
        yformla <- BMisc::toformula("y", BMisc::rhs.vars(xformla))
        dyformla <- BMisc::toformula("dy", BMisc::rhs.vars(xformla))

        dQRt <- quantreg::rq(dyformla, data=untreated.t, tau=u) ## holds by conditional did assumption
        dQRtmin1 <- quantreg::rq(dyformla, data=treated.tmin1, tau=u)
        QRtmin1 <- quantreg::rq(yformla, data=treated.tmin1, tau=u)
        QRtmin2 <- quantreg::rq(yformla, data=treated.tmin2, tau=u)

        ## use csa-type result; exploit that we average over X_i
        n1 <- nrow(treated.t)
        n0 <- nrow(untreated.t)
        QRtmin2F <- predict(QRtmin2, type="Fhat", stepfun=TRUE)
        Ftmin2 <- sapply(1:n1, function(i) QRtmin2F[[i]](treated.tmin2$y[i]))
        QRtmin1Q <- predict(QRtmin1, type="Qhat", stepfun=TRUE)
        Qtmin1 <- sapply(1:n1, function(i) QRtmin1Q[[i]](Ftmin2[i]))

        dQRtmin1F <- predict(dQRtmin1, type="Fhat", stepfun=TRUE)
        dFtmin1 <- sapply(1:n1, function(i) dQRtmin1F[[i]](treated.tmin1$dy[i]))
        dQRtQ <- predict(dQRt, newdata=treated.t, type="Qhat", stepfun=TRUE) ## predict for the treated guys even though estimate with the untreated groups

        dQt <- sapply(1:n1, function(i) dQRtQ[[i]](dFtmin1[i]))

        yvals <- unique( rbind.data.frame(treated.t, untreated.t)$y )
        yvals <- sort(yvals)
        F.treated.t.cf.val <- sapply(yvals, function(yy) mean(1*(dQt + Qtmin1 <= yy)))

        F.treated.t.cf <- makeDist(yvals, F.treated.t.cf.val)

        att <- mean(treated.t$y) - sum(quantile( F.treated.t.cf, probs=u, type=1 ))/length(u)
    } else {
        stop("invalid method supplied")
    }
    

    ##QTE
    F.treated.t <- ecdf(treated.t[,yname])

    qte <- quantile(F.treated.t, probs=probs) -
        quantile(F.treated.t.cf, probs=probs)


    out <- QTE(F.treated.t=F.treated.t,
               F.treated.tmin1=F.treated.tmin1,
               F.treated.tmin2=F.treated.tmin2,
               F.treated.change.tmin1=F.treated.change.tmin1,
               F.untreated.t=F.untreated.t,
               F.untreated.tmin1=F.untreated.tmin1,
               F.untreated.tmin2=F.untreated.tmin2,
               F.untreated.change.t=F.untreated.change.t,
               F.untreated.change.tmin1=F.untreated.change.tmin1,
               F.treated.t.cf=F.treated.t.cf,
               qte=qte, pscore.reg=pscore.reg,  ate=att, probs=probs)
    class(out) <- "QTE"
    return(out)
}


#' @title panel.qtet
#'
#' @description \code{panel.qtet} computes the Quantile Treatment Effect
#' on the Treated (QTET) using the method of Callaway and Li (2015).  This
#' method should be used when the researcher wants to invoke a Difference
#' in Differences assumption to identify the QTET.  Relative to the other
#' Difference in Differences methods available in the \code{qte} package,
#' this method's assumptions are more intuitively similar to the identifying
#' assumptions used in identifying the Average Treatment Effect on the Treated
#' (ATT).
#'
#' Additionally, this method can accommodate covariates in a more
#' flexible way than the other Difference in Differences methods available.
#' In order to accommodate covariates, the user should specify a vector \code{x}
#' of covariate names.  The user also may specify a method for estimating
#' the propensity score.  The default is logit.
#'
#' \code{panel.qtet} can only be used in some situations, however.  The
#' method requires three periods of panel data where individuals
#' are not treated until the last period.  The data should be formatted
#' as a panel; the names of columns containing time periods and ids
#' for each cross sectional unit need to be passed to the method.
#'
#' @param formla The formula y ~ d where y is the outcome and d is the
#'  treatment indicator (d should be binary), d should be equal to one
#'  in all time periods for individuals that are eventually treated
#' @param xformla A optional one sided formula for additional covariates that
#'  will be adjusted for.  E.g ~ age + education.  Additional covariates can
#'  also be passed by name using the x paramater.
#' @param t The 3rd time period in the sample.  Treated individuals should
#'  be treated in this time period and untreated individuals should not be
#'  treated.  The code attempts to enforce this condition, but it is good
#'  try to handle this outside the panel.qtet method.
#' @param tmin1 The 2nd time period in the sample.  This should be a
#'  pre-treatment period for all individuals in the sample.
#' @param tmin2 The 1st time period in the sample.  This should be a
#'  pre-treatment period for all individuals in the sample.
#' @param tname The name of the column containing the time periods
#' @param data A data.frame containing all the variables used
#' @param idname The individual (cross-sectional unit) id name
#' @param probs A vector of values between 0 and 1 to compute the QTET at
#' @param iters The number of iterations to compute bootstrap standard errors.
#'  This is only used if se=TRUE
#' @param alp The significance level used for constructing bootstrap
#'  confidence intervals
#' @param method The method for including covariates, should either be "QR"
#'  for quantile regression or "pscore" for propensity score
#' @param se Boolean whether or not to compute standard errors
#' @param retEachIter Boolean whether or not to return list of results
#'  from each iteration of the bootstrap procedure (default is FALSE).
#'  This is potentially useful for debugging but can cause errors due
#'  to running out of memory.
#' @param pl Whether or not to compute standard errors in parallel
#' @param cores Number of cores to use if computing in parallel
#'
#' @examples
#' ##load the data
#' data(lalonde)
#'
#' ## Run the panel.qtet method on the experimental data with no covariates
#' pq1 <- panel.qtet(re ~ treat, t=1978, tmin1=1975, tmin2=1974, tname="year",
#'  data=lalonde.exp.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(pq1)
#'
#' ## Run the panel.qtet method on the observational data with no covariates
#' pq2 <- panel.qtet(re ~ treat, t=1978, tmin1=1975, tmin2=1974, tname="year",
#'  data=lalonde.psid.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(pq2)
#'
#' ## Run the panel.qtet method on the observational data conditioning on
#' ## age, education, black, hispanic, married, and nodegree.
#' ## The propensity score will be estimated using the default logit method.
#' pq3 <- panel.qtet(re ~ treat, t=1978, tmin1=1975, tmin2=1974, tname="year",
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  data=lalonde.psid.panel, idname="id", se=FALSE, method="pscore",
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(pq3)
#'
#' pq4 <- panel.qtet(re ~ treat, t=1978, tmin1=1975, tmin2=1974, tname="year",
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  data=lalonde.psid.panel, idname="id", se=FALSE, method="qr",
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(pq4)
#' 
#' @references
#' Callaway, Brantly and Tong Li.  ``Quantile Treatment Effects in Difference
#'  in Differences Models with Panel Data.'' Working Paper, 2019.
#'
#' @return \code{QTE} object
#' 
#' @export
panel.qtet <- function(formla, xformla=NULL, t, tmin1, tmin2,
                      tname, data, 
                      idname, probs=seq(0.05,0.95,0.05),
                      iters=100, alp=0.05, method=c("qr","pscore"), se=TRUE,
                      retEachIter=FALSE, pl=FALSE, cores=NULL) {

    method <- method[1]

    data <- panelize.data(data, idname, tname, t, tmin1, tmin2)   

    qp <- QTEparams(formla=formla, xformla=xformla,
                    t=t, tmin1=tmin1, tmin2=tmin2,
                    tname=tname, data=data,
                    idname=idname, probs=probs,
                    iters=iters, alp=alp, method=method,
                    se=se, retEachIter=retEachIter, 
                    pl=pl, cores=cores, panel=TRUE, bootstrapiter=FALSE)
    

    
    ## Do some error handling for common cases...

    ## setup the data as build on this a bit...
    setupData(qp)

    ## do some checking that format of data ok
    panel.checks(qp)

    ##first calculate the actual estimate
    pqte = compute.panel.qtet(qp)


    ## compute standard errors if so desired
    if (se) {

        qp$bootstrapiter <- TRUE

        ##bootstrap the standard errors
        ## the bootstrap method does this generically
        SEobj <- bootstrap(qp, pqte, compute.panel.qtet)

        ## set the results
        out <- QTE(qte=pqte$qte, qte.upper=SEobj$qte.upper,
                   qte.lower=SEobj$qte.lower, ate=pqte$ate,
                   ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                   qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                   c=SEobj$c,
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


######GENERAL HELPER FUNCTIONS#######

##return an SE object
##bootIters should contain ATT as first object in list
#' @title computeDiffSE
#'
#' @description Takes two sets of initial estimates and bootstrap
#'  estimations
#'  (they need to have the same number of iterations) and determines
#'  whether or not the estimates are statistically different from each
#'  other.  It can be used to compare any sets of estimates, but it is
#'  particularly used here to compare estimates from observational methods
#'  with observations from the experimental data (which also have standard
#'  errors because, even though the estimates are cleanly identified, they
#'  are still estimated).
#'
#' @param est1 A QTE object containing the first set of estimates
#' @param bootIters1 A List of QTE objects that have been bootstrapped
#' @param est2 A QTE object containing a second set of estimates
#' @param bootIters2 A List of QTE objects that have been bootstrapped
#'  using the second method
#' @inheritParams panel.qtet
#'
#' @export
computeDiffSE <- function(est1, bootIters1, est2, bootIters2, alp=0.05) {
    iters <- length(bootIters1)
    ate.diff <- est1$ate - est2$ate
    qte.diff <- est1$qte - est2$qte
    ##For now, just plot the qte and att with standard errors
    ##helper function to get the first element out of a list
    getElement <- function(Lst, elemNum) {
        return(as.numeric(unlist((Lst[elemNum])))) #as.numeric is a trick to 
        ##get numerical value of qte
    }
    all.ate1 = unlist(sapply(bootIters1, FUN=getElement,elemNum=2))
    all.ate2 = unlist(sapply(bootIters2, FUN=getElement,elemNum=2))
    all.ate.diff <- all.ate1 - all.ate2
    ##get se
    ate.diff.se <- sd(all.ate.diff)
    ##reorder asc
    all.ate.diff = all.ate.diff[order(all.ate.diff)]
    ate.diff.upper = all.ate.diff[min(iters,round((1-alp/2)*iters))]
    ate.diff.lower = all.ate.diff[max(1,round((alp/2)*iters))]
    
    ##now get CI for qte:
    all.qte1 = lapply(bootIters1, FUN=getElement, elemNum=1)
    all.qte2 = lapply(bootIters2, FUN=getElement, elemNum=1)
    ##all.qte.diff <- all.qte1 - all.qte2
    qte1.mat = do.call(rbind,lapply(all.qte1, FUN=as.numeric, ncol=length(all.qte1[[1]]), byrow=TRUE))
    qte2.mat = do.call(rbind,lapply(all.qte2, FUN=as.numeric, ncol=length(all.qte2[[1]]), byrow=TRUE))
    qte.diff.mat = qte1.mat - qte2.mat
    ##standard error
    qte.diff.se <- apply(qte.diff.mat, FUN=sd, MARGIN=2)
    ##order each column
    sorted.qte.diff.mat = apply(qte.diff.mat, 2, sort)
    qte.diff.upper = sorted.qte.diff.mat[round((1-alp/2)*iters),]
    qte.diff.lower = sorted.qte.diff.mat[max(1,round((alp/2)*iters)),]

    out <- list(ate.diff=ate.diff, qte.diff=qte.diff,
                ate.diff.se=ate.diff.se,
                ate.diff.upper=ate.diff.upper, ate.diff.lower=ate.diff.lower,
                qte.diff.se=qte.diff.se,
                qte.diff.upper=qte.diff.upper, qte.diff.lower=qte.diff.lower)
    class(out) <- "DiffSEObj"
    return(out)
}

##return an SE object
##bootIters should contain ATT as first object in list
#'@title computeSE
#' 
#' @description Computes standard errors from bootstrap results.  This function
#'  is called by several functions in the qte package
#' 
#' @param bootIters List of bootstrap iterations
#' @inheritParams panel.qtet
#'
#' @keywords internal
#'
#' @return SEObj
computeSE <- function(bootIters, qteobj, alp=0.05) {
    ##For now, just plot the qte and att with standard errors
    ##helper function to get the first element out of a list
    qte <- qteobj$qte
    ate <- qteobj$ate
    iters <- length(bootIters)

    getElement <- function(Lst, elemNum) {
        return(as.numeric(unlist((Lst[elemNum])))) #as.numeric is a trick to 
        ##get numerical value of qte
    }
    all.ate = unlist(sapply(bootIters, FUN=getElement,elemNum=2))
    ##get se
    ate.se <- sd(all.ate)
    ate.upper <- ate + qnorm(1-alp/2)*ate.se
    ate.lower <- ate - qnorm(1-alp/2)*ate.se
    ##reorder asc
    ##all.ate = all.ate[order(all.ate)]
    ##ate.upper = all.ate[min(iters,round((1-alp/2)*iters))]
    ##ate.lower = all.ate[max(1,round((alp/2)*iters))]
    
    ##now get CI for qte:
    all.qte = lapply(bootIters, FUN=getElement, elemNum=1)
    qte.mat = do.call(rbind,lapply(all.qte, FUN=as.numeric, ncol=length(all.qte[[1]]), byrow=TRUE))
    ##standard error
    qte.se <- apply(qte.mat, FUN=sd, MARGIN=2)

    sigmahalf <- (apply(qte.mat, 2, function(b) quantile(b, .75, type=1)) -
                  apply(qte.mat, 2, function(b) quantile(b, .25, type=1))) / (qnorm(.75) - qnorm(.25))

    ## this seems to work a bit better in practice when some QTEs are 0
    if (any(sigmahalf==0)) {
        sigmahalf <- apply(qte.mat, 2, sd)
        sigmahalf <- sapply(1:length(qte), function(i) max(sigmahalf[i], .000000001))
    }
    cb <- apply(qte.mat, 1, function(q) max(abs((q-qte)/sigmahalf)))
    c <- quantile(cb, .95, type=1)
    ## qte se by quantiles
    ##sorted.qtemat = apply(qte.mat, 2, sort)
    ##qte.upper = sorted.qtemat[round((1-alp/2)*iters),]
    ##qte.lower = sorted.qtemat[max(1,round((alp/2)*iters)),]
    ## qte se by sd
    qte.upper <- qte + qnorm(1-alp/2)*qte.se
    qte.lower <- qte - qnorm(1-alp/2)*qte.se

    out <- SE(ate.se=ate.se, ate.upper=ate.upper, ate.lower=ate.lower,
              qte.se=qte.se, qte.upper=qte.upper, qte.lower=qte.lower,
              c=c)
    return(out)
}


##summary function for QTE objects
#' @title Summary
#'
#' @description \code{summary.QTE} summarizes QTE objects
#' 
#' @param object A QTE Object
#' @param ... Other params (to work as generic method, but not used)
#' 
#' @export
summary.QTE <- function(object, ...) {
    ##to follow lm, use this function to create a summary.BootQTE object
    ##then that object will have a print method
    ##to check it out for lm, call getAnywhere(print.summary.lm)
    ##and can easily see summary.lm w/o special call
    qte.obj <- object
    
    out <- list(probs=qte.obj$probs, qte=qte.obj$qte,
                       qte.se=qte.obj$qte.se,
                       ate=qte.obj$ate, ate.se=qte.obj$ate.se)
    class(out) <- "summary.QTE"
    return(out)
}

#' @title Print
#' 
#' @description Prints a Summary QTE Object
#' 
#' @param x A summary.QTE object
#' @param ... Other params (required as generic function, but not used)
#' 
#' @export
print.summary.QTE <- function(x, ...) {
    summary.qte.obj <- x
    qte <- summary.qte.obj$qte
    qte.se <- summary.qte.obj$qte.se
    ate <- summary.qte.obj$ate
    ate.se <- summary.qte.obj$ate.se
    probs <- summary.qte.obj$probs
    if (!is.null(qte)) { ##that is we just have att
        ##then, print the qte stuff; otherwise just att stuff
        if (is.null(qte.se)) {
            header <- "QTE"
            body <- qte
        } else {
            header <- c("QTE", "Std. Error")
            body <- cbind(qte, qte.se)
        }
        body <- round(body, digits=3)
        ##colnames(body) <- header
        cat("\n")
        cat("Quantile Treatment Effect:\n")
        cat("\t\t")
        ##cat(header, sep="\t\t")
        cat("\n")
        ##for (i in 1:length(qte)) {
        ##    cat("\t\t")
        ##    cat(format(body[i,],digits=5), sep="\t\t")
        ##    cat("\n")
        ##}
        print.matrix1(body, probs, header=c("tau", header), digits=2, nsmall=2)
        cat("\n")
    }
    cat("Average Treatment Effect:")
    cat("\t")
    cat(format(ate, digits=2, nsmall=2))
    cat("\n")
    if (!is.null(ate.se)) {
        cat("\t Std. Error: \t\t")
        cat(format(ate.se, digits=2, nsmall=2))
        cat("\n")
    }
    ##print(data.frame(body), digits=2)
}

#' @title print.matrix1
#'
#' @description Helper function to print a matrix; used by the print methods
#'
#' @param m Some matrix
#'
#' @keywords internal
print.matrix1 <- function(m, probs=NULL, header=NULL, digits=2, nsmall=2){
    write.table(cbind(probs,
                      format(m, justify="right",
                             digits=digits, nsmall=nsmall)),
                row.names=F, col.names=header, quote=F, sep="\t")
    ##print(m, print.gap=3, right=T)
}

##
#' @title plot.QTE
#' 
#' @description Plots a QTE Object
#' 
#' @param x a QTE Object
#' @param plotate Boolean whether or not to plot the ATE
#' @param plot0 Boolean whether to plot a line at 0
#' @param qtecol Color for qte plot.  Default "black"
#' @param atecol Color for ate plot.  Default "black"
#' @param col0 Color for 0 plot.  Default "black"
#' @param xlab Custom label for x-axis.  Default "tau"
#' @param ylab Custom label for y-axis.  Default "QTE"
#' @param legend Vector of strings to add to legend
#' @param ontreated Boolean whether parameters are "on the treated group"
#' @param ylim The ylim for the plot; if not passed, it will be automatically
#'  set based on the values that the QTE takes
#' @param uselegend Boolean whether or not to print a legend
#' @param legendcol Legend Colors for plotting
#' @param legloc String location for the legend.  Default "topright"
#' @param ... Other parameters to be passed to plot (e.g lwd)
#' 
#' @export
plot.QTE <- function(x, plotate=FALSE, plot0=FALSE,
                         qtecol="black", atecol="black", col0="black",
                         xlab="tau", ylab="QTE",
                         legend=NULL,
                         ontreated=FALSE,
                         ylim=NULL, uselegend=FALSE,
                         legendcol=NULL,
                         legloc="topright", ...) {

    warning("This method is no longer supported.  Try the \"ggqte\" function instead")
    
    qte.obj <- x

    if (is.null(qte.obj$alp)) {
        qte.obj$alp <- .05
    }

    if (!is.null(qte.obj$qte.se)) {
        qte.obj$qte.upper <- qte.obj$qte +
            abs(qnorm(qte.obj$alp/2))*qte.obj$qte.se
        qte.obj$qte.lower <- qte.obj$qte -
            abs(qnorm(qte.obj$alp/2))*qte.obj$qte.se
    }

    if (!is.null(qte.obj$ate.se)) {
        qte.obj$ate.upper <- qte.obj$ate +
            abs(qnorm(qte.obj$alp/2))*qte.obj$ate.se
        qte.obj$ate.lower <- qte.obj$qte -
            abs(qnorm(qte.obj$alp/2))*qte.obj$ate.se
    }

    if (is.null(ylim)) {
        ylim=c(min(qte.obj$qte.lower)-abs(median(qte.obj$qte)),
             max(qte.obj$qte.upper)+abs(median(qte.obj$qte)))
    }
    plot(qte.obj$probs, qte.obj$qte, type="l",
         ylim=ylim,
         xlab=xlab, ylab=ylab, col=qtecol,...)
    lines(qte.obj$probs, qte.obj$qte.lower, lty=2, col=qtecol)
    lines(qte.obj$probs, qte.obj$qte.upper, lty=2, col=qtecol)
    if (plotate) {
        abline(h=qte.obj$ate, col=atecol, ...)
        abline(h=qte.obj$ate.lower, lty=2, col=atecol)
        abline(h=qte.obj$ate.upper, lty=2, col=atecol)
    }
    if (plot0) {
        abline(h=0, col=col0)
    }

    if (uselegend) {
        if (plotate) {
            legend(x=legloc, legend=ifelse(is.null(legend), ifelse(!ontreated, c("QTE", "ATE"), c("QTET", "ATT")), legend), col=ifelse(is.null(legend), c(qtecol, atecol), legendcol), ...)
        } else {
            legend(x=legloc, legend=ifelse(is.null(legend), ifelse(!ontreated, c("QTE"), c("QTET")), legend), col=ifelse(is.null(legend), c(qtecol), legendcol), ...)
        }
    }
}



#####SETUP CLASSES################
#' @title QTE
#'
#' @description Main class of objects.  A \code{QTE} object is returned by
#'  all of the methods that compute the QTE or QTET.
#'
#' @param qte The Quantile Treatment Effect at each value of probs
#' @param qte.se A vector of standard errors for each qte
#' @param qte.upper A vector of upper confidence intervals for each qte (it is
#'  based on the bootstrap confidence interval -- not the se -- so it may not
#'  be symmetric about the qte
#' @param qte.lower A vector of lower confidence intervals for each qte (it is
#'  based on the bootstrap confidence interval -- not the se -- so it may not
#'  be symmyetric about the qte
#' @param ate The Average Treatment Effect (or Average Treatment Effect on
#'  the Treated)
#' @param ate.se The standard error for the ATE
#' @param ate.lower Lower confidence interval for the ATE (it is based on the
#'  bootstrap confidence intervall -- not the se -- so it may not be symmetric
#'  about the ATE
#' @param ate.upper Upper confidence interval for the ATE (it is based on the
#'  bootstrap confidence interval -- not the se -- so it may not be symmetric
#'  about the ATE
#' @param c The critical value from a KS-type statistic used for creating
#'  uniform confidence bands
#' @param pscore.reg The results of propensity score regression, if specified
#' @param probs The values for which the qte is computed
#' @param type Takes the values "On the Treated" or "Population" to indicate
#'  whether the estimated QTE is for the treated group or for the entire
#'  population
#' @param F.treated.t Distribution of treated outcomes for the treated group at
#'  period t
#' @param F.untreated.t Distribution of untreated potential outcomes for the
#'  untreated group at period t
#' @param F.treated.t.cf Counterfactual distribution of untreated potential
#'  outcomes for the treated group at period t
#' @param F.treated.tmin1 Distribution of treated outcomes for the
#'  treated group at period tmin1
#' @param F.treated.tmin2 Distribution of treated outcomes for the
#'  treated group at period tmin2
#' @param F.treated.change.tmin1 Distribution of the change in outcomes for
#'  the treated group between periods tmin1 and tmin2
#' @param F.untreated.change.t Distribution of the change in outcomes for the
#'  untreated group between periods t and tmin1
#' @param F.untreated.change.tmin1 Distribution of the change in outcomes for
#'  the untreated group between periods tmin1 and tmin2
#' @param F.untreated.tmin1 Distribution of outcomes for the untreated group
#'  in period tmin1
#' @param F.untreated.tmin2 Distribution of outcomes for the untreated group
#'  in period tmin2
#' @param condQ.treated.t Conditional quantiles for the treated group in
#'  period t
#' @param condQ.treated.t.cf Counterfactual conditional quantiles for the treated
#'  group in period t
#' @param eachIterList An optional list of the outcome of each bootstrap
#'  iteration
#' @param inffunct The influence function for the treated group;
#'  used for inference when there are multiple
#'  periods and in the case with panel data.  It is needed for computing covariance
#'  terms in the variance-covariance matrix.
#' @param inffuncu The influence function for the untreated group
#'
#' @export
QTE <- function(qte, ate=NULL, qte.se=NULL, qte.lower=NULL,
                qte.upper=NULL, ate.se=NULL, ate.lower=NULL, ate.upper=NULL,
                c=NULL, pscore.reg=NULL, probs, type="On the Treated",
                F.treated.t=NULL, F.untreated.t=NULL, F.treated.t.cf=NULL,
                F.treated.tmin1=NULL, F.treated.tmin2=NULL,
                F.treated.change.tmin1=NULL,
                F.untreated.change.t=NULL,
                F.untreated.change.tmin1=NULL,
                F.untreated.tmin1=NULL,
                F.untreated.tmin2=NULL,
                condQ.treated.t=NULL,
                condQ.treated.t.cf=NULL,
                eachIterList=NULL, inffunct=NULL, inffuncu=NULL) {
    out <- list(qte=qte, ate=ate, qte.se=qte.se, qte.lower=qte.lower,
                qte.upper=qte.upper, ate.se=ate.se, ate.lower=ate.lower,
                ate.upper=ate.upper, c=c,
                pscore.reg=pscore.reg, probs=probs,
                type=type, F.treated.t=F.treated.t, F.untreated.t=F.untreated.t,
                F.treated.t.cf=F.treated.t.cf,
                F.treated.tmin1=F.treated.tmin1,
                F.treated.tmin2=F.treated.tmin2,
                F.treated.change.tmin1=F.treated.change.tmin1,
                F.untreated.change.t=F.untreated.change.t,
                F.untreated.change.tmin1=F.untreated.change.tmin1,
                F.untreated.tmin1=F.untreated.tmin1,
                F.untreated.tmin2=F.untreated.tmin2,
                condQ.treated.t=condQ.treated.t,
                condQ.treated.t.cf=condQ.treated.t.cf,
                eachIterList=eachIterList,
                inffunct=inffunct,
                inffuncu=inffuncu)
    class(out) <- "QTE"
    return(out)
}

#' @title SE
#'
#' @description Class for Standard Error Objects
#'
#' @param qte.se The QTE Standard Error
#' @param ate.se The ATE Standard Error
#' @param qte.upper The QTE upper CI
#' @param qte.lower The QTE lower CI
#' @param ate.upper The ATE upper CI
#' @param ate.lower The ATE lower CI
#' @param c The critical value from a KS-type statistic used for creating
#'  uniform confidence bands
#' @param probs The values at which the QTE is computed
#'
#' @keywords internal
SE <- function(qte.se=NULL, ate.se=NULL, qte.upper=NULL, qte.lower=NULL,
               ate.upper=NULL, ate.lower=NULL, c=NULL, probs=NULL) {

    out <- list(qte.se=qte.se, qte.upper=qte.upper, qte.lower=qte.lower,
                ate.se=ate.se, ate.upper=ate.upper, ate.lower=ate.lower,
                c=c, probs=probs)
    class(out) <- "SE"
    return(out)
}



############## DATA DOCUMENTATION ################
#' @title Lalonde (1986)'s NSW Dataset
#' 
#' @description \code{lalonde} contains data from the National Supported Work
#'  Demonstration.  This program randomly assigned applicants to the job
#'  training program (or out of the job training program).  The dataset is
#'  discussed in Lalonde (1986).  The experimental part of the dataset is
#'  combined with an observational dataset from the Panel Study of Income
#'  Dynamics (PSID).  Lalonde (1986) and many subsequent papers (e.g.
#'  Heckman and Hotz (1989), Dehejia and Wahba (1999), Smith and Todd (2005),
#'  and Firpo (2007) have used this combination to study the effectiveness
#'  of various `observational' methods (e.g. regression, Heckman selection,
#'  Difference in Differences, and propensity score matching) of estimating
#'  the Average Treatment Effect (ATE) of participating in the job training
#'  program.  The idea is that the results from the observational method
#'  can be compared to results that can be easily obtained from the
#'  experimental portion of the dataset.
#'
#'  To be clear, the observational data combines the observations that are
#'  treated from the experimental portion of the data with untreated observations
#'  from the PSID.
#' 
#' @format Four data.frames: (i) lalonde.exp contains a cross sectional version
#'  of the experimental data, (ii) lalonde.psid contains a cross sectional
#'  version of the observational data, (iii) lalonde.exp.panel contains a
#'  panel version of the experimental data, and (iv) lalonde.psid.panel contains
#'  a panel version of the observational data.  Note: the cross sectional
#'  and panel versions of each dataset are identical up to their shape; in
#'  demonstrating each of the methods, it is sometimes convenient to have
#'  one form of the data or the other.
#' @docType data
#' @name lalonde
#' @usage data(lalonde)
#' @references LaLonde, Robert.  ``Evaluating the Econometric Evaluations of
#'  Training Programs with Experimental Data.'' The American Economics Review,
#'  pp. 604-620, 1986.
#'  @source The dataset comes from Lalonde (1986) and has been studied in much
#'  subsequent work.  The \code{qte} package uses a version from the
#'  \code{causalsens} package
#'  (\url{https://CRAN.R-project.org/package=causalsens})
#' @keywords datasets
NULL

#' @title Lalonde's Experimental Dataset
#'
#' @description The cross sectional verion of the experimental part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.exp
#' @keywords datasets
NULL

#' @title Lalonde's Panel Experimental Dataset
#'
#' @description The panel verion of the experimental part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.exp.panel
#' @keywords datasets
NULL

#' @title Lalonde's Observational Dataset
#'
#' @description The cross sectional verion of the observational part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.psid
#' @keywords datasets
NULL

#' @title Lalonde's Experimental Dataset
#'
#' @description The panel verion of the observational part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.psid.panel
#' @keywords datasets
NULL
