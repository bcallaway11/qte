####Bounds with Fan-yu
##Function to implement bounds using the method of Fan and Yu (2012)
#' @title bounds
#' @description \code{bounds} estimates bounds for the Quantile Treatment
#'  Effect on the
#'  Treated (QTET) using the method of Fan and Yu (2012).
#' @inheritParams panel.qtet
#'
#' @examples
#' ## load the data
#' data(lalonde)
#'
#' ## Run the bounds method with no covariates
#' b1 <- bounds(re ~ treat, t=1978, tmin1=1975, data=lalonde.psid.panel,
#'   idname="id", tname="year")
#' summary(b1)
#'
#' @references
#' Fan, Yanqin and Zhengfei Yu.  ``Partial Identification of Distributional
#'  and Quantile Treatment Effects in Difference-in-Differences Models.''
#'  Economics Letters 115.3, pp.511-515, 2012.
#'
#' @return A \code{BoundsObj} object
#' 
#' @export
bounds <- function(formla, xformla=NULL, t, tmin1, tname, data,
                   idname, 
                   probs=seq(0.05,0.95,0.05)) {
    form = stats::as.formula(formla)
    dta = stats::model.frame(stats::terms(form,data=data),data=data) #or model.matrix
    colnames(dta) = c("y","treatment")
    yname="y"
    treat="treatment"
    data=cbind.data.frame(dta,data)

    ## Setup x variables if using formula
    if (!(is.null(xformla))) {
        x <- colnames(stats::model.matrix(terms(stats::as.formula(xformla)), data=data))
        data <- cbind(data[,c(yname,treat,idname,tname)],
                      stats::model.matrix(terms(stats::as.formula(xformla)), data=data))
    }
    
    ##drop the always treated.  Note that this also relies
    ##on no "switchback" or no return to untreated status
    ##after joining treatment.
    
    data = subset(data, (data[,tname]==tmin1 | data[,tname]==t))
    data = makeBalancedPanel(data, idname, tname)
        
    ##just to make sure the factors are working ok
    data = droplevels(data)
    
    ##a) get all the treated (in the last period) observations
    treated.t = data[data[,tname]==t & data[,treat]==1,]

    ##b) set ever.treated to 1 if observation is treated in last period
    ##data$ever.treated = data$treatment
    ##data$ever.treated = 1*(data[,idname] %in% treated.t[,idname])  
    ##ever.treated = "ever.treated"

    ##Setup each of the datasets used below
    ##treated.t = subset(data, data[,treat]==1 & data[,tname]==t)
    ##just get the lagged value of y; otherwise keep the same
    ##dataset.  Note: this will not work if there are x covariates;
    ##well, could follow similar procedure, but as is, would
    ##require some modification.
    treated.tmin1 = data[ data[,tname] == tmin1 & 
        data[,treat] == 1, ]
    ##this is right
    ##untreated.t = subset(data, data[,treat]==0 & data[,tname]==t)
    untreated.t = data[data[,tname]==t & data[,treat]==0,]
    ##get lagged of untreated y
    untreated.tmin1 = data[ data[,tname] == tmin1 &
        data[,treat] == 0, ]
    
    
    ##First, get distribution Y_1t | Dt=1
    F.treated.t = stats::ecdf(treated.t[,yname])
    

    F.treated.tmin1 = stats::ecdf(treated.tmin1[,yname]) #as long as 

    F.treated.change = stats::ecdf(treated.t[,yname]-treated.tmin1[,yname])
    ##Actually -- don't think you need that...
    
    ##2c) Get the distribution of the change in outcomes for the never treated
    F.untreated.change.t = stats::ecdf(untreated.t[,yname]-untreated.tmin1[,yname])

    ##for comparison, compute att first
    att = mean(treated.t[,yname]) - mean(treated.tmin1[,yname]) -
        (mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname]))
    ##2c.1) If there are covariates, then above distribution needs to be changed
    if (!(is.null(xformla))) {
        ##set up the data to do the propensity score re-weighting
        ##we need to bind the datasets back together to estimate pscore
        treated.t$changey = treated.t[,yname] - treated.tmin1[,yname]
        untreated.t$changey = untreated.t[,yname] - untreated.tmin1[,yname]
        pscore.data = rbind(treated.t, untreated.t)
        xmat = pscore.data[,x]
        pscore.reg = stats::glm(pscore.data[,treat] ~ as.matrix(xmat),
            family=stats::binomial(link="logit"))
        pscore = stats::fitted(pscore.reg)
        pscore.data$pscore <- pscore
        pD1 = nrow(treated.t)/nrow(untreated.t)
        pval <- pD1

        ##this contains the support of the change in y
        p.dy.seq = pscore.data$changey #unique(pscore.data$changey)
        ##TODO: What is this?  Need to come up with better name for this variable
        distvals = rep(0,length(p.dy.seq))
        for (dy in p.dy.seq) {
            distvals[which(dy==p.dy.seq)] = mean(1*(pscore.data$changey<=dy)*
                        (1-pscore.data[,treat])*pscore/((1-pscore)*pD1))
        }
        pscore.data$distvals = distvals
        
        pscore.data1 = pscore.data[order(pscore.data$changey),]

        F.untreated.change.t = stats::approxfun(pscore.data1$changey,
            pscore.data1$distvals, method="constant",
            yleft=0, yright=1, f=0, ties="ordered")
        class(F.untreated.change.t) = c("ecdf", "stepfun",
                 class(F.untreated.change.t))
        assign("nobs", length(p.dy.seq), envir = environment(F.untreated.change.t))
        ##att using abadie method
        att = mean(((pscore.data$changey)/pD1)*(pscore.data[,treat] - pscore)/(1-pscore))
    }
    
    ##2c) Get the distribution of outcomes for the newly treated at (t-1)
    ##F.newlytreated.tmin1 <<- ecdf(newly.treated.tmin1[,yname])
    
    ##2d) get the lower bound
    ##make sure that this is right, but we are taking the smallest
    ## over the support (I think) of y
    supy = sort(unique(c(untreated.t$y, treated.tmin1$y, untreated.tmin1$y)))#this should have largest support
    posvals = sort(unique(untreated.t$y - untreated.tmin1$y)) #these are the values to min over; not sure
    ##exactly what they should be, but should cover 0 probably about
    ##as wide as the support of y is in each direction
    ## and should probably be passed into the function
    ## I think that I can pass this as both arguments s, and y below
    ## but maybe should separate them esp. if there are issues
    lbs = vapply(supy,FUN=getlb,FUN.VALUE=1, 
        F.change.treated=F.untreated.change.t,
        F.treated.tmin1=F.treated.tmin1,
        y=posvals)
    F.lb <- stats::approxfun(supy, lbs, method="constant",
                      yleft=0, yright=1, f=0, ties="ordered")
    class(F.lb) = c("ecdf", "stepfun", class(F.lb))
    assign("nobs", length(supy), envir = environment(F.lb))

    ubs = vapply(supy,FUN=getub,FUN.VALUE=1, 
        F.change.treated=F.untreated.change.t,
        F.treated.tmin1=F.treated.tmin1,
        y=posvals)
    F.ub <- stats::approxfun(supy, ubs, method="constant",
                      yleft=0, yright=1, f=0, ties="ordered")
    class(F.ub) = c("ecdf", "stepfun", class(F.ub))
    assign("nobs", length(supy), envir = environment(F.ub))

    
    ##get upper bound quantiles for unobserved untreated observations
    ##these are opposite from lower bound / upper bound on distribution
    ub.quantiles = stats::quantile(F.lb, probs=probs)
    
    ##get lower bound quantiles for unobserved untreated observations
    lb.quantiles = stats::quantile(F.ub, probs=probs)
    
    ##plot bounds on qte
    ## because we are subtracting, the lower bound for the qte
    ## will occur at the upper bound of the quantiles of untreated
    ## distribution, and the upper bound will occur at the lower
    ## bound of the quantiles of the untreated distribution.
    lb.qte = as.numeric(stats::quantile(treated.t[,yname],probs=probs) - 
        ub.quantiles)
    ub.qte = as.numeric(stats::quantile(treated.t[,yname],probs=probs) - 
        lb.quantiles)
    ## if (plot) {
    ##     plot(probs, lb.qte, 
    ##          type="l", lty=2, xlab="tau", ylab="QTE",
    ##          ylim=c(-2.5,2.5))
    ##     graphics::lines(probs, ub.qte, lty=2)
    ##     graphics::abline(a=att, b=0, col="blue")
    ## }    
    return(BoundsObj(lbs=lbs,ubs=ubs, ub.quantiles=ub.quantiles,
                lb.quantiles=lb.quantiles, ub.qte=ub.qte,
                lb.qte = lb.qte, att=att, probs=probs))
}

##summary function for bounds objects
#' @title Summary of BoundsObj
#'
#' @description \code{summary.BoundsObj} is an object for holding
#' \code{bounds} results
#' 
#' @param object A BoundsObj Object
#' @param ... Other params (for consistency as generic S3 method, but not used)
#'
#' @return summary.BoundsObj Object
#' 
#' @export
summary.BoundsObj <- function(object, ...) {
    ##to follow lm, use this function to create a summary.BootQTE object
    ##then that object will have a print method
    ##to check it out for lm, call getAnywhere(print.summary.lm)
    ##and can easily see summary.lm w/o special call
    bounds.obj <- object
    
    out <- list(lbs=bounds.obj$lbs, ubs=bounds.obj$ubs,
                lb.quantiles=bounds.obj$lb.quantiles,
                ub.quantiles=bounds.obj$ub.quantiles,
                lb.qte=bounds.obj$lb.qte,
                ub.qte=bounds.obj$ub.qte,
                att=bounds.obj$att,
                probs=bounds.obj$probs)
    class(out) <- "summary.BoundsObj"
    return(out)
}

#' @title Print a summary.BoundsObj
#' 
#' @description Prints a Summary QTE Object
#' 
#' @param x A summary.BoundsObj
#' @param ... Other objects to pass (not used)
#' 
#' @export
print.summary.BoundsObj <- function(x, ...) {
    summary.bounds.obj <- x
    
    lb.qte <- summary.bounds.obj$lb.qte
    ub.qte <- summary.bounds.obj$ub.qte
    probs <- summary.bounds.obj$probs
    att <- summary.bounds.obj$att
    header <- c("tau", "Lower Bound", "Upper Bound")
    body <- cbind(as.numeric(probs), lb.qte, ub.qte)
    body <- round(body, digits=2)
    colnames(body) <- header
    cat("\n")
    cat("Bounds on the Quantile Treatment Effect on the Treated:\n")
    cat("\t\t")
    ##cat(header, sep="\t\t")
    cat("\n")
    ##for (i in 1:length(qte)) {
    ##    cat("\t\t")
    ##    cat(format(body[i,],digits=5), sep="\t\t")
    ##    cat("\n")
    ##}
    print.matrix2(rbind(header, body), header=header)
    cat("\n")
    cat("Average Treatment Effect on the Treated:")
    cat("\t")
    cat(format(att, digits=3, nsmall=2))
    cat("\n")
    ##print(data.frame(body), digits=2)
}

##
#' @title print.matrix2
#'
#' @description Helper function to print a matrix; used by the print methods
#'
#' @param m Some matrix
#'
#' @keywords internal
print.matrix2 <- function(m, header=NULL, digits=2, nsmall=2){
    write.table(format(m, justify="right",
                       digits=digits, nsmall=nsmall),
                row.names=F, col.names=header, quote=F, sep="\t")
    ##print(m, print.gap=3, right=T)
}


##
#' @title Plot Bounds
#' 
#' @description Plots a BoundObj Object
#'
#' @inheritParams plot.QTE
#' @param x A BoundsObj Object
#' 
#' @export 
plot.BoundsObj <- function(x, plotate=FALSE, plot0=FALSE,
                         qtecol="black", atecol="black", col0="black",
                         ylim=NULL, uselegend=FALSE, legloc="topright", ...) {
    bounds.obj <- x

    if (is.null(ylim)) {
        ylim=c(min(bounds.obj$lb.qte)-abs(median(bounds.obj$lb.qte)),
             max(bounds.obj$ub.qte)+abs(median(bounds.obj$ub.qte)))
    }
    plot(bounds.obj$probs, bounds.obj$lb.qte, type="l",
         ylim=ylim,
         xlab="tau", ylab="QTET", col=qtecol,...)
    graphics::lines(bounds.obj$probs, bounds.obj$ub.qte, col=qtecol)
    if (plotate) {
        graphics::abline(h=bounds.obj$att, col=atecol, ...)
    }
    if (plot0) {
        graphics::abline(h=0, col=col0)
    }

    if (uselegend) {
        if (plotate) {
            legend(x=legloc, legend=c("QTET Bounds", "ATT"), col=c(qtecol, atecol), ...)
        } else {
            legend(x=legloc, legend=c("QTET Bounds"), col=c(qtecol), ...)
        }
    }
}





#####HELPER FUNCTIONS FOR FAN-YU######

#' @title getlb
#'
#' @description Helper function to compute the lower bound in bounds method.
#'  Usually called by vapply function.
#' @param s A particular value of distribution for which to calculate the bound
#' @param F.change.treated ecdf object of distribution of change in outcomes
#'  for the treated group
#' @param F.treated.tmin1 ecdf object of distribution of outcomes in period
#'  t-1 for the treated group
#' @param y a vector of values that observations could take in the previous
#'  period ? 
#' @keywords internal
getlb <- function(s, F.change.treated, F.treated.tmin1, y) {
    return(max(F.change.treated(y) + F.treated.tmin1(s-y) - 1,0))
}

#' @title getub
#'
#' @description Helper function to compute the upper bound in bounds method.
#'  It is usually called by vapply function
#' @inheritParams getlb
#' 
#' @keywords internal
getub <- function(s, F.change.treated, F.treated.tmin1, y) {
    return(1 + min(F.change.treated(y) + F.treated.tmin1(s-y) - 1,0))
}

#' @title BoundsObj
#'
#' @description An object of results from computing bounds
#'
#' @param lbs A vector of the lower bounds for each value in the support
#'  of the outcome
#' @param ubs A vector of the upper bounds for each value in the support
#'  of the outcome
#' @param ub.quantiles A vector of the same length as probs that contains
#'  the upper bound of the quantiles of the counterfactual distribution
#'  of untreated potential outcomes for the treated group
#' @param lb.quantiles A vector of the same length as probs that contains
#'  the lower bound of the quantiles of the counterfactual distribution
#'  of untreated potential outcomes for the treated group
#' @param ub.qte The point estimate of the upper bound for the QTE
#' @param lb.qte The point estimate of the lower bound for the QTE
#' @param att The ATT is point identified under the assumptions required
#'  by the bounds method
#' @inheritParams panel.qtet
#'
#' @keywords internal
BoundsObj <- function(lbs, ubs, ub.quantiles, lb.quantiles, ub.qte,
                      lb.qte, att=NULL, probs) {

    out <- list(lbs=lbs,ubs=ubs, ub.quantiles=ub.quantiles,
                lb.quantiles=lb.quantiles, ub.qte=ub.qte,
                lb.qte = lb.qte, att=att, probs=probs)
    class(out) <- "BoundsObj"
    return(out)
}

