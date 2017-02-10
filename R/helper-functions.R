###makeBalancedPanel is a function to take a dataset
## and make sure that all years are available for 
## all observations.  If some years are not available,
## then that observation is dropped.
#'@title makeBalancedPanel
#' 
#' @description This function drops observations from data.frame
#'  that are not part of balanced panel data set.
#' 
#' @param data data.frame used in function
#' @param idname unique id
#' @param tname time period name
#' 
#' @keywords internal
makeBalancedPanel <- function(data, idname, tname) {
    ## this is the old way
    ##
    ## data=droplevels(data)
    ## allt = unique(data[,tname])
    ## allid = unique(data[,idname])
    
    ## ##loop over each id in the dataset
    ## for (id in allid) {
    ##     ##get the number of time periods for that id
    ##     this.allt = unique(data[data[,idname]==id,tname])
        
    ##     ##check if its equal to the largest number of time
    ##     ##periods in the dataset
    ##     if (!(length(this.allt) == length(allt))) {
    ##         ##if it is fewer, then drop all observations
    ##         ##from that id from the dataset
    ##         data = data[!(data[,idname] == id),]
    ##     }
    ## }
    ## return(data)

    nt <- length(unique(data[,tname]))
    agg <- aggregate(data[,idname], by=list(data[,idname]), length)
    rightids <- agg[,1][agg[,2]==nt]
    bp <- data[data[,idname] %in% rightids,]
    return(bp)

}


#'@title panel2cs
#' 
#' @description panel2cs takes a 2 period dataset and turns it
#'  into a cross sectional dataset.  The default functionality
#'  is to keep all the variables from period 1
#'  and add all the variables listed by name in timevars
#'  from period 2 to those
#' 
#' @param data data.frame used in function
#' @param timevars vector of names of variables to keep
#' @param idname unique id
#' @param tname time period name
#' 
#' @keywords internal
panel2cs <- function(data, timevars, idname, tname) {

    if (length(unique(data[,tname])) != 2) {
        stop("panel2cs only for 2 periods of panel data")
    }

    data <- makeBalancedPanel(data, idname, tname) ## just in case
    data <- data[order(data[,idname], data[,tname]),] ## put everything in the right order
    ## so we can match it easily later on

    tdta <- aggregate(data[,timevars], by=list(data[,idname]), FUN=function(x) { x[2] })

    t1 <- unique(data[,tname])
    t1 <- t1[order(t1)][1]
    retdat <- subset(data, data[,tname]==t1)
    retdat$yt1 <- tdta[,2]
    retdat$dy <- retdat$yt1 - retdat$y
    return(retdat)
   
}

#'@title ids2rownum
#' 
#' @description ids2rownum takes a vector of ids and converts it t the right
#'  row number in the dataset; ids should be unique in the dataset
#'  that is, don't pass the function panel data with multiple same ids
#' 
#' @param ids vector of ids
#' @param data data frame
#' @param idname unique id
#' 
#' @keywords internal
ids2rownum <- function(ids, data, idname) {
    vapply(ids, id2rownum, 1.0, data=data, idname=idname)
}


#'@title ids2rownum
#' 
#' @description id2rownum takes an id and converts it t the right
#'  row number in the dataset; ids should be unique in the dataset
#'  that is, don't pass the function panel data with multiple same ids
#' 
#' @param ids vector of ids
#' @param data data frame
#' @param idname unique id
#' 
#' @keywords internal
id2rownum <- function(id, data, idname) {
    which(data[,idname] == id)
}




#'@title makeDist
#' 
#' @description turn vectors of a values and their distribution function values
#'  into an ecdf.  Vectors should be the same length and both increasing.
#' 
#' @param x vector of values
#' @param Fx vector of the distribution function values
#' @param sorted boolean indicating whether or not x is already sorted;
#'  computation is somewhat faster if already sorted
#' 
#' @keywords internal
makeDist <- function(x, Fx, sorted=FALSE) {
    if (!sorted) {
        tmat <- cbind(x, Fx)
        tmat <- tmat[order(x),]
        x <- tmat[,1]
        Fx <- tmat[,2]
    }
    
    retF <- approxfun(x, Fx, method="constant",
                      yleft=0, yright=1, f=0, ties="ordered")
    class(retF) <- c("ecdf", "stepfun", class(retF))
    assign("nobs", length(x), envir = environment(retF))
    retF
}


#'@title checkfun
#' 
#' @description The check function used for optimizing to get quantiles
#' 
#' @param a vector to compute quantiles for
#' @param tau between 0 and 1, ex. .5 implies get the median
#' 
#' @keywords internal
checkfun <- function(a, tau) {
    return(a*(tau - (1*(a<=0))))
}

#'@title weighted.checkfun
#' 
#' @description Weights the check function
#' 
#' @param q the value to check
#' @param tau between 0 and 1, ex. .5 implies get the median
#' @param weights the weights, weighted.checkfun normalizes the weights
#'  to sum to 1.
#' 
#' @keywords internal
weighted.checkfun = function(q, cvec, tau, weights) {
    w <- weights/(sum(weights))
    retval = sum(w*checkfun(cvec-q, tau))
    return(retval)
}


#'@title getWeightedQuantile
#' 
#' @description Finds the quantile by optimizing the weighted check function
#' 
#' @param tau between 0 and 1, ex. .5 implies get the median
#' @param cvec a vector to compute quantiles for
#' @param weights the weights, weighted.checkfun normalizes the weights
#'  to sum to 1.
#' @param norm normalize the weights so that they have mean of 1, default is
#'  to normalize
#' 
#' @keywords internal
getWeightedQuantile <- function(tau, cvec, weights=NULL, norm=TRUE) {
    if (is.null(weights)) {
        weights <- 1
    }
    mw <- mean(weights)
    if (norm) {
        weights <- weights / mw
    }
    return(optimize(weighted.checkfun, 
                    lower=min(cvec),
                    upper=max(cvec),
                    cvec=cvec, tau=tau, weights=weights)$minimum)
}

#'@title getWeightedQuantile
#' 
#' @description Finds multiple quantiles by repeatedly calling
#'  getWeightedQuantile
#' 
#' @param tau a vector of values between 0 and 1
#' @param cvec a vector to compute quantiles for
#' @param weights the weights, weighted.checkfun normalizes the weights
#'  to sum to 1.
#' @param norm normalize the weights so that they have mean of 1, default is
#'  to normalize
#' 
getWeightedQuantiles <- function(tau, cvec, weights=NULL, norm=TRUE) {
    vapply(tau, getWeightedQuantile, 1.0, cvec=cvec, weights=weights, norm=norm)
    ##wtd.quantile(cvec, weights=weights, probs=tau, normwt=T)
}

#'@title getWeightedMean
#' 
#' @description Get the mean applying some weights
#' 
#' @param y a vector to compute the mean for
#' @param weights the vector of weights, can be NULL, then will just return mean
#' @param norm normalize the weights so that they have mean of 1, default is
#'  to normalize
#' 
#' @keywords internal
getWeightedMean <- function(y, weights=NULL, norm=TRUE) {
    if (is.null(weights)) {
        weights <- 1
    }
    mw <- mean(weights)
    if (norm) {
        weights <- weights/mw
    }
    mean(weights*y)
}

#'@title getWeightedDf
#' 
#' @description Get the mean applying some weights
#' 
#' @param y a vector to compute the mean for
#' @param y.seq an optional vector of values to compute the distribution function
#'  for; the default is to use all unique values of y
#' @param weights the vector of weights, can be NULL, then will just return mean
#' @param norm normalize the weights so that they have mean of 1, default is
#'  to normalize
#' 
#' @keywords internal
getWeightedDf <- function(y, y.seq=NULL, weights=NULL, norm=TRUE) {
    if (is.null(weights)) {
        weights <- 1
    }
    mw <- mean(weights)
    if (norm) {
        weights <- weights/mw
    }
    if (is.null(y.seq)) {
        y.seq <- unique(y)
        y.seq <- y.seq[order(y.seq)]
    }
    dvals <- vapply(y.seq, FUN=function(x) { mean(weights*(y <= x)) }, 1.0)
    makeDist(y.seq, dvals)
}


#'@title qte2mat
#'
#' @description Turn multiple qtes into a matrix for printing
#'
#' @param qteList a list of qte objects
#' @param sset subset of qtes to keep
#' @param se whether or not to include standard errors in the resulting matrix
#' @param rnd how many disgits to round to
#'
#' @return matrix
qtes2mat <- function(qteList, sset=NULL, se=TRUE, rnd=3) {
    if (is.null(sset)) {
        sset <- seq(1,length(qteList[[1]]$probs))
    }
    probs <- qteList[[1]]$probs[sset]
    outmat <- matrix(nrow=2*length(qteList[[1]]$qte[sset]), ncol=length(qteList))
    lqte <- lapply(qteList, function(x) { paste(round(x$qte[sset], rnd)) })
    lqte.se <- lapply(qteList, function(x) { paste("(",round(x$qte.se[sset],rnd),")", sep="") })
    for (i in 1:length(qteList)) {
        outmat[,i] <- c(rbind(lqte[[i]], lqte.se[[i]]))
    }
    probsvals <- c(rbind(paste(probs), rep("", length(probs))))
    outmat <- cbind(probsvals, outmat)
    outmat
}


#'@title cs2panel
#'
#' @description Turn repeated cross sections data into panel data by imposing rank invariance; does not
#'  that the inputs have the same length
#'
#' @param cs1 data frame, the first cross section
#' @param cs2 data frame, the second cross section
#' @param yname the name of the variable to calculate difference for (should be the same in each dataset)
#' 
#' @return the change in outcomes over time
cs2panel <- function(cs1, cs2, yname) {
    nu <- min(nrow(cs2), nrow(cs2))
    if (nu == nrow(cs2)) {
        ut <- cs2[,yname]
        ut <- ut[order(-ut)] ##orders largest to smallest
        ps <- seq(1,0,length.out=length(ut)) ##orders largest to smallest
        utmin1 <- quantile(cs1[,yname], probs=ps, type=1)
        ##F.untreated.change.t <- ecdf(ut-utmin1)
    } else {
        utmin1 <- cs2[,yname]
        utmin1 <- utmin1[order(-utmin1)] ##orders largest to smallest
        ps <- seq(1,0,length.out=length(utmin1)) ##orders largest to smallest
        ut <- quantile(cs1[,yname], probs=ps, type=1)
        ##F.untreated.change.t <- ecdf(ut-utmin1)
    }
    return(ut - utmin1)
}

require(formula.tools)
    
dr.inner <- function(yval, formla, data) {
    y <- lhs(formla)
    x <- rhs(formla)
    lhs(formla) <- substitute(I(y <= yval), list(y=y, yval=yval))
    outreg <- glm(formla, data, family=binomial(link=logit))
}

#'@title dr
#'
#' @description Distribution Regression
#' 
#' @param formla the regression to run
#' @param y.seq the values of y to run the regression on
#' 
#' @return 
dr <- function(formla, data, y.seq) {
    distreg <- lapply(y.seq, dr.inner, formla, data)
    DR(y.seq, distreg)    
}

##function to take in y0 and x and return F(y0|x)
dr.predict.inner <- function(y0, x, drobj) {
    yval <- drobj$yvals[which.min(abs(drobj$yvals-y0))]
    yidx <- which(drobj$yvals==yval)
    predict(drobj$drlist[[yidx]], newdata=x, type="response")
}

##function to take a vector of ys and single x and return vector F(y|x)
dr.predict <- function(y.seq, x, drobj) {
    vapply(y.seq, dr.predict.inner, FUN.VALUE=1.0, x=x, drobj=drobj)
}




