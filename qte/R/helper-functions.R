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
    data=droplevels(data)
    allt = unique(data[,tname])
    allid = unique(data[,idname])
    
    ##loop over each id in the dataset
    for (id in allid) {
        ##get the number of time periods for that id
        this.allt = unique(data[data[,idname]==id,tname])
        
        ##check if its equal to the largest number of time
        ##periods in the dataset
        if (!(length(this.allt) == length(allt))) {
            ##if it is fewer, then drop all observations
            ##from that id from the dataset
            data = data[!(data[,idname] == id),]
        }
    }
    return(data)
}

#'@title makeDist
#' 
#' @description turn vectors of a values and their distribution function values
#'  into an ecdf.  Vectors should be the same length and both increasing.
#' 
#' @param x vector of values
#' @param Fx vector of the distribution function values
#' 
#' @keywords internal
makeDist <- function(x, Fx) {
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
    
