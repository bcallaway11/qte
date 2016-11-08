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
#' 
#' @keywords internal
getWeightedQuantile <- function(tau, cvec, weights) {
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
#' 
getWeightedQuantiles <- function(tau, cvec, weights) {
    vapply(tau, getWeightedQuantile, 1.0, cvec=cvec, weights=weights)
}
