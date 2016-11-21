#' @title bootiter
#'
#' @description \code{bootiter} that computes a single bootstrap iteration
#'
#' @inheritParams bootstrap
#' @param i the bootstrap iteration, this is not used except for calling
#'  by apply-type methods
#' 
#'
#' @keywords internal
#' 
#' @return SE object
bootiter <- function(i, qteparams, func) {
    seedvec <- qteparams$seedvec
    data <- qteparams$data
    n <- nrow(data)
    if(!is.null(seedvec)) {
        set.seed(seedvec[i])
    }
    
    boot.data = data[0,]
    randy = sample(1:n, n, replace=T)
    boot.data <- data[randy,]
    newqp <- qteparams
    newqp$data <- boot.data
    thisIter = func(newqp)
    thisIter$data <- NULL ##drop so that the size doesn't get too big
    return(thisIter)
}



#' @title bootstrap
#'
#' @description \code{bootstrap} is a method that bootstraps standard
#'  errors for methods in the qte package
#' 
#' @param qteparams QTEparams object containing the parameters passed to ciqte
#' @param qteest an initial estimate of the qte
#' @param func the internal function that computes estimates of the qte,
#'  should be the same as what was used to produce qteest
#' 
#' @import parallel
#'
#' @keywords internal
#' 
#' @return SE object
bootstrap <- function(qteparams, qteest, func) {
    ##pass the function to call, qteparams, whether to run in parallel

    env <- parent.frame()
    alp <- qteparams$alp
    iters <- qteparams$iters
    retEachIter <- qteparams$retEachIter
    pl <- qteparams$pl
    cores <- qteparams$cores

    eachIter <- list()

    if (pl) {
        eachIter <- mclapply(1:iters, bootiter, qteparams=qteparams, func=func,
                             mc.cores=cores)
    } else {
        eachIter <- lapply(1:iters, bootiter, qteparams=qteparams, func=func)
    }

    SEobj <- computeSE(eachIter, qteest, alp=alp)

    if(!retEachIter) {
        eachIter <- NULL
    }

    assign("eachIter", eachIter, envir=env) ##this seems kind of hackish

    return(SEobj)
}
