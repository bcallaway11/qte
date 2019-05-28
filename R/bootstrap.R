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
    tname <- qteparams$tname
    t <- qteparams$t
    tmin1 <- qteparams$tmin1
    tmin2 <- qteparams$tmin2
    idname <- qteparams$idname
    data <- qteparams$data
    
    if(!is.null(seedvec)) {
        set.seed(seedvec[i])
    }

    if (qteparams$panel) {
        ids <- sample(unique(data[,idname]), length(unique(data[,idname])), replace=T)
        dta1 <- data[data[,tname] == t, ]
        dta2 <- data[data[,tname] == tmin1, ]
        n <- nrow(dta1)
        if (nrow(dta1) != nrow(dta2)) {
            warning("unexpected unbalanced panel")
        }
        boot.dta1 <- dta1[ids2rownum(ids, dta1, idname),]
        boot.dta2 <- dta2[ids2rownum(ids, dta2, idname),] ##c(ids, n+ids),]
        newids <- sample(seq(1, n), n)
        boot.dta1[,idname] <- newids
        boot.dta2[,idname] <- newids
        boot.dta <- rbind(boot.dta1, boot.dta2)
        if (!is.null(tmin2)) {
            dta3 <- data[data[,tname] == tmin2, ]
            boot.dta3 <- dta3[ids2rownum(ids, dta3, idname),]
            boot.dta3[,idname] <- newids
            boot.dta <- rbind(boot.dta, boot.dta3)
        }

    } else { ## we are not in the panel case
    
        boot.dta = data[0,]
        randy = sample(1:n, n, replace=T)
        boot.dta <- data[randy,]
    
    }

    newqp <- qteparams
    newqp$data <- boot.dta
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

    eachIter <- pblapply(1:iters, bootiter, qteparams=qteparams,
                         func=func, cl=cores)
    ## if (pl) {
    ##     eachIter <- mclapply(1:iters, bootiter, qteparams=qteparams, func=func,
    ##                          mc.cores=cores)
    ## } else {
    ##     eachIter <- lapply(1:iters, bootiter, qteparams=qteparams, func=func)
    ## }

    SEobj <- computeSE(eachIter, qteest, alp=alp)

    if(!retEachIter) {
        eachIter <- NULL
    }

    assign("eachIter", eachIter, envir=env) ##this seems kind of hackish

    return(SEobj)
}
