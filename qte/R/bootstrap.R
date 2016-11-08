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
