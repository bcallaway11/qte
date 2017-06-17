#'@title qtes2mat
#'
#' @description Turn multiple qtes into a matrix for printing
#'
#' @param qteList a list of qte objects
#' @param sset subset of qtes to keep
#' @param se whether or not to include standard errors in the resulting matrix
#' @param rnd how many disgits to round to
#'
#' @return matrix
#' @export
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


#'@title ggqte
#'
#' @description Makes somewhat nicer plots of quantile treatment effects
#'  by using ggplot
#'
#' @import ggplot2
#'
#' @param qteobj a QTE object
#' @param main optional title
#' @param ylab optional y axis label
#' @param ylim optional limits of y axis
#' @param ybreaks optional breaks in y axis
#' @param xbreaks optional breaks in x axis
#' 
#' @return a ggplot object
#' @export
ggqte <- function(qteobj, main="", ylab="", ylim=NULL, ybreaks=NULL, xbreaks=c(.1,.3,.5,.7,.9)) {
    tau <- qteobj$probs
    qte <- qteobj$qte
    qte.se <- qteobj$qte.se
    cmat <- data.frame(tau, qte=qteobj$qte, qte.se=qteobj$qte.se)
    qp <- ggplot2::ggplot(data=cmat, aes(tau, qte, ymax=qte+1.96*qte.se,
                                ymin=qte-1.96*qte.se)) +
        ggplot2::geom_line(aes(tau, qte)) +
        ggplot2::geom_line(aes(tau, qte+1.96*qte.se), linetype="dashed") +
        ggplot2::geom_line(aes(tau, qte-1.96*qte.se), linetype="dashed") +
        ##geom_errorbar(size=.3, width=.02) + 
        ggplot2::geom_hline(yintercept=0) + 
        ggplot2::geom_point(aes(tau, qte)) +
        ggplot2::ggtitle(main) + 
        ggplot2::scale_y_continuous(ylab, limits=ylim, breaks=ybreaks) +
        ggplot2::scale_x_continuous("tau", limits=c(0,1), breaks=xbreaks) + 
        ggplot2::theme_classic() +
        ggplot2::theme(panel.border = element_rect(colour = 'black', size=1,
                                          fill=NA,
                                          linetype='solid'),
                       plot.title = element_text(hjust=0.5))
    qp
}

##functions to get median (or specified quantile)
##requires that qte object has that value of tau
getMedian <- function(qteobj, tau=.5) {
    which.qte <- which(qteobj$probs == tau)
    return(qteobj$qte[which.qte])
}

##functions to get median (or specified quantile) standard error
##requires that qte object has that value of tau
getMedianSE <- function(qteobj, tau=.5) {
    which.qte <- which(qteobj$probs == tau)
    return(qteobj$qte.se[which.qte])
}

##functions to get 80-20 difference or some other difference in quantiles
diffquantiles <- function(qteobj, hightau, lowtau) {
    which.highqte <- which(qteobj$probs == hightau)
    which.lowqte <- which(qteobj$probs == lowtau)
    return(qteobj$qte[which.highqte] - qteobj$qte[which.lowqte])
}


##bootstrap the difference betwen quantiles
##must be called with a qteobj with retEachIter set to true
bootse.diffquantiles <- function(qteobj, hightau, lowtau) {
    bootvals <- lapply(qteobj$eachIterList, diffquantiles, hightau, lowtau)
    se <- stats::sd(unlist(bootvals))
}

#'@title diffQ
#'
#' @description ## takes a single set of quantiles
#' (not qtes as in diffquantiles)
#'  and returns the difference between particular ones
#'
#' @param qvec vector of quantiles
#' @param tauvec vector of tau (should be same length as qvec)
#' @param hightau upper quantile
#' @param lowtau lower quantile
#'
#' @return scalar difference between quantiles
#' @export
diffQ <- function(qvec, tauvec, hightau, lowtau) {
    which.highq <- which(tauvec == hightau)
    which.lowq <- which(tauvec == lowtau)
    return(qvec[which.highq] - qvec[which.lowq])
}



##make tables using R's texreg package
#'@title diffQ
#'
#' @description ## takes a single set of quantiles
#' (not qtes as in diffquantiles)
#'  and returns the difference between particular ones
#'
#' @param qteobj A QTE object
#' @param tau Optional vector of taus to texreg results for
#' @param reportAte Whether or not texreg the ATE (or ATT) as well
#'
qteToTexreg <- function(qteobj, tau=NULL, reportAte=T) {
    if (is.null(tau)) {
        tau <- qteobj$probs
        qte <- qteobj$qte
        qte.se <- qteobj$qte.se
        ate <- qteobj$ate
        ate.se <- qteobj$ate.se        
    } else if(!(all(tau %in% qteobj$probs))) {
        stop("Error not all tau in qte object")
    } else {
        tauloc <- vapply(tau, function(x) { which(x==qteobj$probs) }, 1.0)
        qte <- qteobj$qte[tauloc]
        qte.se <- qteobj$qte.se[tauloc]
        ate <- qteobj$ate
        ate.se <- qteobj$ate.se
    }
    if (reportAte) {
        texreg::createTexreg(c(paste(tau), "ate"),
                     c(qte, ate),
                     c(qte.se, ate.se),
                     2*pnorm(-c(abs(qte/qte.se),
                                abs(ate/ate.se))))
    } else {
        texreg::createTexreg(paste(tau),
                     qte,
                     qte.se,
                     2*pnorm(-c(abs(qte/qte.se))))
    }
}

