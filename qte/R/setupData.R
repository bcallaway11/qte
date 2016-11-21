##
####Setup data for using in compute methods in QTE package####
#' @title setupData
#'
#' @description \code{setupData} sets up the data to use in each
#'  compute method in the QTE package
#' 
#' @inheritParams panel.qtet
#' @param qteParams object holding the function parameters
#' 
#' @keywords internal
#' 
#' @return qteData object holding data to be used in QTE functions
setupData <- function(qteParams) {

    formla <- qteParams$formla
    xformla <- qteParams$xformla
    data <- qteParams$data
    w <- qteParams$w
    probs <- qteParams$probs
    method <- qteParams$method
    
    form <- as.formula(formla)
    dta <- model.frame(terms(form,data=data),data=data) #or model.matrix
    colnames(dta) <- c("y","treatment")
    yname <- "y"
    tname <- qteParams$tname
    treat <- "treatment"
    wname <- "w"
    t <- qteParams$t
    tmin1 <- qteParams$tmin1
    tmin2 <- qteParams$tmin2
    data <- cbind.data.frame(dta,data)

    if (is.null(w)) {
        data$w <- 1
    } else {
        data$w <- w
    }

    if (sum(data$w) != nrow(data)) {
        warning("sum of weights not equal to number of rows of data; something may be wrong...Proceeding")
    }

    ##set up the x variables
    if (!(is.null(xformla))) {
        x <- colnames(model.matrix(terms(as.formula(xformla)), data=data))
        data <- cbind(data[,c(yname,treat,wname)],
                      model.matrix(terms(as.formula(xformla)), data=data))
        ##if (!is.null(weights)) {
        ##    data <- cbind(data, weights=ww)
        ##} else{
        ##    data$weights <- 1
        ##}
    } else {
        x <- NULL
    }

    ##setup the data
    treated <- data[data[,treat]==1,]
    untreated <- data[data[,treat]==0,]

    env <- parent.frame()

    assign("formla", formla, envir=env)
    assign("xformla", xformla, envir=env)
    assign("x", x, envir=env)
    assign("yname", "y", envir=env)
    assign("treat", "treatment", envir=env)
    assign("wname", "w", envir=env)    
    assign("data", data, envir=env)
    assign("treated", treated, envir=env)
    assign("untreated", untreated, envir=env)
    assign("probs", probs, envir=env)
    assign("method", method, envir=env)
    assign("tname", tname, envir=env)
    assign("t", t, envir=env)
    assign("tmin1", tmin1, envir=env)
    assign("tmin2", tmin2, envir=env)
    
}
