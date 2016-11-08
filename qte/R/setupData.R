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
    weights <- qteParams$weights
    probs <- qteParams$probs
    method <- qteParams$method
    
    form <- as.formula(formla)
    dta <- model.frame(terms(form,data=data),data=data) #or model.matrix
    colnames(dta) <- c("y","treatment")
    yname <- "y"
    treat <- "treatment"
    wname <- "weights"
    data <- cbind.data.frame(dta,data)
    
    if (!is.null(weights)) {
        ww <- data$weights / sum(data$weights)
    }

    ##set up the x variables
    if (!(is.null(xformla))) {
        x <- colnames(model.matrix(terms(as.formula(xformla)), data=data))
        data <- cbind(data[,c(yname,treat)],
                      model.matrix(terms(as.formula(xformla)), data=data))
        if (!is.null(weights)) {
            data <- cbind(data, weights=ww)
        } else{
            data$weights <- 1/nrow(data)
        }
    } else {
        x <- NULL
        data$weights <-  1/nrow(data)
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
    assign("wname", "weights", envir=env)    
    assign("data", data, envir=env)
    assign("treated", treated, envir=env)
    assign("untreated", untreated, envir=env)
    assign("probs", probs, envir=env)
    assign("method", method, envir=env)
    
}
