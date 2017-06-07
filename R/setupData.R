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
#' @return qteData object holding data to be used in QTE functions
#' @export
setupData <- function(qteParams) {

    formla <- qteParams$formla
    xformla <- qteParams$xformla
    data <- qteParams$data
    w <- qteParams$w
    probs <- qteParams$probs
    method <- qteParams$method
    
    form <- stats::as.formula(formla)
    dta <- stats::model.frame(stats::terms(form,data=data),data=data) #or model.matrix
    colnames(dta) <- c("y","treatment")
    yname <- "y"
    idname <- qteParams$idname
    tname <- qteParams$tname
    treat <- "treatment"
    wname <- "w"
    t <- qteParams$t
    tmin1 <- qteParams$tmin1
    tmin2 <- qteParams$tmin2
    panel <- qteParams$panel
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
        ##data <- cbind(data[,c(yname,treat,wname,tname)],
        ##              model.matrix(terms(as.formula(xformla)), data=data))
    } else {
        x <- NULL
    }

    ##setup the data
    treated <- data[data[,treat]==1,]
    untreated <- data[data[,treat]==0,]


    env <- parent.frame()

    ##with multiple periods
    if (!is.null(tname)) {

        if (panel) {
            data <- makeBalancedPanel(data, idname, tname)
        }

        treated.t <- data[ data[,tname] == t & data[,treat]==1, ]
        treated.tmin1 <- data[ data[,tname] == tmin1 & data[,treat]==1, ]
        untreated.t <- data[ data[,tname] == t & data[,treat]==0, ]
        untreated.tmin1 <- data[ data[,tname] == tmin1 & data[,treat]==0, ]
      

        ##Sort data and rely on having a balanced panel to get the change
        ## distributions right
        if (panel) {
            treated.t <- treated.t[order(treated.t[,idname]),]
            treated.tmin1 <- treated.tmin1[order(treated.tmin1[,idname]),]
            untreated.t <- untreated.t[order(untreated.t[,idname]),]
            untreated.tmin1 <- untreated.tmin1[order(untreated.tmin1[,idname]),]
        }
        
        ##3) Get the distributions that we need below
        
        ##a) Get distribution of y0.tmin2 | Dt=1
        F.treated.t <- stats::ecdf(treated.t[,yname])
        F.treated.tmin1 <- stats::ecdf(treated.tmin1[,yname])
        F.untreated.t <- stats::ecdf(untreated.t[,yname])
        F.untreated.tmin1 <- stats::ecdf(untreated.tmin1[,yname])

        
        assign("treated.t", treated.t, envir=env)
        assign("treated.tmin1", treated.tmin1, envir=env)
        assign("untreated.t", untreated.t, envir=env)
        assign("untreated.tmin1", untreated.tmin1, envir=env)
        assign("F.treated.t", F.treated.t, envir=env)
        assign("F.treated.tmin1", F.treated.tmin1, envir=env)
        assign("F.untreated.t", F.untreated.t, envir=env)
        assign("F.untreated.tmin1", F.untreated.tmin1, envir=env)

        if (!is.null(tmin2)) {
            treated.tmin2 <- data[ data[,tname] == tmin2 & data[,treat]==1, ]
            untreated.tmin2 <- data[ data[,tname] == tmin2 & data[,treat]==0, ]
            if (panel) {
                treated.tmin2 <- treated.tmin2[order(treated.tmin2[,idname]),]
                untreated.tmin2 <- untreated.tmin2[order(untreated.tmin2[,idname]),]
            }
            F.treated.tmin2 <- ecdf(treated.tmin2[,yname])
            F.untreated.tmin2 <- ecdf(untreated.tmin2[,yname])
            
            assign("treated.tmin2", treated.tmin2, envir=env)
            assign("untreated.tmin2", untreated.tmin2, envir=env)
            assign("F.treated.tmin2", F.treated.tmin2, envir=env)
            assign("F.untreated.tmin2", F.untreated.tmin2, envir=env)
        }

    }


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
    assign("panel", panel, envir=env)
    
}
