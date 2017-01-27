#' @title DR
#'
#' @description A distribution regression object
#'
#' @param yvals A vector of values that y can take
#' @param drlist A list where for each value of y, a distribution regression
#' 
#' @export
DR <- function(yvals, drlist) {
    out <- list(yvals=yvals, drlist=drlist)
    class(out) <- "DR"
    return(out)
}
