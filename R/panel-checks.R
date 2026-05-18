# =============================================================================
# Title: Panel data validation and formatting utilities
# Description: panel.checks validates that treated group composition is stable
#   across time periods and warns if covariates vary over time.
#   panelize.data subsets and balances data for use with panel estimators.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

#' @title panel.checks
#'
#' @description Validates that the panel data structure is compatible with
#'  the panel estimators in the qte package. Specifically: checks that the
#'  treated group composition is constant across all time periods, warns if
#'  covariates appear to vary over time, and stops with an informative message
#'  if the treated group is empty or changes.
#'
#' @inheritParams compute.panel.qtet
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects
#'  (warnings and errors).
#'
#' @keywords internal
panel.checks <- function(qp) {
  setupData(qp)

  ## 1) give warning if covariates appear to vary over time
  if (is.null(xformla)) xformla1 <- ~1 else xformla1 <- xformla
  dft <- model.frame(xformla1, data = rbind.data.frame(treated.t, untreated.t))
  dftmin1 <- model.frame(xformla1, data = rbind.data.frame(treated.tmin1, untreated.tmin1))

  if (is.null(tmin2)) {
    ## this is just hack way to use the same checks for 2 period and 3 period case
    dftmin2 <- dftmin1
    treated.tmin2 <- treated.tmin1
  } else {
    dftmin2 <- model.frame(xformla1, data = rbind.data.frame(treated.tmin2, untreated.tmin2))
  }

  if (!(identical(dftmin2, dftmin1) && identical(dftmin2, dft))) {
    warning("covariates appear to vary over time...\n  only conditioning on first period covariates...\n  this is recommended practice, but worth noting...")
  }

  ## 2) some basic error handling about treated group being constant over time
  if ((nrow(treated.tmin2) == 0 || nrow(treated.tmin1) == 0)) {
    stop("Treatment status should be equal to 1 for all individuals in the treated group -- that is, individuals that ever become treated")
  }

  if (!(all.equal(treated.t[, idname], treated.tmin1[, idname]) &&
    all.equal(treated.t[, idname], treated.tmin2[, idname]))) {
    stop("The composition of the treated group is changing over time...\n  The treated group should consistent of observations that are first treated in the last period...\n  The treatment variable should be set equal to 1 in all time periods for this group... \n Individuals that are treated before the last time period should be removed from the dataset as treatment effects are not identified for this group...")
  }
}

#' @title panelize.data
#'
#' @description get data in correct format for using panel methods in qte package
#' @inheritParams panel.qtet
#'
#' @return data.frame
#'
#' @export
panelize.data <- function(data, idname, tname, t, tmin1, tmin2 = NULL) {
  ## drop obs that are not in period t, tmin1, tmin2
  pren <- nrow(data)
  data <- data[data[, tname] %in% c(t, tmin1, tmin2), ]
  if (!(nrow(data) == pren)) {
    warning(paste0("dropping ", pren - nrow(data), " observations that are not in period: ", t, ", ", tmin1, ", ", tmin2, "..."))
  }

  ## make a balanced panel
  pren <- nrow(data)
  data <- setDF(makeBalancedPanel(data, idname, tname))
  if (!(nrow(data) == pren)) {
    warning(paste0("enforcing balanced panel condition...\n  this drops ", (pren - nrow(data)) / 3, " observations..."))
  }

  return(data)
}
