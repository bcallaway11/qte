# =============================================================================
# Title: Mean Difference-in-Differences (MDiD)
# Description: Implements the MDiD estimator: the counterfactual quantile is
#   constructed by shifting treated pre-period quantiles by the mean DiD
#   (rather than the full distributional shift used in QDiD/CiC).
#   Note: compute.MDiD does its own data setup rather than using setupData(),
#   which is inconsistent with other compute.* functions in the package.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

# --- Internal compute function -----------------------------------------------

### Mean Difference-in-Differences
## Note that you need to pass in data where treated status is noted in
## every period.  Data is form of (year-individual-outcome-x-evertreated)
#' @title compute.MDiD
#'
#' @description Internal function for computing the MDiD point estimate.
#'  Unlike other \code{compute.*} functions in the package, this function
#'  performs its own data setup rather than calling \code{setupData}.
#'
#' @inheritParams panel.qtet
#'
#' @return QTE object
#'
#' @keywords internal
compute.MDiD <- function(formla, xformla = NULL, t, tmin1, tname, data,
                         panel = FALSE,
                         idname = NULL, probs = seq(0.05, 0.95, 0.05)) {
  form <- as.formula(formla)
  dta <- model.frame(terms(form, data = data), data = data) # or model.matrix
  colnames(dta) <- c("y", "treatment")
  yname <- "y"
  treat <- "treatment"
  data <- cbind.data.frame(dta, data)

  ## Setup x variables if using formula
  x <- NULL
  if (!(is.null(xformla))) {
    ## in this case, we need to drop the intercept
    x <- colnames(model.matrix(terms(as.formula(xformla)), data = data))[-1]
    data <- cbind(
      data[, c(yname, treat, idname, tname)],
      model.matrix(terms(as.formula(xformla)), data = data)
    )[, c(1:4, 6:(5 + length(x)))]
  }

  ## drop the always treated.  Note that this also relies
  ## on no "switchback" or no return to untreated status
  ## after joining treatment.
  ## first line gets the correct two years of data
  data <- subset(data, (data[, tname] == tmin1 | data[, tname] == t))

  if (panel) {
    data <- makeBalancedPanel(data, idname, tname)
  }

  ## just to make sure the factors are working ok
  data <- droplevels(data)

  ## adjust for covariates
  ## after adjustment then everything should proceed as before
  if (!(is.null(x))) {
    cov.data <- data
    cov.data$group <- as.factor(paste(cov.data[, treat],
      cov.data[, tname],
      sep = "-"
    ))
    group <- "group"
    xmat <- cov.data[, x]
    first.stage <- lm(cov.data[, yname] ~ -1 + cov.data[, group] +
      as.matrix(xmat))
    ## get residuals not including group dummies
    bet <- coef(first.stage)[5:length(coef(first.stage))]
    yfit <- cov.data[, yname] - as.matrix(xmat) %*% bet
    data[, yname] <- yfit
  }


  ## Setup each of the datasets used below
  ## a) get all the treated (in the last period) observations
  treated.t <- data[data[, tname] == t & data[, treat] == 1, ]

  ## b) set ever.treated to 1 if observation is treated in last period
  ## Try not to use this b/c won't work in the case with repeated cross sections
  ## data$ever.treated = data$treatment
  ## data$ever.treated = 1*(data[,idname] %in% treated.t[,idname])
  ## ever.treated = "ever.treated"

  ## Setup each of the datasets used below
  ## treated.t = subset(data, data[,treat]==1 & data[,tname]==t)
  ## just get the lagged value of y; otherwise keep the same
  ## dataset.  Note: this will not work if there are x covariates;
  ## well, could follow similar procedure, but as is, would
  ## require some modification.
  treated.tmin1 <- data[data[, tname] == tmin1 & data[, treat] == 1, ]
  ## this is right
  untreated.t <- data[data[, tname] == t & data[, treat] == 0, ]
  ## get lagged of untreated y
  untreated.tmin1 <- data[data[, tname] == tmin1 & data[, treat] == 0, ]

  ## 5) Compute Quantiles
  ## a) Quantiles of observed distribution
  q1 <- quantile(treated.t[, yname], probs = probs)
  q0 <- quantile(treated.tmin1[, yname], probs = probs) + mean(untreated.t[, yname]) - mean(untreated.tmin1[, yname])


  ## 7) Estimate ATT using MDID
  att <- mean(treated.t[, yname]) - (mean(treated.tmin1[, yname]) + mean(untreated.t[, yname]) - mean(untreated.tmin1[, yname]))


  out <- QTE(ate = att, qte = (q1 - q0), probs = probs)

  return(out)
}

# --- User-facing wrapper -----------------------------------------------------

#' @title Mean Difference in Differences
#'
#' @description \code{MDiD} computes the Quantile Treatment Effect on the
#' Treated (QTET) using a Mean Difference-in-Differences approach. The
#' counterfactual quantile is formed by shifting the pre-treatment treated
#' distribution by the mean DiD (i.e., the parallel-trends DiD in means),
#' rather than the full distributional shift used in \code{QDiD} or \code{CiC}.
#'
#' The method can accommodate covariates by first residualizing outcomes on a
#' linear model with group-time fixed effects and covariates, then applying
#' the MDiD estimator to the quasi-residuals.
#'
#' @inheritParams panel.qtet
#' @inheritParams CiC
#'
#' @examples
#' ## load the data
#' data(lalonde)
#'
#' ## Run the Mean Difference in Differences method conditioning on
#' ## age, education, black, hispanic, married, and nodegree
#' md1 <- MDiD(re ~ treat,
#'   t = 1978, tmin1 = 1975, tname = "year",
#'   xformla = ~ age + I(age^2) + education + black + hispanic + married + nodegree,
#'   data = lalonde.psid.panel, idname = "id", se = FALSE,
#'   probs = seq(0.05, 0.95, 0.05)
#' )
#' summary(md1)
#'
#' @references
#' Athey, Susan and Guido Imbens.  ``Identification and Inference in Nonlinear
#'  Difference-in-Differences Models.'' Econometrica 74.2, pp. 431-497,
#'  2006.
#'
#' Thuysbaert, Bram.  ``Distributional Comparisons in Difference in Differences
#'  Models.'' Working Paper, 2007.
#'
#' @param biters Number of bootstrap iterations; alias for \code{iters}
#'   matching the \code{did}/\code{ptetools} naming convention. If both are
#'   supplied, \code{biters} takes precedence. Note: \code{MDiD} does not
#'   support parallel computation.
#' @return A \code{QTE} object
#'
#' @export
MDiD <- function(formla, xformla = NULL, t, tmin1, tname, data,
                 panel = FALSE, se = TRUE,
                 idname = NULL,
                 alp = 0.05, probs = seq(0.05, 0.95, 0.05), iters = 100,
                 retEachIter = FALSE,
                 biters = NULL) {
  # biters is an alias matching the did/ptetools API (MDiD has no parallel support)
  if (!is.null(biters)) iters <- biters
  form <- as.formula(formla)
  dta <- model.frame(terms(form, data = data), data = data) # or model.matrix
  colnames(dta) <- c("y", "treatment")
  yname <- "y"
  treat <- "treatment"
  data <- cbind.data.frame(dta, data)

  ## drop the always treated.  Note that this also relies
  ## on no "switchback" or no return to untreated status
  ## after joining treatment.
  ## first line gets the correct two years of data
  data <- subset(data, (data[, tname] == tmin1 | data[, tname] == t))

  if (panel) {
    if (is.null(idname)) {
      stop("Must provide idname when using panel option")
    }
    data <- makeBalancedPanel(data, idname, tname)
  }


  ## just to make sure the factors are working ok
  data <- droplevels(data)

  ## Setup each of the datasets used below
  ## a) get all the treated (in the last period) observations
  treated.t <- data[data[, tname] == t & data[, treat] == 1, ]
  treated.tmin1 <- data[data[, tname] == tmin1 & data[, treat] == 1, ]
  untreated.t <- data[data[, tname] == t & data[, treat] == 0, ]
  untreated.tmin1 <- data[data[, tname] == tmin1 & data[, treat] == 0, ]

  ## first calculate the actual estimate
  mdid <- compute.MDiD(
    formla, xformla, t, tmin1, tname, data,
    panel, idname, probs
  )

  if (se) {
    ## now calculate the bootstrap confidence interval
    eachIter <- list()
    ## Need to build dataset by sampling individuals, and then
    ## taking all of their time periods
    ## when it's a panel make draws by individual
    if (panel) {
      ## all.ids = unique(data[,idname])
      ## here we rely on having a balanced panel to get the right obs.
      treated.t <- treated.t[order(treated.t[, idname]), ]
      treated.tmin1 <- treated.tmin1[order(treated.tmin1[, idname]), ]
      untreated.t <- untreated.t[order(untreated.t[, idname]), ]
      untreated.tmin1 <- untreated.tmin1[order(untreated.tmin1[, idname]), ]
      nt <- nrow(treated.t)
      nu <- nrow(untreated.t)
      ## out.bootdatalist <<- list()
      for (i in 1:iters) {
        ## reset boot.data
        ## boot.data = data[0,]
        randy.t <- sample(seq_len(nt), nt, replace = TRUE)
        randy.u <- sample(seq_len(nu), nu, replace = TRUE)
        boot.data.treated.t <- treated.t[randy.t, ]
        boot.data.treated.tmin1 <- treated.tmin1[randy.t, ]
        boot.data.untreated.t <- untreated.t[randy.u, ]
        boot.data.untreated.tmin1 <- untreated.tmin1[randy.u, ]
        boot.data <- rbind(
          boot.data.treated.t, boot.data.untreated.t,
          boot.data.treated.tmin1,
          boot.data.untreated.tmin1
        )
        thisIter <- compute.MDiD(formla, xformla, t, tmin1, tname,
          boot.data,
          panel = FALSE, idname, probs
        )
        ## already have a balanced panel so can increase speed by calling
        ## with panel option set to F.
        eachIter[[i]] <- QTE(
          ate = thisIter$ate, qte = thisIter$qte,
          probs = probs
        )
      }
    } else { # make draws within each sample
      treated.t <- data[data[, tname] == t & data[, treat] == 1, ]
      treated.tmin1 <- data[data[, tname] == tmin1 & data[, treat] == 1, ]
      untreated.t <- data[data[, tname] == t & data[, treat] == 0, ]

      untreated.tmin1 <- data[data[, tname] == tmin1 & data[, treat] == 0, ]

      for (i in 1:iters) {
        n <- nrow(treated.t)
        ran <- sample(seq_len(n), n, replace = TRUE)
        boot.treated.t <- treated.t[ran, ]

        n <- nrow(treated.tmin1)
        ran <- sample(seq_len(n), n, replace = TRUE)
        boot.treated.tmin1 <- treated.tmin1[ran, ]

        n <- nrow(untreated.t)
        ran <- sample(seq_len(n), n, replace = TRUE)
        boot.untreated.t <- untreated.t[ran, ]

        n <- nrow(untreated.tmin1)
        ran <- sample(seq_len(n), n, replace = TRUE)
        boot.untreated.tmin1 <- untreated.tmin1[ran, ]

        boot.data <- rbind(
          boot.treated.t, boot.untreated.t,
          boot.treated.tmin1, boot.untreated.tmin1
        )
        thisIter <- compute.MDiD(
          formla, xformla, t, tmin1, tname,
          boot.data,
          panel, idname, probs
        )
        eachIter[[i]] <- QTE(
          ate = thisIter$ate, qte = thisIter$qte,
          probs = probs
        )
      }
    }

    SEobj <- computeSE(eachIter, mdid, alp = alp)

    if (!retEachIter) {
      eachIter <- NULL
    }

    out <- QTE(
      qte = mdid$qte, qte.upper = SEobj$qte.upper,
      qte.lower = SEobj$qte.lower, ate = mdid$ate,
      ate.upper = SEobj$ate.upper, ate.lower = SEobj$ate.lower,
      qte.se = SEobj$qte.se, ate.se = SEobj$ate.se,
      alp = alp,
      eachIterList = eachIter,
      probs = probs
    )
    return(out)
  } else {
    return(mdid)
  }
}
