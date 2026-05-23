# =============================================================================
# Title: Dataset documentation
# Description: Roxygen2 documentation for the lalonde datasets bundled with
#   the qte package. Actual data objects live in data/lalonde.RData.
# Author: Brant Callaway
# Last update: 2026-05-22
# Date created: 2026-05-22
# =============================================================================

#' @title Lalonde (1986)'s NSW Dataset
#'
#' @description \code{lalonde} contains data from the National Supported Work
#'  Demonstration.  This program randomly assigned applicants to the job
#'  training program (or out of the job training program).  The dataset is
#'  discussed in Lalonde (1986).  The experimental part of the dataset is
#'  combined with an observational dataset from the Panel Study of Income
#'  Dynamics (PSID).  Lalonde (1986) and many subsequent papers (e.g.
#'  Heckman and Hotz (1989), Dehejia and Wahba (1999), Smith and Todd (2005),
#'  and Firpo (2007) have used this combination to study the effectiveness
#'  of various `observational' methods (e.g. regression, Heckman selection,
#'  Difference in Differences, and propensity score matching) of estimating
#'  the Average Treatment Effect (ATE) of participating in the job training
#'  program.  The idea is that the results from the observational method
#'  can be compared to results that can be easily obtained from the
#'  experimental portion of the dataset.
#'
#'  To be clear, the observational data combines the observations that are
#'  treated from the experimental portion of the data with untreated observations
#'  from the PSID.
#'
#' @format Four data.frames: (i) lalonde.exp contains a cross sectional version
#'  of the experimental data, (ii) lalonde.psid contains a cross sectional
#'  version of the observational data, (iii) lalonde.exp.panel contains a
#'  panel version of the experimental data, and (iv) lalonde.psid.panel contains
#'  a panel version of the observational data.  Note: the cross sectional
#'  and panel versions of each dataset are identical up to their shape; in
#'  demonstrating each of the methods, it is sometimes convenient to have
#'  one form of the data or the other.
#' @docType data
#' @name lalonde
#' @usage data(lalonde)
#' @references LaLonde, Robert.  ``Evaluating the Econometric Evaluations of
#'  Training Programs with Experimental Data.'' The American Economics Review,
#'  pp. 604-620, 1986.
#'  @source The dataset comes from Lalonde (1986) and has been studied in much
#'  subsequent work.  The \code{qte} package uses a version from the
#'  \code{causalsens} package
#'  (\url{https://CRAN.R-project.org/package=causalsens})
#' @keywords datasets
NULL

#' @title Lalonde's Experimental Dataset
#'
#' @description The cross sectional verion of the experimental part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.exp
#' @keywords datasets
NULL

#' @title Lalonde's Panel Experimental Dataset
#'
#' @description The panel verion of the experimental part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.exp.panel
#' @keywords datasets
NULL

#' @title Lalonde's Observational Dataset
#'
#' @description The cross sectional verion of the observational part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.psid
#' @keywords datasets
NULL

#' @title Lalonde's Panel Observational Dataset
#'
#' @description The panel verion of the observational part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.psid.panel
#' @keywords datasets
NULL
