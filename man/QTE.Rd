% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R, R/qte.R
\docType{package}
\name{QTE}
\alias{QTE}
\title{qte: A package for computating quantile treatment effects}
\usage{
QTE(
  qte,
  ate = NULL,
  qte.se = NULL,
  qte.lower = NULL,
  qte.upper = NULL,
  ate.se = NULL,
  ate.lower = NULL,
  ate.upper = NULL,
  c = NULL,
  alp = 0.05,
  pscore.reg = NULL,
  probs,
  type = "On the Treated",
  F.treated.t = NULL,
  F.untreated.t = NULL,
  F.treated.t.cf = NULL,
  F.treated.tmin1 = NULL,
  F.treated.tmin2 = NULL,
  F.treated.change.tmin1 = NULL,
  F.untreated.change.t = NULL,
  F.untreated.change.tmin1 = NULL,
  F.untreated.tmin1 = NULL,
  F.untreated.tmin2 = NULL,
  condQ.treated.t = NULL,
  condQ.treated.t.cf = NULL,
  eachIterList = NULL,
  inffunct = NULL,
  inffuncu = NULL
)
}
\arguments{
\item{qte}{The Quantile Treatment Effect at each value of probs}

\item{ate}{The Average Treatment Effect (or Average Treatment Effect on
the Treated)}

\item{qte.se}{A vector of standard errors for each qte}

\item{qte.lower}{A vector of lower confidence intervals for each qte (it is
based on the bootstrap confidence interval -- not the se -- so it may not
be symmyetric about the qte}

\item{qte.upper}{A vector of upper confidence intervals for each qte (it is
based on the bootstrap confidence interval -- not the se -- so it may not
be symmetric about the qte}

\item{ate.se}{The standard error for the ATE}

\item{ate.lower}{Lower confidence interval for the ATE (it is based on the
bootstrap confidence intervall -- not the se -- so it may not be symmetric
about the ATE}

\item{ate.upper}{Upper confidence interval for the ATE (it is based on the
bootstrap confidence interval -- not the se -- so it may not be symmetric
about the ATE}

\item{c}{The critical value from a KS-type statistic used for creating
uniform confidence bands}

\item{alp}{The significance level}

\item{pscore.reg}{The results of propensity score regression, if specified}

\item{probs}{The values for which the qte is computed}

\item{type}{Takes the values "On the Treated" or "Population" to indicate
whether the estimated QTE is for the treated group or for the entire
population}

\item{F.treated.t}{Distribution of treated outcomes for the treated group at
period t}

\item{F.untreated.t}{Distribution of untreated potential outcomes for the
untreated group at period t}

\item{F.treated.t.cf}{Counterfactual distribution of untreated potential
outcomes for the treated group at period t}

\item{F.treated.tmin1}{Distribution of treated outcomes for the
treated group at period tmin1}

\item{F.treated.tmin2}{Distribution of treated outcomes for the
treated group at period tmin2}

\item{F.treated.change.tmin1}{Distribution of the change in outcomes for
the treated group between periods tmin1 and tmin2}

\item{F.untreated.change.t}{Distribution of the change in outcomes for the
untreated group between periods t and tmin1}

\item{F.untreated.change.tmin1}{Distribution of the change in outcomes for
the untreated group between periods tmin1 and tmin2}

\item{F.untreated.tmin1}{Distribution of outcomes for the untreated group
in period tmin1}

\item{F.untreated.tmin2}{Distribution of outcomes for the untreated group
in period tmin2}

\item{condQ.treated.t}{Conditional quantiles for the treated group in
period t}

\item{condQ.treated.t.cf}{Counterfactual conditional quantiles for the treated
group in period t}

\item{eachIterList}{An optional list of the outcome of each bootstrap
iteration}

\item{inffunct}{The influence function for the treated group;
used for inference when there are multiple
periods and in the case with panel data.  It is needed for computing covariance
terms in the variance-covariance matrix.}

\item{inffuncu}{The influence function for the untreated group}
}
\description{
Main class of objects.  A \code{QTE} object is returned by
all of the methods that compute the QTE or QTET.
}
