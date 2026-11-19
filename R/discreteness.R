# =============================================================================
# Title: Discrete-outcome check
# Description: Quantile treatment effects are built for continuous outcomes.
#   outcome_discreteness() summarizes how discrete an outcome is (distinct
#   values, largest weighted mass point); check_discrete_outcome() warns once
#   per estimator call when the outcome looks discrete (issue #17). The
#   estimators call it on the whole outcome column, across all groups and
#   periods, before estimation; the *_gt() functions never warn, so the
#   bootstrap can't flood the console.
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-11-19
# =============================================================================

# --- outcome_discreteness ----------------------------------------------------

# Summarize how discrete an outcome is.
#
# Args:
#   y  outcome vector (missing values are dropped, with their weights)
#   w  optional sampling weights, same length as y
#
# Returns a list: n (observations used), n_unique (distinct values),
# max_mass (largest weighted share of observations at a single value) and
# mass_at (that value).
outcome_discreteness <- function(y, w = NULL) {
  if (is.null(w)) w <- rep(1, length(y))
  keep <- !is.na(y) & !is.na(w)
  y <- y[keep]
  w <- w[keep]
  if (length(y) == 0L || sum(w) <= 0) {
    return(list(n = length(y), n_unique = length(unique(y)),
                max_mass = NA_real_, mass_at = NA_real_))
  }
  mass <- tapply(w, y, sum) / sum(w)
  top <- which.max(mass)
  list(n = length(y), n_unique = length(mass),
       max_mass = as.numeric(mass[top]),
       mass_at = as.numeric(names(mass)[top]))
}

# --- check_discrete_outcome --------------------------------------------------

# Warn (once) when the outcome looks discrete, and return the diagnostics.
#
# The outcome counts as discrete if its largest mass point holds more than
# getOption("qte.mass_point_tol", 0.10) of the (weighted) observations, or if
# it takes at most getOption("qte.discrete_n_unique", 20) distinct values.
# The 10% default flags zero-inflated outcomes and binary or Likert outcomes,
# but not log counts with a mean around 20 (largest mass about 9%) or
# continuous data rounded to one or two decimals (at most about 4%). Setting
# both options to 1 and 0 turns the check off.
#
# The warning has class "qte_discrete_outcome", so it can be caught or
# muffled by class.
#
# Args:
#   data         the estimator's data
#   yname        outcome column
#   weightsname  weights column or NULL
#   warn         FALSE to only compute the diagnostics (e.g. ATT runs of
#                estimators whose ATT doesn't depend on the outcome's
#                distribution)
#   advice       closing sentence of the warning
#
# Returns the outcome_discreteness() list plus `discrete` (TRUE/FALSE),
# invisibly.
check_discrete_outcome <- function(data, yname, weightsname = NULL,
                                   warn = TRUE,
                                   advice = paste0(
                                     "Make sure this is what you intend; ",
                                     "cdf_did() is identified for discrete ",
                                     "outcomes."
                                   )) {
  y <- data[[yname]]
  w <- if (is.null(weightsname)) NULL else data[[weightsname]]
  out <- outcome_discreteness(y, w)
  tol <- getOption("qte.mass_point_tol", 0.10)
  max_unique <- getOption("qte.discrete_n_unique", 20)
  out$discrete <- isTRUE(out$max_mass > tol) || out$n_unique <= max_unique

  if (warn && out$discrete) {
    found <- if (out$n_unique <= max_unique) {
      sprintf("it takes only %d distinct values", out$n_unique)
    } else {
      sprintf("%.0f%% of observations equal %s",
              100 * out$max_mass, format(out$mass_at, digits = 6))
    }
    msg <- paste0(
      "The outcome `", yname, "` looks discrete: ", found, ". ",
      "Quantile treatment effects usually require a continuous outcome. ",
      "With mass points, quantiles jump and effects at quantiles inside a ",
      "mass point are 0 by construction; rank- and quantile-based ",
      "identification (e.g. change-in-changes) then only bounds the ",
      "effects; and bootstrap inference for quantiles can be unreliable. ",
      advice
    )
    warning(warningCondition(msg, class = "qte_discrete_outcome"))
  }
  invisible(out)
}
