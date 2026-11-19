# =============================================================================
# Title: Fast smoke tests — panel_qtt() distributional panel QTT estimator
# Description: Regression checks for panel_qtt() covering panel ATT and QTT
#   modes. Uses a synthetic 4-period staggered panel (1 group treated at t=3)
#   so that pre2 = period 1 and pre1 = period 2 are always available.
#   biters kept low (20) so the suite runs quickly. Reference values computed
#   from the first verified run. Also checks pre_copula = "long" with
#   unequally spaced calendar years (lalonde.psid.panel: 1974, 1975, 1978).
# Author: Brant Callaway
# Last update: 2026-11-19
# Date created: 2026-05-20
# =============================================================================

# Synthetic 4-period panel: never-treated (group 0) + one group treated at t=3.
# Periods 1..4 give pre2=1, pre1=2, post={3,4} for the treated group.
make_panel_qtt_data <- function(seed = 42) {
  set.seed(seed)
  n <- 30L; periods <- 1:4; groups <- c(0, 3)
  do.call(rbind, lapply(groups, function(g) {
    do.call(rbind, lapply(seq_len(n), function(i) {
      id <- (match(g, groups) - 1L) * n + i
      data.frame(id = id, period = periods, group = g,
                 y = rnorm(4L) + ifelse(g > 0 & periods >= g, 1, 0))
    }))
  }))
}

# Panel ATT (pre_copula = "short"): 2 ATT(g,t) cells (1 group x 2 post periods),
# no NAs, overall ATT near 1 (true treatment effect in the synthetic data).
# Using "short" here so both post-treatment cells are valid with a 4-period dataset.
test_that("panel_qtt panel ATT returns valid pte_emp_boot with stable estimate", {
  dat <- make_panel_qtt_data()
  set.seed(42)
  res <- panel_qtt(yname = "y", gname = "group", tname = "period", idname = "id",
                   data = dat, biters = 20, pre_copula = "short")
  expect_s3_class(res, "pte_emp_boot")
  expect_true(is.numeric(res$overall_results$att))
  expect_false(anyNA(res$attgt_results$att))
  expect_equal(nrow(res$attgt_results), 2L)
  expect_equal(res$overall_results$att, 1.057125, tolerance = 1e-4)
})

# QTT mode (pre_copula = "long", the default): only e=0 is valid for a 4-period
# dataset with group=3, so aggregation is over a single cell. pte_qtt class,
# overall curve has one row per prob, all CIs valid, median QTT near 1.
test_that("panel_qtt panel QTT returns valid pte_qtt with correct structure", {
  dat <- make_panel_qtt_data()
  set.seed(42)
  res <- panel_qtt(yname = "y", gname = "group", tname = "period", idname = "id",
                   data = dat, biters = 20, gt_type = "qtt", probs = c(0.25, 0.5, 0.75))
  expect_s3_class(res, "pte_qtt")
  expect_true(all(c("probs", "qtt", "se", "lower_pw", "upper_pw") %in% names(res$overall)))
  expect_equal(nrow(res$overall), 3L)
  expect_false(anyNA(res$overall$qtt))
  expect_true(all(res$overall$lower_pw < res$overall$upper_pw))
  # exact CDF mixture (#19); 0.9732522 on the old 1,000-point grid
  expect_equal(res$overall$qtt[res$overall$probs == 0.5], 0.9690978, tolerance = 1e-4)
})

# pre_copula = "long" with unequally spaced calendar years. With periods
# 1974/1975/1978 and g = 1978, the only cell is e = 0, where the long and
# short copulas coincide (pre2 = 1974, pre1 = 1975). panel_qtt_long_agg()
# used to check pre2 in calendar units (2 * 1978 - 1978 - 2 = 1976, not in the
# data) and errored with "no valid post-treatment cells". Point estimates must
# not depend on whether time is coded in calendar years or as 1..T.
test_that("panel_qtt long-copula QTT works with unequally spaced years", {
  d <- lalonde.psid.panel
  d$g <- ifelse(d$treat == 1, 1978, 0)
  d$year_idx <- match(d$year, c(1974, 1975, 1978))
  d$g_idx <- ifelse(d$treat == 1, 3, 0)
  fit <- function(gname, tname, pre_copula) {
    set.seed(1)
    suppressWarnings(panel_qtt(yname = "re", gname = gname, tname = tname,
                               idname = "id", data = d, biters = 2,
                               gt_type = "qtt", pre_copula = pre_copula,
                               probs = c(0.5, 0.7, 0.9)))
  }
  long_years <- fit("g", "year", "long")
  expect_s3_class(long_years, "pte_qtt")
  expect_false(anyNA(long_years$overall$qtt))
  expect_equal(long_years$overall$qtt, fit("g_idx", "year_idx", "long")$overall$qtt)
  expect_equal(long_years$overall$qtt, fit("g", "year", "short")$overall$qtt)
})

# Anticipation is in units of the time variable, as in did::att_gt() and
# ptetools: with unequally spaced years and anticipation = 2, the base period
# (pre1) is the last year at least 2 years before treatment, cells are
# post-treatment from g - 2 on, and not-yet-treated controls must not be
# treated or anticipating treatment by the cell's latest period.
test_that("three_period_subset() uses anticipation in calendar time", {
  set.seed(1)
  n <- 300
  yrs <- c(2001, 2002, 2004, 2007, 2011)
  grp <- sample(c(0, 2007, 2011), n, replace = TRUE)
  d <- data.frame(id = rep(seq_len(n), each = 5), year = rep(yrs, n),
                  grp = rep(grp, each = 5))
  d$y <- rnorm(nrow(d))
  ptep <- suppressWarnings(ptetools::setup_pte(
    yname = "y", gname = "grp", tname = "year", idname = "id", data = d,
    anticipation = 2, required_pre_periods = 2
  ))
  idx <- function(year) match(year, yrs)
  cell <- function(g, tp, pre_copula) {
    gt <- three_period_subset(ptep$data, idx(g), idx(tp), anticipation = 2,
                              pre_copula = pre_copula)$gt_data
    list(periods = yrs[sort(unique(gt$.period))],
         controls = sort(unique(gt$grp[gt$D == 0])))
  }
  # 2007 cohort: base period 2004; 2007 and 2011 are post-treatment
  expect_equal(cell(2007, 2007, "long")$periods, c(2002, 2004, 2007))
  expect_equal(cell(2007, 2011, "long")$periods, c(2001, 2004, 2011))
  expect_equal(cell(2007, 2011, "short")$periods, c(2002, 2004, 2011))
  # the 2011 cohort is a valid control only through 2007 + 2 < 2011
  expect_equal(cell(2007, 2007, "long")$controls, c(0, 2011))
  expect_equal(cell(2007, 2011, "long")$controls, 0)
  # 2011 cohort: 2007 is still pre-treatment (2007 < 2011 - 2)
  expect_equal(cell(2011, 2007, "long")$periods, c(2002, 2004, 2007))
  expect_equal(cell(2011, 2011, "long")$periods, c(2004, 2007, 2011))
})
