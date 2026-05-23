# =============================================================================
# Title: Panel QTT (Callaway-Li 2019)
# Description: Implements the Callaway and Li (2019) panel QTT estimator.
#   three_period_subset is the custom ptetools subset function that pulls three
#   periods (pre2, pre1, post) for each (g,t) cell. panel_qtt_gt is the
#   workhorse attgt_fun implementing the Rosenblatt + copula-stability
#   identification. panel_qtt_long_agg is the custom aggregation function for
#   pre_copula = "long" (two-step weight renormalization over valid cells).
#   panel_qtt is the user-facing wrapper.
#   All three attgt_fun/subset_fun functions follow the ptetools patterns used
#   by cic_gt, qdid_gt, mdid_gt, and ddid_gt.
# Author: Brant Callaway
# Last update: 2026-05-21
# Date created: 2026-05-20
# =============================================================================

# --- 3-period subset function -------------------------------------------------

#' @title Three-Period Subset for Panel QTT
#'
#' @description Custom \code{subset_fun} for \code{ptetools::pte} that pulls
#'   three time periods for each (g,t) cell: pre2, pre1 = g - anticipation - 1,
#'   and post = tp. How pre2 is chosen depends on \code{pre_copula}:
#'
#'   \describe{
#'     \item{\code{"long"} (default)}{pre2 = 2*(g - anticipation) - tp - 2, so
#'       the copula-transfer window matches the event horizon e = tp - (g -
#'       anticipation). Cells where this pre2 does not exist in the data return
#'       empty data and \code{panel_qtt_gt} returns \code{NA}.}
#'     \item{\code{"short"}}{pre2 = g - anticipation - 2 (fixed at treatment
#'       onset regardless of e).}
#'   }
#'
#'   For pre-treatment placebo cells (\code{tp < g - anticipation}), the base
#'   always slides backwards regardless of \code{pre_copula}:
#'   pre2 = tp - 2, pre1 = tp - 1.  \code{panel_qtt()} sets
#'   \code{required_pre_periods = 2} so that the minimum tp in \code{tlist}
#'   is the third period, guaranteeing tp - 2 is always a valid period.
#'
#' @param data A data frame prepared by \code{ptetools::setup_pte}.
#' @param g The current group (first treatment period).
#' @param tp The current time period.
#' @param control_group \code{"notyettreated"} (default) or
#'   \code{"nevertreated"}.
#' @param anticipation Number of periods of anticipation. Default \code{0}.
#' @param pre_copula \code{"long"} (default) or \code{"short"}. See above.
#' @param ... Additional arguments; silently ignored.
#'
#' @return A list with elements \code{gt_data}, \code{n1}, \code{disidx}.
#'
#' @keywords internal
three_period_subset <- function(data,
                                g,
                                tp,
                                control_group = "notyettreated",
                                anticipation  = 0,
                                pre_copula    = "long",
                                ...) {
  effective_g <- g - anticipation

  if (tp >= effective_g) {
    # post-treatment cell: pre1 fixed at treatment onset; pre2 depends on pre_copula
    pre1 <- effective_g - 1L
    pre2 <- if (pre_copula == "long") {
      # extend copula window to match event horizon e = tp - effective_g
      2L * effective_g - tp - 2L
    } else {
      # anchor both base periods at treatment onset
      effective_g - 2L
    }
  } else {
    # pre-test placebo cell: slide base backwards regardless of pre_copula
    pre1 <- tp - 1L
    pre2 <- tp - 2L
  }

  if (control_group == "notyettreated") {
    this.data <- subset(data, G == g | G > tp | G == 0)
  } else {
    this.data <- subset(data, G == g | G == 0)
  }

  this.data <- subset(this.data, period %in% c(pre2, pre1, tp))

  this.data$name <- ifelse(this.data$period == tp,  "post",
                    ifelse(this.data$period == pre1, "pre1", "pre2"))
  this.data$D    <- 1L * (this.data$G == g)
  this.data      <- ptetools::gt_data_frame(this.data)

  n1     <- length(unique(this.data$id))
  disidx <- unique(data$id) %in% unique(this.data$id)

  list(gt_data = this.data, n1 = n1, disidx = disidx)
}


# --- 3-period workhorse -------------------------------------------------------

#' @title Panel QTT: group-time estimator (Callaway-Li 2019)
#'
#' @description Computes the panel QTT and counterfactual outcome distribution
#'   for a single (g,t) cell using the three-period copula-stability estimator
#'   of Callaway and Li (2019). Serves directly as the \code{attgt_fun}
#'   argument to \code{ptetools::pte}.
#'
#'   \strong{Identification.} Under copula stability
#'   \deqn{C(F_{Y^{\text{pre2}}|D=1},\, F_{\Delta Y^{\text{pre}}|D=1}) =
#'         C(F_{Y^{\text{pre1}}(0)|D=1},\, F_{\Delta Y^{\text{post}}(0)|D=1})}
#'   and distributional parallel trends on changes
#'   \deqn{F_{\Delta Y^{\text{post}}(0)|D=1} = F_{\Delta Y^{\text{ctrl}}}}
#'   the counterfactual outcome for each treated unit \eqn{i} is
#'   \eqn{kcf_i = L_i + C_i,}
#'   where \eqn{L_i = Q_{Y^{\text{pre1}}|D=1}(u_i)} is the treated pre1
#'   quantile at rank \eqn{u_i = F_{Y^{\text{pre2}}|D=1}(Y^{\text{pre2}}_i)}
#'   (the Rosenblatt transform), and
#'   \eqn{C_i = Q_{\Delta Y^{\text{ctrl}}}(v_i)} is the control change at rank
#'   \eqn{v_i = F_{\Delta Y^{\text{pre}}|D=1}(\Delta Y^{\text{pre}}_i)}.
#'   The goal is distributional: \eqn{\{kcf_i\}} is a sample from
#'   \eqn{F_{Y(0)^{\text{post}}|D=1}}, not individual counterfactuals.
#'
#'   \strong{Panel data required for both groups.} Three periods for treated
#'   (pre2, pre1, post) and two for control (pre1, post).
#'
#' @param gt_data A \code{gt_data_frame} with columns \code{name}
#'   (\code{"pre2"}, \code{"pre1"}, \code{"post"}), \code{D}, \code{Y},
#'   \code{id}, \code{.w}, and any covariate columns.
#' @param xformula One-sided formula for covariates. Default \code{~1}.
#' @param ... Additional arguments passed through by \code{ptetools}.
#'
#' @return A \code{ptetools::attgt_noif} object with \code{att} and, in
#'   \code{extra_gt_returns}, \code{F0} (counterfactual ECDF), \code{F1}
#'   (observed post ECDF), and \code{Fte} (individual-effect ECDF).
#'
#' @references
#'   Callaway, Brantly and Tong Li. ``Quantile Treatment Effects in
#'   Difference-in-Differences Models with Panel Data.'' Quantitative
#'   Economics 10(4), pp. 1579-1618, 2019.
#'
#' @importFrom quantreg rq
#'
#' @export
panel_qtt_gt <- function(gt_data, xformula = ~1, ...) {
  gt_data <- droplevels(gt_data)

  # extract the five cells; pre2_ctrl is not needed for identification
  pre2_trt  <- gt_data[gt_data$name == "pre2" & gt_data$D == 1, ]
  pre1_trt  <- gt_data[gt_data$name == "pre1" & gt_data$D == 1, ]
  post_trt  <- gt_data[gt_data$name == "post" & gt_data$D == 1, ]
  pre1_ctrl <- gt_data[gt_data$name == "pre1" & gt_data$D == 0, ]
  post_ctrl <- gt_data[gt_data$name == "post" & gt_data$D == 0, ]

  # guard: all five cells must be non-empty
  if (nrow(pre2_trt)  == 0L || nrow(pre1_trt)  == 0L ||
      nrow(post_trt)  == 0L || nrow(pre1_ctrl) == 0L ||
      nrow(post_ctrl) == 0L) {
    return(ptetools::attgt_noif(attgt = NA))
  }

  # sort by id — panel structure enforced by panel_qtt() via panel = TRUE
  pre2_trt  <- pre2_trt[order(pre2_trt$id), ]
  pre1_trt  <- pre1_trt[order(pre1_trt$id), ]
  post_trt  <- post_trt[order(post_trt$id), ]
  pre1_ctrl <- pre1_ctrl[order(pre1_ctrl$id), ]
  post_ctrl <- post_ctrl[order(post_ctrl$id), ]

  # uppercase Y objects follow econometric notation # nolint: object_name_linter
  Y_pre2_trt  <- pre2_trt$Y  # nolint: object_name_linter
  Y_pre1_trt  <- pre1_trt$Y  # nolint: object_name_linter
  Y_post_trt  <- post_trt$Y  # nolint: object_name_linter
  Y_pre1_ctrl <- pre1_ctrl$Y # nolint: object_name_linter
  Y_post_ctrl <- post_ctrl$Y # nolint: object_name_linter

  # normalise sampling weights within each (period x treatment) cell
  w_pre2_trt  <- pre2_trt$.w  / sum(pre2_trt$.w)  # nolint: object_name_linter
  w_pre1_trt  <- pre1_trt$.w  / sum(pre1_trt$.w)  # nolint: object_name_linter
  w_post_trt  <- post_trt$.w  / sum(post_trt$.w)  # nolint: object_name_linter
  w_pre1_ctrl <- pre1_ctrl$.w / sum(pre1_ctrl$.w) # nolint: object_name_linter

  # observed changes: treated pre-period change and ctrl post-pre1 change
  dY_trt  <- Y_pre1_trt - Y_pre2_trt   # nolint: object_name_linter
  dY_ctrl <- Y_post_ctrl - Y_pre1_ctrl  # nolint: object_name_linter

  # Helper: weighted type-1 quantile (smallest y with weighted CDF >= p).
  wquant <- function(y_vec, w_vec, probs) {
    ord <- order(y_vec)
    y_s <- y_vec[ord]
    cdf <- cumsum(w_vec[ord])
    y_s[vapply(probs, function(p) {
      idx <- which(cdf >= p)
      if (length(idx) == 0L) length(y_s) else idx[1L]
    }, integer(1L))]
  }

  # No-covariate path:
  #
  # Step 1 (Rosenblatt): transfer each treated unit's pre2 rank into a pre1 level.
  #   u_i = F_{Y^pre2|D=1}(Y^pre2_i)  — rank in the pre2 treated marginal
  #   L_i = Q_{Y^pre1|D=1}(u_i)       — pre1 level at rank u_i
  u <- vapply(Y_pre2_trt,
              function(y) sum(w_pre2_trt * (Y_pre2_trt <= y)),
              numeric(1L))
  L <- wquant(Y_pre1_trt, w_pre1_trt, u) # nolint: object_name_linter

  # Step 2 (copula stability + distributional PT on changes):
  #   v_i = F_{dY^pre|D=1}(dY^pre_i)  — rank of pre-treatment change
  #   C_i = Q_{dY^ctrl}(v_i)          — ctrl change at rank v_i
  v <- vapply(dY_trt,
              function(dy) sum(w_pre2_trt * (dY_trt <= dy)),
              numeric(1L))
  C <- wquant(dY_ctrl, w_pre1_ctrl, v) # nolint: object_name_linter

  kcf <- L + C

  att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(kcf, w_post_trt)
  F0  <- BMisc::weighted_ecdf(kcf,        weights = w_post_trt) # nolint: object_name_linter
  F1  <- BMisc::weighted_ecdf(Y_post_trt, weights = w_post_trt) # nolint: object_name_linter
  # Fte identified because kcf is unit-indexed (treated group)
  Fte <- BMisc::weighted_ecdf(Y_post_trt - kcf, weights = w_post_trt) # nolint: object_name_linter

  # Covariate adjustment: four QR models replace the unconditional rank transforms.
  # Sampling weights are passed to rq() throughout.
  if (length(BMisc::rhs_vars(xformula)) > 0) {
    u_seq   <- seq(0.01, 0.99, 0.01)
    n1      <- nrow(pre2_trt)
    yformla <- BMisc::toformula("Y", BMisc::rhs_vars(xformula)) # nolint: object_name_linter

    # --- Rosenblatt step (conditioned on pre2 covariates throughout) ---
    # u_i = F(Y^pre2_i | X_pre2_i, D=1) via QR on pre2 treated
    # L_i = Q(Y^pre1 | u_i, X_pre2_i, D=1) via QR on pre1 treated,
    #        predicted at pre2 covariates (copula transfers conditional rank
    #        across pre2 → pre1 marginals holding X fixed)
    QR_Y_pre2_trt <- suppressWarnings(rq(yformla, data = pre2_trt, tau = u_seq, weights = .w)) # nolint: object_name_linter
    QR_Y_pre1_trt <- suppressWarnings(rq(yformla, data = pre1_trt, tau = u_seq, weights = .w)) # nolint: object_name_linter

    QR_Y_pre2F <- predict(QR_Y_pre2_trt, newdata = pre2_trt, type = "Fhat", stepfun = TRUE) # nolint: object_name_linter
    u <- sapply(seq_len(n1), function(i) QR_Y_pre2F[[i]](pre2_trt$Y[i]))

    QR_Y_pre1Q <- predict(QR_Y_pre1_trt, newdata = pre2_trt, type = "Qhat", stepfun = TRUE) # nolint: object_name_linter
    L <- sapply(seq_len(n1), function(i) QR_Y_pre1Q[[i]](u[i])) # nolint: object_name_linter

    # --- Copula-stability step ---
    # QR_dY_trt: fit on treated units with pre2 covariates and dY_trt outcome
    # v_i = F(dY_trt_i | X_pre2_i, D=1)
    # QR_dY_ctrl: fit on ctrl units with pre1 covariates and dY_ctrl outcome
    # C_i = Q(dY_ctrl | v_i, X_pre1_i, D=1) evaluated at treated pre1 X
    dY_trt_data    <- pre2_trt   # nolint: object_name_linter
    dY_trt_data$Y  <- dY_trt
    dY_ctrl_data   <- pre1_ctrl  # nolint: object_name_linter
    dY_ctrl_data$Y <- dY_ctrl

    QR_dY_trt  <- suppressWarnings(rq(yformla, data = dY_trt_data,  tau = u_seq, weights = .w)) # nolint: object_name_linter
    QR_dY_ctrl <- suppressWarnings(rq(yformla, data = dY_ctrl_data, tau = u_seq, weights = .w)) # nolint: object_name_linter

    QR_dY_trtF <- predict(QR_dY_trt, newdata = dY_trt_data, type = "Fhat", stepfun = TRUE) # nolint: object_name_linter
    v <- sapply(seq_len(n1), function(i) QR_dY_trtF[[i]](dY_trt[i]))

    # QR_dY_ctrl evaluated at treated pre1 covariates (matching ctrl model's X)
    QR_dY_ctrlQ <- predict(QR_dY_ctrl, newdata = pre1_trt, type = "Qhat", stepfun = TRUE) # nolint: object_name_linter
    C <- sapply(seq_len(n1), function(i) QR_dY_ctrlQ[[i]](v[i])) # nolint: object_name_linter

    kcf <- L + C

    att <- weighted.mean(Y_post_trt, w_post_trt) - weighted.mean(kcf, w_post_trt)
    F0  <- BMisc::weighted_ecdf(kcf,             weights = w_post_trt) # nolint: object_name_linter
    Fte <- BMisc::weighted_ecdf(Y_post_trt - kcf, weights = w_post_trt) # nolint: object_name_linter
  }

  ptetools::attgt_noif(attgt = att, extra_gt_returns = list(F0 = F0, F1 = F1, Fte = Fte))
}


# --- Aggregation function for pre_copula = "long" ----------------------------

#' @title QTT Aggregation for pre_copula = "long"
#'
#' @description Custom \code{aggregation_fun} for \code{ptetools::pte} that
#'   handles the two-step weight renormalization needed when
#'   \code{pre_copula = "long"}: some (g,t) cells are structurally
#'   uncomputable (their pre2 does not exist in the data), and those cells
#'   must be excluded from aggregation with proper weight redistribution.
#'
#'   \strong{Two-step weights (overall):}
#'   \enumerate{
#'     \item Identify surviving groups (those with at least one valid
#'       post-treatment cell); renormalize \eqn{p_g} among those groups only.
#'     \item Within each surviving group, assign equal weight
#'       \eqn{p_g / |\text{valid cells for } g|} to each valid post-treatment
#'       cell.
#'   }
#'
#'   \strong{Dynamic (event-study) weights:} for each event time e, renormalize
#'   \eqn{p_g} among groups where (g, g+e) is valid.
#'
#'   \strong{Error behaviour:} if a structurally valid cell returns \code{NA}
#'   att, the function errors immediately — such NAs indicate a bug in
#'   \code{panel_qtt_gt}, not a structural limitation.
#'
#' @param attgt.list List of per-cell attgt results from \code{compute.pte}.
#' @param ptep ptetools parameter environment.
#' @param extra_gt_returns List of per-cell extra returns (F0, F1, Fte ECDFs).
#'
#' @return A list with \code{overall_results}, \code{dyn_results},
#'   \code{group_results}, \code{attgt_results}, \code{F0_overall},
#'   \code{F1_overall} — matching the structure returned by
#'   \code{ptetools::qtt_pte_aggregations}.
#'
#' @keywords internal
panel_qtt_long_agg <- function(attgt.list, ptep, extra_gt_returns) {
  probs        <- if (is.null(ptep$probs)) seq(0.05, 0.95, 0.05) else ptep$probs
  anticipation <- if (!is.null(ptep$anticipation)) ptep$anticipation else 0L
  data         <- ptep$data
  yname        <- ptep$yname

  # --- resolve to original time periods (mirrors attgt_pte_aggregations) ---
  time.periods <- ptep$tlist
  groups       <- ptep$glist
  original_tp  <- sort(unique(data[, ptep$tname]))

  # must capture the flag before converting time.periods, so the same
  # conversion is applied consistently to cell_g/cell_tp below
  needs_t2orig <- !all(time.periods %in% original_tp)
  if (needs_t2orig) {
    time.periods <- sapply(time.periods, BMisc::t2orig, original_tp)
    groups       <- sapply(groups,       BMisc::t2orig, original_tp)
  }

  # --- per-cell metadata ---
  cell_g   <- vapply(attgt.list, function(r) r$group,       numeric(1L))
  cell_tp  <- vapply(attgt.list, function(r) r$time.period, numeric(1L))
  cell_att <- vapply(attgt.list, function(r) r$att,         numeric(1L))

  if (needs_t2orig) {
    cell_g  <- sapply(cell_g,  BMisc::t2orig, original_tp)
    cell_tp <- sapply(cell_tp, BMisc::t2orig, original_tp)
  }

  effective_g <- cell_g - anticipation
  cell_e      <- cell_tp - effective_g
  is_post     <- cell_e >= 0L

  # --- validity mask ---
  # post-treatment: valid iff pre2 = 2*effective_g - tp - 2 is in original_tp
  # pre-test: always valid (sliding base, guaranteed by required_pre_periods = 2)
  pre2_long <- 2L * effective_g - cell_tp - 2L
  is_valid  <- ifelse(is_post, pre2_long %in% original_tp, TRUE)

  # --- error on valid cells with unexpected NA att ---
  if (any(is_valid & is.na(cell_att))) {
    n_bad <- sum(is_valid & is.na(cell_att))
    stop(sprintf(
      paste0("panel_qtt_long_agg: %d structurally valid cell(s) returned NA att;",
             " check panel_qtt_gt for bugs"),
      n_bad
    ))
  }

  # --- group sizes (count treated units at first time period) ---
  n_group <- sapply(groups, function(gg) {
    nrow(subset(data, data[, ptep$gname] == gg &
                      data[, ptep$tname] == original_tp[1L]))
  })
  names(n_group) <- as.character(groups)

  # --- F0/F1 for all cells (NULL for invalid/NA cells) ---
  F0_gt <- lapply(extra_gt_returns, function(egr) egr$extra_gt_returns$F0) # nolint: object_name_linter
  F1_gt <- lapply(extra_gt_returns, function(egr) egr$extra_gt_returns$F1) # nolint: object_name_linter

  y.seq <- quantile(data[, yname], probs = seq(0, 1, length.out = 1000))

  # --- overall QTT (two-step weights over valid post-treatment cells) ---
  valid_post <- is_valid & is_post
  surviving_g <- unique(cell_g[valid_post])

  if (length(surviving_g) == 0L) {
    stop("panel_qtt_long_agg: no valid post-treatment cells; cannot compute overall QTT")
  }

  # step 1: renormalize pg among surviving groups
  pg <- n_group[as.character(surviving_g)]
  pg <- pg / sum(pg)

  # step 2: within-group equal weight per valid post-treatment cell
  n_valid_post_g <- vapply(as.character(surviving_g),
                           function(g) sum(valid_post & as.character(cell_g) == g),
                           integer(1L))

  overall_w <- rep(0, length(attgt.list))
  for (i in seq_along(surviving_g)) {
    g   <- surviving_g[i]
    idx <- which(valid_post & cell_g == g)
    overall_w[idx] <- pg[i] / n_valid_post_g[i]
  }

  nonzero    <- which(overall_w > 0)
  F0_overall <- BMisc::combine_ecdfs(y.seq, ecdflist = F0_gt[nonzero], weights = overall_w[nonzero]) # nolint: object_name_linter
  F1_overall <- BMisc::combine_ecdfs(y.seq, ecdflist = F1_gt[nonzero], weights = overall_w[nonzero]) # nolint: object_name_linter
  overall_qtt <- quantile(F1_overall, probs = probs, type = 1) -
                 quantile(F0_overall, probs = probs, type = 1)
  overall_results <- data.frame(probs = probs, qtt = overall_qtt)

  # --- dynamic QTT (renormalize pg by group size at each e) ---
  eseq <- sort(unique(cell_e[is_valid]))
  dyn_results_list <- lapply(eseq, function(this_e) {
    e_idx <- which(is_valid & cell_e == this_e)
    if (length(e_idx) == 0L) return(NULL)
    g_at_e <- cell_g[e_idx]
    pg_e   <- n_group[as.character(g_at_e)]
    pg_e   <- pg_e / sum(pg_e)
    F0_e <- BMisc::combine_ecdfs(y.seq, ecdflist = F0_gt[e_idx], weights = pg_e) # nolint: object_name_linter
    F1_e <- BMisc::combine_ecdfs(y.seq, ecdflist = F1_gt[e_idx], weights = pg_e) # nolint: object_name_linter
    qtt_e <- quantile(F1_e, probs = probs, type = 1) -
              quantile(F0_e, probs = probs, type = 1)
    data.frame(e = this_e, probs = probs, qtt = qtt_e)
  })
  dyn_results <- do.call(rbind, Filter(Negate(is.null), dyn_results_list))

  # --- group QTT (uniform weight over valid post-treatment cells per group) ---
  group_results_list <- lapply(surviving_g, function(g) {
    g_idx <- which(valid_post & cell_g == g)
    if (length(g_idx) == 0L) return(NULL)
    w_g <- rep(1 / length(g_idx), length(g_idx))
    F0_g <- BMisc::combine_ecdfs(y.seq, ecdflist = F0_gt[g_idx], weights = w_g) # nolint: object_name_linter
    F1_g <- BMisc::combine_ecdfs(y.seq, ecdflist = F1_gt[g_idx], weights = w_g) # nolint: object_name_linter
    qtt_g <- quantile(F1_g, probs = probs, type = 1) -
              quantile(F0_g, probs = probs, type = 1)
    data.frame(group = g, probs = probs, qtt = qtt_g)
  })
  group_results <- do.call(rbind, Filter(Negate(is.null), group_results_list))

  # --- attgt_results (scalar ATTs for valid cells; for SE tables) ---
  valid_idx <- which(is_valid)
  attgt_results <- data.frame(
    group       = cell_g[valid_idx],
    time.period = cell_tp[valid_idx],
    att         = cell_att[valid_idx]
  )

  list(
    overall_results = overall_results,
    dyn_results     = dyn_results,
    group_results   = group_results,
    attgt_results   = attgt_results,
    F0_overall      = F0_overall, # nolint: object_name_linter
    F1_overall      = F1_overall  # nolint: object_name_linter
  )
}


# --- Unified user-facing function ---------------------------------------------

#' @title Panel QTT (Callaway-Li 2019)
#'
#' @description Computes Quantile Treatment Effects on the Treated (QTT) and
#'   the Average Treatment Effect on the Treated (ATT) using the three-period
#'   panel identification strategy of Callaway and Li (2019). Handles
#'   two-period (one post + two pre-treatment periods) and staggered treatment
#'   adoption via \code{ptetools}. Requires panel data for both groups.
#'
#'   \strong{Identification.} Under a copula stability assumption and
#'   distributional parallel trends on changes, the counterfactual distribution
#'   \eqn{F_{Y(0)^{\text{post}}|D=1}} is recovered from three periods of panel
#'   data without requiring rank invariance.
#'
#'   \strong{pre_copula.} Controls which pre-treatment periods are used as the
#'   copula-transfer base for post-treatment cells (g, g+e):
#'   \describe{
#'     \item{\code{"long"} (default)}{pre2 = 2*(g - anticipation) - (g+e) - 2,
#'       matching the copula window to the event horizon e. Cells where pre2
#'       does not exist in the data are excluded from aggregation with two-step
#'       weight renormalization.}
#'     \item{\code{"short"}}{pre2 = g - anticipation - 2, always anchored at
#'       treatment onset.}
#'   }
#'   Pre-test placebo cells always use the sliding base regardless of
#'   \code{pre_copula}.
#'
#' @param yname Name of the outcome variable.
#' @param gname Name of the treatment group variable (first treatment period;
#'   0 for never-treated units).
#' @param tname Name of the time period variable.
#' @param idname Name of the unit id variable (required).
#' @param data A data frame.
#' @param xformula One-sided formula for covariates. Default \code{~1}.
#' @param weightsname Name of the sampling weights column. Default \code{NULL}.
#' @param control_group \code{"notyettreated"} (default) or
#'   \code{"nevertreated"}.
#' @param anticipation Number of periods of anticipation. Default \code{0}.
#' @param alp Significance level. Default \code{0.05}.
#' @param cband Logical; uniform confidence band if \code{TRUE} (default).
#' @param biters Number of bootstrap iterations. Default \code{100}.
#' @param cl Number of parallel clusters. Default \code{1}.
#' @param gt_type \code{"att"} (default) or \code{"qtt"}.
#' @param probs Quantile grid for \code{gt_type = "qtt"}. Default
#'   \code{seq(0.05, 0.95, 0.05)}.
#' @param pre_copula \code{"long"} (default) or \code{"short"}. See Details.
#'
#' @return For \code{gt_type = "att"}, a \code{pte_emp_boot} object. For
#'   \code{gt_type = "qtt"}, a \code{pte_qtt} object.
#'
#' @references
#'   Callaway, Brantly and Tong Li. ``Quantile Treatment Effects in
#'   Difference-in-Differences Models with Panel Data.'' Quantitative
#'   Economics 10(4), pp. 1579-1618, 2019.
#'
#' @seealso \code{\link{cic}}, \code{\link{qdid}}, \code{\link{mdid}},
#'   \code{\link{ddid}}
#'
#' @examples
#' \donttest{
#' data(mpdta, package = "did")
#'
#' ## Panel QTT with rolling pre-period copula (default pre_copula = "long")
#' res <- panel_qtt(yname = "lemp", gname = "first.treat", tname = "year",
#'                  idname = "countyreal", data = mpdta,
#'                  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
#' summary(res)
#' }
#'
#' @export
panel_qtt <- function(yname,
                      gname,
                      tname,
                      idname,
                      data,
                      xformula      = ~1,
                      weightsname   = NULL,
                      control_group = "notyettreated",
                      anticipation  = 0,
                      alp           = 0.05,
                      cband         = TRUE,
                      biters        = 100,
                      cl            = 1,
                      gt_type       = "att",
                      probs         = NULL,
                      pre_copula    = "long") {

  pre_copula <- match.arg(pre_copula, c("long", "short"))

  # Select aggregation_fun based on gt_type and pre_copula.
  # For ATT mode, attgt_pte_aggregations already drops NA cells (complete.cases)
  # and applies group-size weights over surviving groups — correct two-step
  # behaviour for both pre_copula variants.
  # For QTT mode with "long", custom two-step aggregation is required because
  # qtt_pte_aggregations does not filter invalid cells before combine_ecdfs.
  aggregation_fun <- if (gt_type == "qtt" && pre_copula == "long") {
    panel_qtt_long_agg
  } else if (gt_type == "qtt") {
    ptetools::qtt_pte_aggregations
  } else {
    # wrap to match the 3-argument aggregation_fun signature
    function(al, p, eg) ptetools::attgt_pte_aggregations(al, p)
  }

  # required_pre_periods = 2: drop first two periods from tlist so every tp
  # satisfies tp >= 3, guaranteeing the sliding base (pre2 = tp-2) is valid
  # for pre-test cells.
  ptetools::pte(
    yname            = yname,
    gname            = gname,
    tname            = tname,
    idname           = idname,
    data             = data,
    panel            = TRUE,
    setup_pte_fun    = ptetools::setup_pte,
    subset_fun       = three_period_subset,
    attgt_fun        = panel_qtt_gt,
    aggregation_fun  = aggregation_fun,
    xformula         = xformula,
    weightsname      = weightsname,
    control_group    = control_group,
    anticipation     = anticipation,
    cband            = cband,
    alp              = alp,
    boot_type        = "empirical",
    biters           = biters,
    cl               = cl,
    gt_type          = gt_type,
    probs            = probs,
    pre_copula       = pre_copula,
    required_pre_periods = 2L
  )
}
