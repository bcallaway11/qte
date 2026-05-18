# =============================================================================
# Title: Package startup and global variable declarations
# Description: Declares global variable names used by setupData via
#   parent.frame() assignment, suppressing R CMD CHECK NOTEs.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

# Variables injected into caller environments by setupData() via assign()
utils::globalVariables(c("yname", "treat", "panel", "treated.t", "treated.tmin1", "untreated.t", "untreated.tmin1", "F.treated.t", "F.treated.tmin1", "F.untreated.t", "F.untreated.tmin1", "xformla", "data", "x", "wname", "probs", "method", "eachIter", "F.untreated.tmin2", "F.treated.tmin2", "untreated.tmin2", "treated.tmin2", "untreated.change.t", "F.untreated.change.t", "untreated.change.tmin1", "F.untreated.change.tmin1", "treated.change.tmin1", "F.treated.change.tmin1", "tmin1", "tmin2", "idname"))
