# =============================================================================
# Title: Build the pkgdown site
# Description: Wrapper around pkgdown::build_site() that keeps the LLM agent
#   instruction files (CLAUDE.md, AGENTS.md) off the public site. pkgdown
#   renders every top-level *.md file as a page (and indexes it in
#   sitemap.xml and search.json), with no option to exclude one, so the
#   files are moved aside for the build and restored afterwards. Deleting
#   the pages after the build is not enough: search.json would still hold
#   their full text. Also drops pkgdown's `qte.html` redirect, which
#   Dropbox renames to "qte (Case Conflict).html" because it differs from
#   `QTE.html` only in case (only that renamed copy is removed: deleting
#   "qte.html" would remove QTE.html on a case-insensitive file system).
#   Run from the package root: Rscript pkgdown/build_site.R
# Author: Brant Callaway
# Last update: 2026-11-23
# Date created: 2026-11-23
# =============================================================================

# Inside a function so on.exit() restores the files even if the build fails
# (on.exit() has no effect at the top level of a script).
build_site_hidden <- function(hidden = c("CLAUDE.md", "AGENTS.md")) {
  hidden <- hidden[file.exists(hidden)]
  aside <- paste0(hidden, ".pkgdown-hidden")
  on.exit(file.rename(aside[file.exists(aside)], hidden[file.exists(aside)]),
          add = TRUE)
  stopifnot(all(file.rename(hidden, aside)))

  # pages left over from earlier builds
  unlink(file.path("docs", c(paste0(sub("\\.md$", "", hidden), ".html"),
                             hidden)))

  pkgdown::build_site(preview = FALSE)

  unlink(file.path("docs", "reference", "qte (Case Conflict).html"))
}

build_site_hidden()
