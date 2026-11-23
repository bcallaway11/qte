-- =============================================================================
-- Title: Quarto filter for figure alt text in the gfm README
-- Description: Quarto (1.6) writes a chunk's `fig-alt` as `data-fig-alt` in
--   gfm output, which GitHub and pkgdown do not read as alt text (pkgdown
--   reports "Missing alt-text in 'README.md'"). This filter, run after
--   Quarto's own filters (see README.qmd), turns it into a real `alt`
--   attribute.
-- Author: Brant Callaway
-- Last update: 2026-11-23
-- Date created: 2026-11-23
-- =============================================================================

function Image(img)
  local alt = img.attributes["fig-alt"]
  if alt then
    img.attributes["fig-alt"] = nil
    img.attributes["alt"] = alt
  end
  return img
end
