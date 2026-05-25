# Package index

## Cross-sectional estimators

Estimators for the QTE and QTT under unconfoundedness (no panel data
required).

- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  : unc_qte

## Panel and repeated cross-section estimators

Estimators for the QTT using panel or repeated cross-section data. All
support staggered treatment adoption via ptetools.

- [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) : Change
  in Changes
- [`cic_gt()`](https://bcallaway11.github.io/qte/reference/cic_gt.md) :
  Change in Changes: group-time estimator
- [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md) :
  Quantile Difference-in-Differences
- [`qdid_gt()`](https://bcallaway11.github.io/qte/reference/qdid_gt.md)
  : Quantile Difference-in-Differences: group-time estimator
- [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  : Panel QTT (Callaway-Li 2019)
- [`panel_qtt_gt()`](https://bcallaway11.github.io/qte/reference/panel_qtt_gt.md)
  : Panel QTT: group-time estimator (Callaway-Li 2019)
- [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) :
  Distributional Difference-in-Differences
- [`ddid_gt()`](https://bcallaway11.github.io/qte/reference/ddid_gt.md)
  : Distributional DiD: group-time estimator
- [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) : Mean
  Difference-in-Differences
- [`mdid_gt()`](https://bcallaway11.github.io/qte/reference/mdid_gt.md)
  : Mean Difference-in-Differences: group-time estimator
- [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  : Lagged Outcome Unconfoundedness QTT

## Package overview

- [`qte`](https://bcallaway11.github.io/qte/reference/qte-package.md)
  [`qte-package`](https://bcallaway11.github.io/qte/reference/qte-package.md)
  : qte: A package for computing quantile treatment effects

## Output and plotting

The QTE S3 class returned by unc_qte(), and methods for printing,
summarising, and plotting.

- [`QTE()`](https://bcallaway11.github.io/qte/reference/QTE.md) : QTE
- [`summary(`*`<QTE>`*`)`](https://bcallaway11.github.io/qte/reference/summary.QTE.md)
  : Summary
- [`print(`*`<summary.QTE>`*`)`](https://bcallaway11.github.io/qte/reference/print.summary.QTE.md)
  : Print summary.QTE
- [`plot(`*`<QTE>`*`)`](https://bcallaway11.github.io/qte/reference/plot.QTE.md)
  : plot.QTE
- [`autoplot(`*`<QTE>`*`)`](https://bcallaway11.github.io/qte/reference/autoplot.QTE.md)
  : autoplot.QTE
- [`ggqte()`](https://bcallaway11.github.io/qte/reference/ggqte.md) :
  ggqte

## Data

Datasets bundled with the package.

- [`lalonde`](https://bcallaway11.github.io/qte/reference/lalonde.md) :
  Lalonde (1986)'s NSW Dataset
- [`lalonde.exp`](https://bcallaway11.github.io/qte/reference/lalonde.exp.md)
  : Lalonde's Experimental Dataset
- [`lalonde.exp.panel`](https://bcallaway11.github.io/qte/reference/lalonde.exp.panel.md)
  : Lalonde's Panel Experimental Dataset
- [`lalonde.psid`](https://bcallaway11.github.io/qte/reference/lalonde.psid.md)
  : Lalonde's Observational Dataset
- [`lalonde.psid.panel`](https://bcallaway11.github.io/qte/reference/lalonde.psid.panel.md)
  : Lalonde's Panel Observational Dataset
