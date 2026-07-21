# qte: A package for computing quantile treatment effects

Provides several methods for computing the Quantile Treatment Effect
(QTE) and Quantile Treatment Effect on the Treated (QTT). The main cases
covered are (i) treatment is randomly assigned, (ii) treatment is as
good as randomly assigned after conditioning on covariates (selection on
observables) using the methods of Firpo (2007)
[doi:10.1111/j.1468-0262.2007.00738.x](https://doi.org/10.1111/j.1468-0262.2007.00738.x)
, and (iii) identification is based on a Difference in Differences
assumption, with support for several varieties including Athey and
Imbens (2006)
[doi:10.1111/j.1468-0262.2006.00668.x](https://doi.org/10.1111/j.1468-0262.2006.00668.x)
, Callaway and Li (2019)
[doi:10.3982/QE935](https://doi.org/10.3982/QE935) , and Callaway, Li,
and Oka (2018)
[doi:10.1016/j.jeconom.2018.06.008](https://doi.org/10.1016/j.jeconom.2018.06.008)
. Version 2.0 adds a unified staggered treatment adoption API (built on
'ptetools') for all DiD-based estimators, as well as a new
lagged-outcome unconfoundedness estimator ('lou_qtt').

## See also

Useful links:

- <https://bcallaway11.github.io/qte/>

- <https://github.com/bcallaway11/qte>

- Report bugs at <https://github.com/bcallaway11/qte/issues>

## Author

**Maintainer**: Brantly Callaway <brantly.callaway@uga.edu>

Authors:

- Brantly Callaway <brantly.callaway@uga.edu>
