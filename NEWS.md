# evmissing 1.0.1

The research paper Simpson, E. J. and Northrop, P. J. (2026) Accounting for missing data when modelling block maxima will appear in [Environmetrics in 2026](https://onlinelibrary.wiley.com/loi/1099095x/year/2026).

## Bug fixes and minor improvements

* When the argument `epsilon` to `confint.evmissing()` or `confint.return_level()` is negative, monotonic cubic spline interpolation is now used instead of quadratic interpolation. Quadratic interpolation could fail in some cases.
* Fixed bugs in `evmissing-internal.R` to avoid storing some slightly incorrect parameter values near the confidence limits in the attribute `"for_plot"` in the returned object.
