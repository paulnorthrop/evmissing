
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Appveyor Build
status](https://ci.appveyor.com/api/projects/status/99jojhgk9t4agdmv/branch/main?svg=true)](https://ci.appveyor.com/project/paulnorthrop/evmiss/branch/main)
[![R-CMD-check](https://github.com/paulnorthrop/evmiss/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paulnorthrop/evmiss/actions/workflows/R-CMD-check.yaml)
[![Coverage
Status](https://codecov.io/github/paulnorthrop/evmiss/coverage.svg?branch=master)](https://app.codecov.io/github/paulnorthrop/evmiss?branch=master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/evmiss)](https://cran.r-project.org/package=evmiss)
[![Downloads
(monthly)](https://cranlogs.r-pkg.org/badges/evmiss?color=brightgreen)](https://cran.r-project.org/package=evmiss)
[![Downloads
(total)](https://cranlogs.r-pkg.org/badges/grand-total/evmiss?color=brightgreen)](https://cran.r-project.org/package=evmiss)

## Extreme Value Analyses with Missing Data

Performs likelihood-based extreme value inferences with adjustment for
the presence of missing values. A GEV distribution is fitted to block
maxima using maximum likelihood estimation, with the GEV location and
scale parameter reflecting the numbers of non-missing raw values in each
block. A Bayesian version is also provided. For the purposes of
comparison, there are options to make no adjustment for missing values
or to discard any block maximum for which greater than a percentage of
the underlying raw values are missing. A plot method provides a set of
standard model diagnostic plots, with appropriate adjustment made for
the presence of missing values. Example datasets containing missing
values are provided.

The `evmiss` package was created to accompany the research paper
[Simpson, E. S. and Northrop, P. J. (2025) Accounting for missing data
when modelling block maxima](https://arxiv.org/abs/2512.15429).

## An example

The main function in `evmiss` is `gev_mle()`, which fits a GEV
distribution to block maxima using maximum likelihood estimation, with
the option to make an adjustment for the numbers of non-missing raw
values in each block.

The main adjustment is based on the strong assumption that missing
values occur completely at random. We suppose that a block maximum $M_n$
based on a full block of length $n$ has a GEV($\mu$, $\sigma$, $\xi$)
distribution, with distribution function $G(x)$. Let $n_i$ be the number
of non-missing values in block $i$ and let $M_{n_i}$ denote the block
maximum of such a block. We suppose that $M_{n_i}$ has a GEV($\mu(n_i)$,
$\sigma(n_i)$, $\xi$) distribution, where
$\mu(n_i) = \mu + \sigma [(n_i/n)^\xi -1] / \xi$,
$\sigma(n_i) = \sigma (n_i/n)^\xi$.

### Sea surge height data

We illustrate this using the data `BrestSurgeMaxima`, which is a data
frame containing annual maximum sea surge heights (in cm) at Brest,
France for the years 1846-2007 inclusive.

``` r
library(evmiss)
head(BrestSurgeMaxima)
#>      maxima notNA   n block
#> 1846 59.987   361 365     1
#> 1847 58.873   344 365     2
#> 1848 59.749   366 366     3
#> 1849 49.547   365 365     4
#> 1850 55.512   365 365     5
#> 1851 69.422   365 365     6
```

In addition to the annual maxima this data frame includes `notNA`, the
number of days in each year for which raw data were available, and `n`,
the number of days in the year, that is, the block length. The function
`block_maxima()` can be used to create data of this format from a raw
time series and information about how blocks are defined. Alternatively,
we can provide a raw time series directly to `gev_mle()`.

### Model fitting

``` r
# Make the adjustment
fit_adjust <- gev_mle(BrestSurgeMaxima)
summary(fit_adjust)
#> 
#> Call:
#> gev_mle(data = BrestSurgeMaxima)
#> 
#>       Estimate Std. Error
#> mu    52.89000     1.0650
#> sigma 11.84000     0.7361
#> xi    -0.02375     0.0445

# Make no adjustment
fit_no_adjust <- gev_mle(BrestSurgeMaxima, adjust = FALSE)
summary(fit_no_adjust)
#> 
#> Call:
#> gev_mle(data = BrestSurgeMaxima, adjust = FALSE)
#> 
#>       Estimate Std. Error
#> mu    52.27000    1.07300
#> sigma 12.09000    0.76170
#> xi    -0.03005    0.04388
```

The most obvious difference between these model fits is that the
estimated location parameter with adjustment (52.27 cm) is smaller than
that when the adjustment is used (52.89 cm). This is as expected because
if some of the data that could contribute to an annual maximum are
missing then the observed annual maximum is stochastically smaller than
the unobserved maximum based on full data.

Alternatively, if the argument `adjust` to `gev_mle()` is a numeric
scalar then any block maximum for which greater than `adjust` percent of
the underlying raw values are missing is discarded before fitting a GEV
distribution, with no further adjustment made.

### Model diagnostic plots

The plot method for an object returned from `gev_mle` creates visual
model diagnostics like those described in Coles (2001) and implemented
in Heffernan and Stephenson (2018), where the values plotted have been
adjusted, where necessary, for the presence of missing values.

``` r
plot(fit_adjust)
```

<p align="center">

<img src="man/figures/README-plots-1.png" width="60%" style="display: block; margin: auto;" />
</p>

We see that overall the GEV model fits these data well, although the
largest annual maximum sea surge height lies outside the profile-based
95% confidence interval that is relevant to this observation.

## Installation

To install this development version from GitHub:

``` r
remotes::install_github("paulnorthrop/evmiss")
```

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Coles2001" class="csl-entry">

Coles, Stuart G. 2001. *An Introduction to Statistical Modeling of
Extreme Values*. London: Springer.
<https://doi.org/10.1007/978-1-4471-3675-0>.

</div>

<div id="ref-ismev" class="csl-entry">

Heffernan, Janet E., and Alec G. Stephenson. 2018.
*<span class="nocase">i</span>smev: An Introduction to Statistical
Modeling of Extreme Values*.
<https://doi.org/10.32614/CRAN.package.ismev>.

</div>

</div>
