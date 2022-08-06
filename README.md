
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`tmle3adapt`

[![Travis-CI Build
Status](https://travis-ci.org/tlverse/tmle3adapt.svg?branch=master)](https://travis-ci.org/tlverse/tmle3adapt)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/tlverse/tmle3adapt?branch=master&svg=true)](https://ci.appveyor.com/project/podTockom/tmle3adapt)
[![Coverage
Status](https://img.shields.io/codecov/c/github/tlverse/tmle3adapt/master.svg)](https://codecov.io/github/tlverse/tmle3adapt?branch=master)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

> Targeted Sequential Designs for Learning Optimal Treatment Rules or
> Average Treatment Effects from Optimal Surrogates

**Authors:** [Ivana Malenica](https://github.com/podTockom) and [Mark
van der Laan](https://vanderlaan-lab.org)

------------------------------------------------------------------------

## Description

The `tmle3adapt` is an adapter/extension R package in the `tlverse`
ecosystem, that estimates the mean outcome under the following regimes
in a sequential adaptive trial:

1)  Optimal Individualized Treatment for binary treatment
2)  Average Treatment Effect with/without surrogate outcomes.

------------------------------------------------------------------------

## Installation

You can install the development version of `tmle3adapt` from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with

``` r
devtools::install_github("tlverse/tmle3adapt")
```

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/tlverse/tmle3adapt/issues).

------------------------------------------------------------------------

## Citation

After using the tmle3adapt R package, please cite the following:

``` r
@software{malenica2021tmle3adapt,
      author = {Malenica, Ivana and {van der Laan}, Mark J},
      title = {{tmle3adapt}: Targeted Sequential Designs for Learning Optimal Treatment Rules or Average Treatment Effects from Optimal Surrogates},
      year  = {2021},
      doi = {},
      url = {https://github.com/imalenica/tmle3adapt},
      note = {R package version 0.0.1}
    }
```

## License

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

------------------------------------------------------------------------

## References
