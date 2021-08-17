# netcom <img src="man/figures/netcom_sticker.png" align="right" width="120" />

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/netcom)](https://cran.r-project.org/package=netcom)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/netcom)](https://cran.r-project.org/package=netcom)

`netcom` is an R package to infer system functioning by empirically comparing networks to each other. There are many uses of this general approach to understanding network data. The vignette covers several common kinds of inference. Once installed, this can be opened by running:

```R
vignette("tutorial", package = "netcom")
```

You can learn more about comparative inference in the following:

- Langendorf, R. E. & Burgess, M. G. (2020) Empirically classifying network mechanisms. arXiv preprint arXiv:2012.15863.

- Langendorf, R. E. (2020). Using Structural Comparisons to Measure the Behavior of Complex Systems. In Evaluating Climate Change Impacts (pp. 123-137). Chapman and Hall/CRC.

- Langendorf, R. E. & Goldberg, D. S. (2019) Aligning statistical dynamics captures biological network functioning. arXiv preprint arXiv:1912.12551.

You can install the **netcom** package two main ways:

1. A release version of the package can be installed from CRAN (the Comprehensive R Archive Network): https://cran.r-project.org/package=netcom.

```R
install.packages("netcom").
```

2. Alternatively, the (usually) more recent development version can be installed from GitHub: https://github.com/langendorfr/netcom. This can be accomplished with the **devtools** package. We recommend new users install the other version, from CRAN, which is has less functioning but has been more reliably tested.

```R
install.packages("devtools")
devtools::install_github("langendorfr/netcom")
```
