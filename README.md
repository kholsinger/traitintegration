`traitintegration` is a package that facilitates analyses of trait
integration. It provides functions to use the output of an analysis by
`brms` to

- Visualize the pairwise correlations among traits in the analysis
- Calculate, summarize, and visualize the posterior distribution of
  three statistics used to summarize the amount of trait integration
	- The maximum eigenvalue of the correlation matrix.
	- The variance of eigenvalues of the correlation matrix.
	- entR, the correlation coefficient of a correlation matrix that
      has a constant value, r, on the off-diagonal and the same
      entropy as the observed correlation matrix (assuming that the
      traits follow a multivariate normal distribution)
      
## Installation

To install `traitintegration` you will need the `devtools` package. If
you don't already have it, simply `install.packages("devtools")`. Once
you've done that, simply
`devtools::install_github("kholsinger/traitintegration",
build_vignettes = TRUE)`. It will churn for a while, and you may be
asked to install some new packages, upgrade some existing packages, or
both. 

Assuming that the installations succeeds, simply
`library(traitintegration)` and `vignette("Introduction", package =
"traitintegration")` to get an overview of how to use
`traitintegration`. You can also look at the vignette here simply by
clicking this link: [Using `brms` and
`traitintegration`](https://kholsinger.github.io/traitintegration/vignettes/Introduction.html) 
