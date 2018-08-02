
<!-- README.md is generated from README.Rmd. Please edit that file -->
malariaModelFit
===============

[![Travis build status](https://travis-ci.org/mrc-ide/malariaModelFit.svg?branch=master)](https://travis-ci.org/mrc-ide/malariaModelFit)

Rcpp package containing code for fitting the malaria model by MCMC.

Installation
------------

*malariaModelFit* relies on the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) package, which requires the following OS-specific steps:

-   Windows
    -   Download and install the appropriate version of [Rtools](https://cran.rstudio.com/bin/windows/Rtools/) for your version of R. On installation, ensure you check the box to arrange your system PATH as recommended by Rtools
-   Mac OS X
    -   Download and install [XCode](http://itunes.apple.com/us/app/xcode/id497799835?mt=12)
    -   Within XCode go to Preferences : Downloads and install the Command Line Tools
-   Linux (Debian/Ubuntu)
    -   Install the core software development utilities required for R package development as well as LaTeX by executing

            sudo apt-get install r-base-dev texlive-full

Next, in R, ensure that you have the [devtools](https://www.rstudio.com/products/rpackages/devtools/) package installed by running

``` r
install.packages("devtools", repos='http://cran.us.r-project.org')
```

Then install the *malariaModelFit* package directly from GitHub by running

``` r
devtools::install_github("mrc-ide/malariaModelFit")
```

Finally, we need to load the package:

``` r
library(malariaModelFit)
```
