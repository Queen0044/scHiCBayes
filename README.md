
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scHiCBayes

<!-- badges: start -->

<!-- badges: end -->

**scHiCBayes** (Xie, Han, Jin, and Lin, 2021) is a Bayesian hierarchy
model that goes beyond data quality improvement by also identifying
observed zeros that are in fact structural zeros. scHiCBayes takes
spatial dependencies of scHi-C 2D data structure into account while also
borrowing information from similar single cells and bulk data, when such
are available.

## Installation

The scHiCBayes package has the following R-package dependencies: Rcpp,
RcppArmadillo, parallel, Rtsne, ggplot2, ggpubr, and mclust. The
dependent packages will be automatically installed along with
scHiCBayes. You can use the following commands to install scHiCBayes
from GitHub.

``` r
# Install and load "devtools" package. 
install.packages("devtools")
library("devtools")

# Install "scHiCBayes" package from github.
install_github("Queen0044/scHiCBayes")
```

If you are Windows user, please install **Rtools40**
(<https://cran.r-project.org/bin/windows/Rtools/>) first, and restart R
to install scHiCBayes package. Note that Rtools40 is for R 4.0.0+ so
that you might have to update your R version.

If you have OneDrive backing-up “C:\\User\\Your\_user\_name\\Documents”,
the installation may fail. You can download the zip file from Github and
install scHiCBayes manually.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scHiCBayes)
#> 
#> Attaching package: 'scHiCBayes'
#> The following object is masked from 'package:stats':
#> 
#>     heatmap
#data("K562_T1_7k")
#data("K562_bulk")
#single=K562_T1_7k
#T1_7k_res=MCMCImpute(niter=100000,burnin=5000,single=K562_T1_7k,bulk=K562_bulk,
#startval=c(100,100,10,8,10,0.1,900,0.2,0,replicate(dim(single)[2],8)),n=61,mc.cores = 1,cutoff=0.5)
```

For more information of functions, please read the vignettes.
