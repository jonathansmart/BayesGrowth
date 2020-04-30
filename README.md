
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesGrowth

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/jonathansmart/BayesGrowth.svg?branch=master)](https://travis-ci.org/jonathansmart/BayesGrowth)
[![DOI](https://zenodo.org/badge/223878432.svg)](https://zenodo.org/badge/latestdoi/223878432)
<!-- badges: end -->

BayesGrowth combines length-at-age modelling for fisheries with MCMC
implemented using JAGS and the
[rjags](https://cran.r-project.org/web/packages/rjags/index.html)
package. Growth modelling using models such as the von Bertalanffy
growth model involves three parameters: \(L_{\infty}\), \(k\) and either
\(L_{0}\) or \(t_{0}\). Two of these parameters: \(L_{0}\) and
\(L_{\infty}\) have direct biological meaning as the size-at-birth and
maximum length, respectively. This package provides the tools to run an
MCMC model with these two parameters treated as size-at-birth and
maximum length using a JAGS model. This MCMC model is pre-specified and
built into wrapper functions.

The user can therefore run an MCMC growth model using knowledge of
species length-at-birth and maximum size as priors.

## Installation

This package provides a series of wrapper functions to the rjags package
which will run JAGS (“Just Another Gibbs Sampler”) MCMC models.
Therefore, **JAGS must be installed before this package is installed**.
To install JAGS, visit:
<https://sourceforge.net/projects/mcmc-jags/files/>

You can install the released version of BayesGrowth from
[Github](https://github.com/jonathansmart/BayesGrowth) with:

``` r
devtools::install_github("jonathansmart/BayesGrowth")
```

## Usage

The main `BayesGrowth` function is `Estimate_MCMC_Growth` which is the
wrapper function around an rjags model. It requires a data input that
includes columns that can be identified “Age” and “Length”, the model
needs to be specified (several options are available) and the priors
must be specified. Priors include the max size with an error,
length-at-birth with an error and upper limits for \(k\) and \(\sigma\).
These latter two parameters have no informative priors and only require
sensible upper bounds. Many fish species (including this example) have a
size at birth of zero. Therefore, this can value can be used as a prior
along with a very small error to indicate high certainty of this prior.
The `L0.se` argument cannot be zero, but the model is specified to
truncate \(L_{0}\) at zero and keep growth positive.

``` r
library(BayesGrowth)

data("example_data")
## Biological info - lengths in mm
max_size <- 440
max_size_se <- 5
birth_size <- 0
birth_size_se <- 0.001 # an se cannot be zero

# Use the function to estimate the JAGS model
results <- Estimate_MCMC_Growth(example_data, Model = "VB" ,
                     iter = 10000,
                     Linf = max_size,
                     Linf.se = max_size_se,
                     L0 = birth_size,
                     sigma.max = 100,
                     L0.se = birth_size_se,
                     k.max = 1)
```

The function returns the rjags outputs which is an object of class
‘mcmc.list’

``` r
head(results)
#> [[1]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 4.488882e-05 325.6767 0.5960316 26.59236
#> [2,] 7.517356e-05 326.6002 0.5986476 25.38877
#> [3,] 7.385076e-05 324.2891 0.5983142 25.59835
#> [4,] 4.121038e-04 325.8852 0.5922294 25.69103
#> [5,] 1.641973e-05 326.0010 0.5934298 27.12353
#> [6,] 4.553976e-04 327.8349 0.5934005 25.19848
#> [7,] 1.578843e-03 326.0242 0.6058683 26.54532
#> 
#> [[2]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0017147504 326.8731 0.5991217 26.32451
#> [2,] 0.0002984886 325.8409 0.5981972 26.96292
#> [3,] 0.0011963286 325.8721 0.6098430 25.53892
#> [4,] 0.0009314804 322.8726 0.6100467 25.00895
#> [5,] 0.0018254837 324.4057 0.6064083 25.90006
#> [6,] 0.0007660022 325.2780 0.5982529 25.27672
#> [7,] 0.0009750646 325.8695 0.6000366 25.25688
#> 
#> [[3]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 1.662549e-04 324.4075 0.6065113 25.17975
#> [2,] 8.069489e-04 325.5753 0.6096032 25.85848
#> [3,] 1.256401e-03 323.6330 0.6161591 23.89699
#> [4,] 6.643710e-04 323.1181 0.6187880 24.88484
#> [5,] 3.401003e-05 321.0391 0.6208219 24.71695
#> [6,] 1.345030e-03 322.6422 0.6236229 24.83958
#> [7,] 1.342958e-03 323.6096 0.6241476 25.59646
#> 
#> [[4]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0003766009 318.1324 0.6740486 24.93456
#> [2,] 0.0014715359 314.3714 0.6703935 24.87597
#> [3,] 0.0020864968 316.0134 0.6830210 24.68304
#> [4,] 0.0012641278 313.3945 0.6883822 24.58046
#> [5,] 0.0016279314 316.1042 0.6809096 23.23320
#> [6,] 0.0003830199 315.6487 0.6773715 25.33736
#> [7,] 0.0002557087 317.1293 0.6666446 23.74806
#> 
#> attr(,"class")
#> [1] "mcmc.list"
```

Therefore, all of the diagnostics from the rjags library can be used.

``` r
summary(results)
#> 
#> Iterations = 1001:11000
#> Thinning interval = 1 
#> Number of chains = 4 
#> Sample size per chain = 10000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>            Mean       SD  Naive SE Time-series SE
#> L0    8.003e-04 0.000603 3.015e-06      3.004e-06
#> Linf  3.180e+02 4.146394 2.073e-02      1.424e-01
#> k     6.604e-01 0.034451 1.723e-04      1.171e-03
#> sigma 2.431e+01 0.868429 4.342e-03      2.221e-02
#> 
#> 2. Quantiles for each variable:
#> 
#>            2.5%       25%       50%       75%     97.5%
#> L0    3.341e-05 3.211e-04 6.779e-04 1.153e-03 2.234e-03
#> Linf  3.108e+02 3.151e+02 3.177e+02 3.204e+02 3.272e+02
#> k     5.901e-01 6.384e-01 6.605e-01 6.836e-01 7.264e-01
#> sigma 2.274e+01 2.371e+01 2.427e+01 2.486e+01 2.614e+01
plot(results,density = T, smooth = F)
```

<img src="man/figures/README-Diagnostics-1.png" width="100%" />

Additional `BayesGrowth` functions are available that help the user
manipulate the returned `Estimate_MCMC_Growth` object. The
`Calculate_MCMC_growth_curve` function will provide confidence intervals
around the growth curve based on MCMC parameter percentiles. This is
essentially a wrapper around the `tidybayes::mean_qi()` function which
means it can be passed straight into a ggplot with the
`tidybayes::geom_line_ribbon` function.

``` r

# Return a growth curve with 50th and 95th percentiles
growth_curve <- Calculate_MCMC_growth_curve(results, Model = "VB",
                                            max.age = max(example_data$Age), probs = c(.5,.95))
library(tidybayes)
library(ggplot2)

ggplot(growth_curve, aes(Age, LAA))+
  geom_point(data = example_data, aes(Age, Length), alpha = .3)+
  geom_lineribbon( size = .8) +
  labs(y = "Total Length (mm)", x = "Age (yrs)")+
  scale_fill_brewer(palette="BuPu", direction=-1)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,13,1))+
  theme_bw()+
  theme(text = element_text(size = 14))
```

<img src="man/figures/README-plot-1.png" width="100%" />

This represents a much improved fit over a standard non-linear estimated
model, even if the length-at-birth were fixed at zero. Here the fit is
compared using an nls model fit using the [AquaticLifeHistory
package](https://github.com/jonathansmart/AquaticLifeHistory).

<img src="man/figures/README-compare_plot-1.png" width="100%" />
