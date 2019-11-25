
<!-- README.md is generated from README.Rmd. Please edit that file -->
BayesGrowth
===========

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/jonathansmart/BayesGrowth.svg?branch=master)](https://travis-ci.org/jonathansmart/BayesGrowth) <!-- badges: end -->

BayesGrowth combines length-at-age modelling for fisheries with MCMC implemented using JAGS and the [rjags](https://cran.r-project.org/web/packages/rjags/index.html) package. Growth modelling using models such as the von Bertalanffy growth model involves three parameters: *L*<sub>∞</sub>, *k* and either *L*<sub>0</sub> or *t*<sub>0</sub>. Two of these parameters: *L*<sub>0</sub> and *L*<sub>∞</sub> have direct biological meaning as the size-at-birth and maximum length, respectively. This package provides the tools to run an MCMC model with these two parameters treated as size-at-birth and maximum length using a JAGS model. This MCMC model is pre-specified and built into wrapper functions.

The user can therefore run an MCMC growth model using knowledge of species length-at-birth and maximum size as priors.

Installation
------------

This package provides a series of wrapper functions to the rjags package which will run JAGS ("Just Another Gibbs Sampler") MCMC models. Therefore, **JAGS must be installed before this package is installed**. To install JAGS, visit: <https://sourceforge.net/projects/mcmc-jags/files/>

You can install the released version of BayesGrowth from [Github](https://github.com/jonathansmart/BayesGrowth) with:

``` r
devtools::install_github("jonathansmart/BayesGrowth")
```

Usage
-----

The main `BayesGrowth` function is `Estimate_MCMC_Growth` which is the wrapper function around an rjags model. It requires a data input that includes columns that can be identified "Age" and "Length", the model needs to be specified (several options are available) and the priors must be specified. Priors include the max size with an error, length-at-birth with an error and upper limits for *k* and *σ*. These latter two parameters have no informative priors and only require sensible upper bounds. Many fish species (including this example) have a size at birth of zero. Therefore, this can value can be used as a prior along with a very small error to indicate high certainty of this prior. The `L0.se` argument cannot be zero, but the model is specified to truncate *L*<sub>0</sub> at zero and keep growth positive.

``` r
library(BayesGrowth)

data("example_data")
## Biological info - lengths in mm
max_size <- 440
max_size_se <- 4.4
birth_size <- 0
birth_size_se <- 0.001 # an se cannot be zero

# Use the function to estimate the JAGS model
results <- Estimate_MCMC_Growth(example_data, Model = "VB" ,
                     iter = 10000,
                     Linf = 440,
                     Linf.se = 4.4,
                     L0 = 0,
                     sigma.max = 100,
                     L0.se = .001,
                     k.max = 1)
```

The function returns the rjags outputs which is an object of class 'mcmc.list'

``` r
head(results)
#> [[1]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0006845955 387.9871 0.3825990 33.51139
#> [2,] 0.0006323733 388.1990 0.3828767 33.35809
#> [3,] 0.0012184070 388.1364 0.3728779 35.09864
#> [4,] 0.0005579542 390.4181 0.3781114 33.96910
#> [5,] 0.0011577562 390.8885 0.3740543 34.00062
#> [6,] 0.0006962200 392.4660 0.3698290 34.10827
#> [7,] 0.0002979405 391.2133 0.3808038 33.44632
#> 
#> [[2]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 6.377257e-04 382.5258 0.3995279 33.71099
#> [2,] 5.080692e-04 382.1476 0.3945358 33.08543
#> [3,] 8.282502e-04 383.5932 0.3993649 33.04936
#> [4,] 4.679651e-05 381.0961 0.4013889 32.23014
#> [5,] 3.517429e-04 381.5760 0.3980532 32.51476
#> [6,] 1.697605e-03 385.1239 0.3866437 33.23379
#> [7,] 1.613641e-03 386.5134 0.3898650 32.88913
#> 
#> [[3]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0001014281 388.7340 0.3792396 34.02489
#> [2,] 0.0008298235 387.5507 0.3829330 34.08608
#> [3,] 0.0014248159 387.0859 0.3857916 33.44174
#> [4,] 0.0006557494 386.7858 0.3848339 33.18208
#> [5,] 0.0025026851 386.3485 0.3903799 33.35219
#> [6,] 0.0007971464 384.8473 0.3864583 33.26862
#> [7,] 0.0001179350 388.1401 0.3834866 33.90521
#> 
#> [[4]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 5.222799e-04 386.7902 0.3894941 33.58785
#> [2,] 8.565135e-04 384.8060 0.3878843 32.50063
#> [3,] 1.722052e-06 385.2598 0.3843889 34.04545
#> [4,] 8.038023e-04 386.4523 0.3870630 32.72811
#> [5,] 6.076832e-05 387.1748 0.3854848 34.72827
#> [6,] 1.502280e-03 387.1375 0.3857406 33.72246
#> [7,] 1.922586e-04 387.3524 0.3835895 34.48434
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
#> L0    7.971e-04 0.000609 3.045e-06      3.045e-06
#> Linf  3.881e+02 3.140774 1.570e-02      6.050e-02
#> k     3.818e-01 0.008597 4.298e-05      1.696e-04
#> sigma 3.354e+01 0.769629 3.848e-03      6.897e-03
#> 
#> 2. Quantiles for each variable:
#> 
#>            2.5%       25%       50%       75%     97.5%
#> L0    2.867e-05 3.144e-04 6.680e-04 1.151e-03 2.266e-03
#> Linf  3.820e+02 3.859e+02 3.880e+02 3.901e+02 3.943e+02
#> k     3.654e-01 3.759e-01 3.817e-01 3.875e-01 3.989e-01
#> sigma 3.208e+01 3.301e+01 3.353e+01 3.405e+01 3.510e+01
plot(results,density = T, smooth = F)
```

<img src="man/figures/README-Diagnostics-1.png" width="100%" />

Additional `BayesGrowth` functions are available that help the user manipulate the returned `Estimate_MCMC_Growth` object. The `Calculate_MCMC_growth_curve` function will provide confidence intervals around the growth curve based on MCMC parameter percentiles. This is essentially a wrapper around the `tidybayes::mean_qi()` function which means it can be passed straight into a ggplot with the `tidybayes::geom_line_ribbon` function.

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

This represents a much improved fit over a standard non-linear estimated model, even if the length-at-birth were fixed at zero. Here the fit is compared using an nls model fit using the [AquaticLifeHistory package](https://github.com/jonathansmart/AquaticLifeHistory).

<img src="man/figures/README-compare_plot-1.png" width="100%" />
