
<!-- README.md is generated from README.Rmd. Please edit that file -->
BayesGrowth
===========

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/jonathansmart/BayesGrowth.svg?branch=master)](https://travis-ci.org/jonathansmart/BayesGrowth) <!-- badges: end -->

BayesGrowth combines length-at-age modelling for fisheries with MCMC implemented using JAGS and the [rjags](https://cran.r-project.org/web/packages/rjags/index.html) package. Growth modelling using models such as the von Bertalanffy growth model involves three parameters: *L*<sub>∞</sub>, *k* and either *L*<sub>0</sub> or *t*<sub>0</sub>. Two of these parameters: *L*<sub>0</sub> and *L*<sub>∞</sub> have direct biological meaning as the size-at-birth and maximum length, respectively. This package provides the tools to run an MCMC model with these two parameters treated as size-at-birth and maximum length using a JAGS model. This MCMC model is pre-specified and built into wrapper functions.

The user can therefore run an MCMC growth model using knowledge of species length-at-birth and maximum size ast priors.

Installation
------------

You can install the released version of BayesGrowth from [Github](https://github.com/jonathansmart/BayesGrowth) with:

``` r
devtools::install_github("jonathansmart/BayesGrowth")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(BayesGrowth)

data("example_data")
## Biological info - lengths in mm
max_size <- 440
max_size_se <- 4.4
birth_size <- 0
birth_size_se <- 0.001 # an se cannot be zero so

# Use the function to estimate the JAGS model
results <- Estimate_MCMC_Growth(example_data, Model = "VB" ,
                     iter = 10000,
                     thin = 1,
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
#> [1,] 0.0005620452 392.0610 0.3723220 33.94318
#> [2,] 0.0001369667 389.9172 0.3719000 34.37843
#> [3,] 0.0024029083 391.4290 0.3718815 33.94910
#> [4,] 0.0007986814 390.5713 0.3740651 34.59104
#> [5,] 0.0003338447 392.0651 0.3714542 34.41875
#> [6,] 0.0014637191 391.7149 0.3684055 34.31933
#> [7,] 0.0003701148 392.6997 0.3682633 32.82473
#> 
#> [[2]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 1.437831e-03 391.1706 0.3686934 33.68340
#> [2,] 5.876792e-04 391.6915 0.3714747 34.09151
#> [3,] 7.586228e-04 392.5083 0.3685232 34.56320
#> [4,] 1.450904e-05 389.9673 0.3697185 32.49913
#> [5,] 9.742589e-04 391.8381 0.3727185 34.91729
#> [6,] 6.155311e-04 392.1416 0.3698109 34.80104
#> [7,] 1.829126e-03 394.7320 0.3689887 34.33722
#> 
#> [[3]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 3.278961e-04 384.0970 0.3948504 33.13292
#> [2,] 1.368021e-03 383.1382 0.3948640 32.31931
#> [3,] 1.792267e-04 382.7616 0.3968951 32.05192
#> [4,] 1.236768e-03 382.3860 0.3968147 32.80638
#> [5,] 3.561401e-05 382.7382 0.3933473 32.01858
#> [6,] 1.241996e-03 383.8103 0.3935411 32.53760
#> [7,] 1.299801e-05 380.5846 0.3941646 32.11056
#> 
#> [[4]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0008450929 385.2614 0.3947052 33.65721
#> [2,] 0.0007748658 382.6058 0.3948355 33.64765
#> [3,] 0.0007460382 384.5935 0.3940808 33.63621
#> [4,] 0.0015213388 383.0320 0.3921592 32.67531
#> [5,] 0.0019646179 385.3097 0.3902396 33.14265
#> [6,] 0.0001802971 384.8634 0.3911948 33.54316
#> [7,] 0.0002729824 387.5196 0.3891440 32.95159
#> 
#> attr(,"class")
#> [1] "mcmc.list"
```

Therefore, all of the diagnostics of the rjags library can be used.

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
#> L0    7.965e-04 0.000604 3.020e-06      3.057e-06
#> Linf  3.880e+02 3.141438 1.571e-02      6.146e-02
#> k     3.818e-01 0.008558 4.279e-05      1.701e-04
#> sigma 3.354e+01 0.758228 3.791e-03      6.948e-03
#> 
#> 2. Quantiles for each variable:
#> 
#>            2.5%       25%       50%       75%     97.5%
#> L0    3.113e-05 3.145e-04 6.735e-04 1.148e-03   0.00223
#> Linf  3.820e+02 3.860e+02 3.880e+02 3.901e+02 394.39122
#> k     3.652e-01 3.761e-01 3.817e-01 3.875e-01   0.39887
#> sigma 3.210e+01 3.302e+01 3.353e+01 3.405e+01  35.05989
plot(results,density = T, smooth = F)
```

<img src="man/figures/README-Diagnostics-1.png" width="100%" />

Additional `BayesGrowth` functions are available that help the user manipulate the objects. The `Calculate_MCMC_growth_curve` function will provide confidence intervals around the growth curve based on MCMC percentiles. This is essentially a wrapper around the `tidybayes::mean_qi()` which means it can be passed straight into a ggplot with the `tidybayes::geom_line_ribbon`.

``` r

growth_curve <- Calculate_MCMC_growth_curve(results, Model = "VB",
                                            max.age = max(example_data$Age), probs = .95)
library(tidybayes)
library(ggplot2)

ggplot(growth_curve, aes(Age, LAA))+
  geom_point(data = example_data, aes(Age, Length), alpha = .3)+
  geom_lineribbon(size = 1.5, fill = "steelblue") +
  labs(y = "Total Length (mm)", x = "Age (yrs)")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,13,1))+
  theme_bw()+
  theme(text = element_text(size = 14))
```

<img src="man/figures/README-plot-1.png" width="100%" />

This represents a much improved fit over a standard non-linear estimated model, even if the length-at-birth were fixed at zero. Here the fit is compared using an nls model fit using the [AquaticLifeHistory package](https://github.com/jonathansmart/AquaticLifeHistory).

<img src="man/figures/README-compare plot-1.png" width="100%" />
