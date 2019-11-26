
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
#> [1,] 9.569273e-04 390.1958 0.3848924 33.85797
#> [2,] 4.067566e-04 387.0640 0.3834602 33.34804
#> [3,] 1.413548e-04 389.4191 0.3749205 32.69723
#> [4,] 4.676487e-04 390.5338 0.3762362 33.01676
#> [5,] 6.423364e-05 390.5068 0.3774113 33.06027
#> [6,] 4.414119e-04 389.8423 0.3789968 32.41419
#> [7,] 1.287751e-03 387.0774 0.3803874 32.67245
#> 
#> [[2]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0002511498 385.9061 0.3899029 33.57525
#> [2,] 0.0012376698 386.5598 0.3883286 33.40353
#> [3,] 0.0012621650 383.8936 0.3944821 33.08603
#> [4,] 0.0003852762 383.3767 0.3925615 32.69749
#> [5,] 0.0019371831 383.5260 0.3937610 33.82304
#> [6,] 0.0015517513 381.7044 0.3995990 33.05927
#> [7,] 0.0014588565 381.7352 0.3988780 34.05518
#> 
#> [[3]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 1.769096e-03 391.2649 0.3757260 33.11727
#> [2,] 5.758551e-04 389.8574 0.3779173 34.59689
#> [3,] 8.458344e-04 391.0381 0.3741445 35.18781
#> [4,] 3.349123e-05 392.1897 0.3727302 32.44165
#> [5,] 6.751491e-04 391.1179 0.3709289 32.27036
#> [6,] 1.079350e-03 390.8434 0.3763794 35.27099
#> [7,] 1.512822e-03 390.0623 0.3726670 32.71365
#> 
#> [[4]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0002970924 383.3299 0.3942570 33.33363
#> [2,] 0.0008734722 384.7989 0.3919450 33.05928
#> [3,] 0.0001001651 384.8165 0.3873583 32.98044
#> [4,] 0.0008350172 385.2169 0.3872892 34.19510
#> [5,] 0.0002610241 386.6675 0.3882188 32.18001
#> [6,] 0.0007794328 387.7779 0.3881649 34.64965
#> [7,] 0.0021415681 385.6544 0.3871649 33.03116
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
#>            Mean        SD  Naive SE Time-series SE
#> L0    7.988e-04 0.0006075 3.037e-06      3.035e-06
#> Linf  3.880e+02 3.1719584 1.586e-02      6.236e-02
#> k     3.819e-01 0.0086671 4.334e-05      1.752e-04
#> sigma 3.354e+01 0.7682464 3.841e-03      7.245e-03
#> 
#> 2. Quantiles for each variable:
#> 
#>            2.5%       25%       50%       75%     97.5%
#> L0    3.233e-05 3.154e-04 6.724e-04   0.00115 2.274e-03
#> Linf  3.818e+02 3.859e+02 3.880e+02 390.16533 3.943e+02
#> k     3.652e-01 3.760e-01 3.817e-01   0.38755 3.994e-01
#> sigma 3.206e+01 3.302e+01 3.353e+01  34.04812 3.509e+01
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
