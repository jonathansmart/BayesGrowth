
<!-- README.md is generated from README.Rmd. Please edit that file -->
BayesGrowth
===========

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/jonathansmart/BayesGrowth.svg?branch=master)](https://travis-ci.org/jonathansmart/BayesGrowth) <!-- badges: end -->

BayesGrowth combines length-at-age modelling for fisheries with MCMC implemented using JAGS and the [rjags](https://cran.r-project.org/web/packages/rjags/index.html) package. Growth modelling using models such as the von Bertalanffy growth model involves three parameters: *L*<sub>∞</sub>, *k* and either *L*<sub>0</sub> or *t*<sub>0</sub>. Two of these parameters: *L*<sub>0</sub> and *L*<sub>∞</sub> have direct biological meaning as the size-at-birth and maximum length, respectively. This package provides the tools to run an MCMC model with these two parameters treated as size-at-birth and maximum length using a JAGS model. This MCMC model is pre-specified and built into wrapper functions.

The user can therefore run an MCMC growth model using knowledge of species length-at-birth and maximum size ast priors.

Installation
------------

This package provides a series of wrapper functions to the rjags package which will run JAGS ("Just Another Gibbs Sampler") MCMC models. Therefore, *JAGS must be installed before this package is installed*. To install JAGS, visit: <https://sourceforge.net/projects/mcmc-jags/files/>

You can install the released version of BayesGrowth from [Github](https://github.com/jonathansmart/BayesGrowth) with:

``` r
devtools::install_github("jonathansmart/BayesGrowth")
```

Usage
-----

The main `BayesGrowth` function is `Estimate_MCMC_Growth` which is the wrapper function around an rjags model. It requires a data input that includes Age and Length, a model needs to be specified (several options are available) and the priors must be specified. Priors include the max size with an error, length-at-birth with an error, an upper value for *k* and *σ*. These latter two parameters have no informative priors so just need senscible upper bounds. Many fish species (including this example) have a size at birth of zero. Therefore, this can be used along with a very small error to indicate that this is a certain value. The `L0.se` argument cannot be zero, although the model is specified to truncate it to zero and keep growth positive.

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
#> [1,] 0.0007140421 389.2083 0.3812362 33.88953
#> [2,] 0.0009684128 388.0148 0.3774248 33.74464
#> [3,] 0.0000225300 389.3610 0.3833927 33.02893
#> [4,] 0.0003189595 385.8569 0.3813557 33.42119
#> [5,] 0.0001301606 387.0420 0.3821805 33.38645
#> [6,] 0.0008697618 387.2100 0.3810786 34.72908
#> [7,] 0.0001933690 389.9525 0.3807542 32.30058
#> 
#> [[2]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0009373289 388.8947 0.3742611 33.74613
#> [2,] 0.0006973102 390.0837 0.3778313 33.04571
#> [3,] 0.0019347654 390.5374 0.3731896 32.99147
#> [4,] 0.0009321508 392.9994 0.3666824 33.76493
#> [5,] 0.0001182634 393.4887 0.3645631 34.17591
#> [6,] 0.0013784802 395.0423 0.3648289 33.78332
#> [7,] 0.0019840112 393.2772 0.3741709 33.88302
#> 
#> [[3]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 1.524884e-03 382.0567 0.3973708 32.04994
#> [2,] 3.744324e-04 383.6138 0.3974500 32.37718
#> [3,] 1.855561e-03 382.7780 0.3979980 32.60266
#> [4,] 6.220523e-05 382.1310 0.3978527 32.35666
#> [5,] 4.363249e-04 383.0975 0.3899973 33.08628
#> [6,] 5.341089e-04 385.7289 0.3898472 33.88497
#> [7,] 7.557750e-04 385.3176 0.3869741 33.16633
#> 
#> [[4]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 1.589315e-03 387.0945 0.3808870 33.14831
#> [2,] 1.094638e-03 386.6254 0.3860141 33.50220
#> [3,] 4.755118e-04 388.3166 0.3865569 32.66871
#> [4,] 5.591126e-05 386.3635 0.3879419 33.04041
#> [5,] 3.037665e-04 386.4624 0.3811012 33.92760
#> [6,] 1.787088e-04 387.6920 0.3800311 31.29473
#> [7,] 8.820132e-05 388.4971 0.3803537 31.60521
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
#> L0    7.980e-04 0.0006064 3.032e-06      3.017e-06
#> Linf  3.881e+02 3.1672541 1.584e-02      6.144e-02
#> k     3.816e-01 0.0086143 4.307e-05      1.680e-04
#> sigma 3.355e+01 0.7696899 3.848e-03      7.130e-03
#> 
#> 2. Quantiles for each variable:
#> 
#>            2.5%       25%       50%       75%     97.5%
#> L0    3.149e-05 3.145e-04 6.732e-04 1.151e-03 2.254e-03
#> Linf  3.820e+02 3.860e+02 3.881e+02 3.902e+02 3.945e+02
#> k     3.648e-01 3.758e-01 3.816e-01 3.873e-01 3.989e-01
#> sigma 3.210e+01 3.303e+01 3.354e+01 3.406e+01 3.512e+01
plot(results,density = T, smooth = F)
```

<img src="man/figures/README-Diagnostics-1.png" width="100%" />

Additional `BayesGrowth` functions are available that help the user manipulate the objects. The `Calculate_MCMC_growth_curve` function will provide confidence intervals around the growth curve based on MCMC percentiles. This is essentially a wrapper around the `tidybayes::mean_qi()` function which means it can be passed straight into a ggplot with the `tidybayes::geom_line_ribbon` function.

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

<img src="man/figures/README-compare_plot-1.png" width="100%" />
