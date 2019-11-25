
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
#> [1,] 0.0006652715 390.9521 0.3717089 34.23016
#> [2,] 0.0002653097 391.8202 0.3741431 34.12191
#> [3,] 0.0008244993 391.3167 0.3752415 33.85880
#> [4,] 0.0002247448 390.2134 0.3762566 33.34104
#> [5,] 0.0009638235 390.3530 0.3751728 32.35000
#> [6,] 0.0004318568 389.5652 0.3738423 32.36783
#> [7,] 0.0000784418 389.7947 0.3724705 32.89823
#> 
#> [[2]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0003172034 383.8720 0.3879421 33.79979
#> [2,] 0.0003126417 389.0027 0.3790549 33.81706
#> [3,] 0.0003552461 387.1678 0.3852045 32.89205
#> [4,] 0.0010176001 386.8890 0.3853119 33.15583
#> [5,] 0.0023008715 387.9941 0.3819780 32.42553
#> [6,] 0.0005101382 388.5263 0.3830371 33.42809
#> [7,] 0.0001943342 387.2390 0.3861711 34.75159
#> 
#> [[3]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0003541779 388.3304 0.3739821 33.64270
#> [2,] 0.0017186522 389.0593 0.3693600 33.48490
#> [3,] 0.0014169990 391.0227 0.3716809 33.82420
#> [4,] 0.0005503210 391.0124 0.3766435 33.44774
#> [5,] 0.0016794807 388.7599 0.3800068 33.19777
#> [6,] 0.0007007575 389.1321 0.3823088 33.41521
#> [7,] 0.0007860959 386.6312 0.3818066 33.26740
#> 
#> [[4]]
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1001 
#> End = 1007 
#> Thinning interval = 1 
#>                L0     Linf         k    sigma
#> [1,] 0.0003068650 383.7703 0.3970474 33.50448
#> [2,] 0.0006980269 384.1806 0.3917652 32.61900
#> [3,] 0.0006432910 382.1401 0.3915466 32.32490
#> [4,] 0.0004167164 386.0858 0.3919686 34.80228
#> [5,] 0.0002497434 383.3461 0.3948591 34.43699
#> [6,] 0.0008830061 383.6027 0.3943774 34.14556
#> [7,] 0.0007979783 384.3037 0.3958615 34.08158
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
#> L0    7.979e-04 0.0005993 2.997e-06      2.985e-06
#> Linf  3.880e+02 3.1929815 1.596e-02      6.364e-02
#> k     3.820e-01 0.0087000 4.350e-05      1.765e-04
#> sigma 3.354e+01 0.7759623 3.880e-03      7.147e-03
#> 
#> 2. Quantiles for each variable:
#> 
#>            2.5%       25%       50%       75%     97.5%
#> L0    3.056e-05 3.194e-04 6.814e-04 1.149e-03 2.231e-03
#> Linf  3.818e+02 3.858e+02 3.880e+02 3.901e+02 3.943e+02
#> k     3.655e-01 3.760e-01 3.818e-01 3.878e-01 3.994e-01
#> sigma 3.206e+01 3.301e+01 3.353e+01 3.406e+01 3.510e+01
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

<img src="man/figures/README-compare plot-1.png" width="100%" />
