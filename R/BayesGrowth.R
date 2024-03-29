#' BayesGrowth
#'
#' @description A package to estimate fish length-at-age models using MCMC analysis with 'rstan' models.
#'
#' @docType package
#' @name BayesGrowth
#' @description
#' Estimate fish length-at-age models using MCMC analysis with
#' 'rstan' models. This package allows a multimodel approach to growth fitting
#' to be applied to length-at-age data and is supported by further analyses to
#' determine model selection and result presentation. The core methods of this
#' package are presented in Smart and Grammer (2021) "Modernising fish and
#' shark growth curves with Bayesian length-at-age models". PLOS ONE 16(2):
#'   e0246734.
#' @aliases BayesGrowth
#' @useDynLib BayesGrowth, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @import bayesplot
#' @import ggplot2
#' @import AquaticLifeHistory
#' @importFrom RcppParallel CxxFlags RcppParallelLibs
#' @docType package
#' @references To cite the BayesGrowth package in publications, type citation('BayesGrowth').
#' The Stan software should also be referenced:
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. <https://mc-stan.org>
#'
NULL

globalVariables(c( "Age", "LAA", "Parameter", "Penalised deviance", "Value", "quantile", "sd", "sigma"))

#' Example length-at-age data
#'
#' A dataset used as examples in the vignettes and for users to test code with
#'
#' \itemize{
#'   \item Age. Number of growth bands determined from vertebral analysis
#'   \item Length. Total Length in mm determined via back-calculation
#' }
#'
#' @docType data
#' @keywords datasets
#' @name example_data
#' @usage data(example_data)
#' @format A data frame with 509 rows and 2 variables
NULL


#' Example MCMC model outputs
#'
#' The results of an MCMC model used as examples in the vignettes and for users to test code with
#'
#' @docType data
#' @keywords datasets
#' @name MCMC_example_results
#' @usage data(MCMC_example_results)
#' @format An 'rstan' model with the class stan.fit
NULL

#' Example MCMC model comparison outputs outputs
#'
#' The results of an MCMC model comparison used as examples in the vignettes and for users to test code with
#'
#' @docType data
#' @keywords datasets
#' @name Looic_example_results
#' @usage data(Looic_example_results)
#' @format An dataframe with results of LooIC
NULL
