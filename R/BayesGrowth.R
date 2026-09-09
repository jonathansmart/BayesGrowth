#' BayesGrowth: Bayesian Fish Length-at-Age Models
#'
#' Estimate fish length-at-age models using MCMC analysis with
#' 'rstan' models. This package allows a multimodel approach to growth fitting
#' to be applied to length-at-age data and is supported by further analyses to
#' determine model selection and result presentation. The core methods of this
#' package are presented in Smart and Grammer (2021) "Modernising fish and
#' shark growth curves with Bayesian length-at-age models". PLOS ONE 16(2):
#'   e0246734.
#'
#' @useDynLib BayesGrowth, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @import bayesplot
#' @import ggplot2
#' @import AquaticLifeHistory
#' @importFrom RcppParallel CxxFlags RcppParallelLibs
#' @references To cite the BayesGrowth package in publications, type citation('BayesGrowth').
#' The Stan software should also be referenced:
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. <https://mc-stan.org>
#'
"_PACKAGE"

globalVariables(c( "Age", "LAA", "Parameter", "Penalised deviance", "Value", "quantile", "sd", "sigma"))

#' Example length-at-age data
#'
#' A dataset used as examples in the vignettes and for users to test code with
#'
#' \itemize{
#'   \item Age. Number of growth bands determined from otolith ageing
#'   \item Length. Total Length in mm determined
#' }
#'
#' @docType data
#' @keywords datasets
#' @name example_data
#' @usage data(example_data)
#' @format A data frame with 509 rows and 2 variables
NULL

#' Example back-calculated length-at-age data
#'
#' A dataset used as examples in the vignettes and for users to test code with
#'
#' \itemize{
#'   \item Age. Number of growth bands determined from vertebral analysis
#'   \item Length. Total Length in mm determined via back-calculation
#'   \item Sex. F for female and M for male
#'   \item Tag. ID of each individual for back-calculated length-at-age
#' }
#'
#' @docType data
#' @keywords datasets
#' @name example_BC_data
#' @usage data(example_BC_data)
#' @format A data frame with 294 rows and 4 variables
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
