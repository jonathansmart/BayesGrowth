#' BayesGrowth
#'
#' @description A package to estimate fish length-at-age models using MCMC analysis with rstan models.
#'
#' @docType package
#' @name BayesGrowth
#' @aliases BayesGrowth
#' @useDynLib BayesGrowth, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @import RcppParallel
#' @import bayesplot
#' @import ggplot2
#' @import AquaticLifeHistory
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
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
#' @format An rstan model with the class stan.fit
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
