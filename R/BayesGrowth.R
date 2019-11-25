#' BayesGrowth
#'
#' @name BayesGrowth
#' @author Jonathan Smart
#' @description A package to estimate length-at-age using MCMC using rjags models
#' @docType package
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


