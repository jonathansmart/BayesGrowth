#' Estimate_MCMC_Growth
#' @description A wrapper function that creates a JAGS model using the rjags package.
#'     The data and priors provided are combined into an rjags model that estimates a length-at-age
#'     model with a normal distribution. Three different growth models can be used: a von Bertalanffy
#'     model, Gompertz model or a logistic model. The prior on Linf and L0 are normally distributed
#'     and determined through the user providing a mean and se for eacd parameter. The growth completion
#'     parameter or any model (k) has a uniform prior which only requires an upper bound with the
#'     lower bound set at zero. Sigma is the residual variance of the data around the model and
#'     is set up in the same manner as 'k'. The growth estimates in the model are truncated to
#'     remain above zero so negative growth cannot occur.
#' @param data A data.frame that contains columns named 'Age' and "Length'. The function can
#'     detect columns with similar names. If age and length columns cannot be determined then
#'     an error will occur. The dataset can have additional columns which will be ignored by
#'     the function
#' @param Model Which growth model should be run? Must be one of "VB", "Gom" or "Log" for von
#'     Bertalanffy, Gompertz or Logistic models, respectively
#' @param Linf The prior for asymptotic length. MUst be in the same unit (i.e. cm or mm) as the data.
#'     This should be based off of maximum size for the species.
#' @param Linf.se The prior for normally distributed standard error around
#'     asymptotic length. MUst be in the same unit (i.e. cm or mm) as the data. Cannot be zero.
#' @param L0 The prior for length-at-birth. MUst be in the same unit (i.e. cm or mm) as the data.
#'     This should be based off of minimum size for the species.
#' @param L0.se The prior for normally distributed standard error around
#'     length-at-birth. MUst be in the same unit (i.e. cm or mm) as the data. Cannot be zero.
#' @param k.max The maximum value to consider for the growth completion parameter 'k'. In
#'     the Gompertz and Logistic models, this parameter is often notated as 'g' instead of 'k'.
#' @param sigma.max The maximum value to consider for sigma. This is the variance around the
#'     length-at-age residuals.
#' @param iter How many MCMC iterations should be run? Default is 10000 but fewer can be useful to
#'     avoid longer run times when testing code or data
#' @param BurnIn The number of iterations at the begining of each chain to discard ('Burn in') to
#'     avoid biased values from starting values that do not resemble the target distribution.
#'     Default is 1000.
#' @param n.chains Number of MCMC chains to be run. Default is 4.
#' @param thin The thinning of the MCMC simulations. Deafault is 1 which means no thinning occurs.
#'     Thinning is generally only necessary for complicated models as it increases run time.
#' @import dplyr rjags
#'
#' @return An object of class 'mcmc.list'. This is the rjags object that returns a list of
#'     MCMC outputs from the rjags model. Each list element is the results of each thinned
#'     simulation from each chain.
#'
#' @export

Estimate_MCMC_Growth <- function(data,  Model = NULL, Linf = NULL, Linf.se = NULL,
                                 L0 = NULL, L0.se = NULL, k.max = NULL, sigma.max = NULL,
                                 iter = 10000, BurnIn = 1000,
                                 n.chains = 4, thin = 1){

  if(any(is.null(c(Linf, Linf.se, L0, L0.se, k.max, sigma.max)))) stop("At least one parameter or its error are not correctly specified")
  if(length(Model) != 1) stop("Only one growth model can be used in each function call")
  if(is.null(Model))stop("Growth model has not been specified")


  age_col <- grep("age", substr(tolower(names(data)),1,3))
  if(length(age_col) <1) stop("Age column heading could not be distinguished ")
  if(length(age_col) >1) stop("Multiple age columns detected. Remove unecessary variables or rename desired column to 'Age' ")

  len_col <- grep("len|tl|lt|siz", substr(tolower(names(data)),1,3))
  if(length(len_col) <1) stop("Length column heading could not be distinguished ")
  if(length(len_col) >1) stop("Multiple length columns detected. Remove unecessary variables or rename desired column to 'Length' ")


  if(Linf.se == 0 | L0.se == 0) stop("L0 and Linf standard error priors cannot be zero")
  if(any(is.na(data))) stop("data contains NA's")


  Age <- data[,age_col]
  Length <- data[,len_col]


  if(Model == "VB"){
    Growth_model <- paste0(
      "model{
      # Likelihood model for Length[i]
      for(i in 1:length(Length)) {
      Length[i] ~ dnorm(Lt[i], sigma^(-2))
      Lt[i] <- Linf-(Linf-L0)*exp(-k*Age[i])
      }

      # Prior models for a, b, s
      Linf ~ dnorm(",Linf,",",Linf.se,"^(-2))
      k ~ dunif(0,",k.max,")#,
      # L0_dist,
      L0 ~ dnorm(",L0,",",L0.se,"^(-2)) T(0,",max(Length),")
      sigma ~ dunif(0,",sigma.max,")
  }")
  } else if(Model == "Gom"){
    Growth_model <- paste0(
      "model{
      # Likelihood model for Length[i]
      for(i in 1:length(Length)) {
      Length[i] ~ dnorm(Lt[i], sigma^(-2))
      Lt[i] <-L0*exp(log(Linf/L0)*(1-exp(-k*Age[i])))
      }

      # Prior models for a, b, s
      Linf ~ dnorm(",Linf,",",Linf.se,"^(-2))
      k ~ dunif(0.0001,",k.max,")
      L0 ~ dnorm(",L0,",",L0.se,"^(-2)) T(0,",max(Length),")
      sigma ~ dunif(0,",sigma.max,")
}")
  } else if(Model == "Log"){
    Growth_model <- paste0(
      "model{
      # Likelihood model for Length[i]
      for(i in 1:length(Length)) {
      Length[i] ~ dnorm(Lt[i], sigma^(-2))
      Lt[i] <-(Linf*L0*exp(k*Age[i]))/(Linf+L0*(exp(k*Age[i])-1))
      }

      # Prior models for a, b, s
      Linf ~ dnorm(",Linf,",",Linf.se,"^(-2))
      k ~ dunif(0.0001,",k.max,")
      L0 ~ dnorm(",L0,",",L0.se,"^(-2)) T(0,",max(Length),")
      sigma ~ dunif(0,",sigma.max,")
      }")
  } else{
    stop("Model must be specified as either'VB', 'Log' or 'Gom'")
  }

  Growth_jags <- rjags::jags.model(textConnection(Growth_model),
                                   data = list(Age = Age, Length = Length),
                                   n.chains = n.chains, n.adapt = BurnIn,
                                   quiet=TRUE)


  Growth_sim <- rjags::coda.samples(model = Growth_jags,
                                    variable.names = c("Linf", "k", "L0", "sigma"),
                                    n.iter = iter,
                                    thin = thin)


  return(Growth_sim)

}


#' Calculate_DIC
#' @description Deviance Information Criteria (DIC) are calculated for three growth models
#'     (von Bertalanffy, Gompertz and Logistic) using the same prior parameters for each. These
#'     are used for model selection.
#' @param data A data.frame that contains columns named 'Age' and "Length'. The function can
#'     detect columns with similar names. If age and length columns cannot be determined then
#'     an error will occur. The dataset can have additional columns which will be ignored by
#'     the function
#' @param Linf The prior for asymptotic length. MUst be in the same unit (i.e. cm or mm) as the data.
#'     This should be based off of maximum size for the species.
#' @param Linf.se The prior for normally distributed standard error around
#'     asymptotic length. MUst be in the same unit (i.e. cm or mm) as the data. Cannot be zero.
#' @param L0 The prior for length-at-birth. MUst be in the same unit (i.e. cm or mm) as the data.
#'     This should be based off of minimum size for the species.
#' @param L0.se The prior for normally distributed standard error around
#'     length-at-birth. MUst be in the same unit (i.e. cm or mm) as the data. Cannot be zero.
#' @param k.max The maximum value to consider for the growth completion parameter 'k'. In
#'     the Gompertz and Logistic models, this parameter is often notated as 'g' instead of 'k'.
#' @param sigma.max The maximum value to consider for sigma. This is the variance around the
#'     length-at-age residuals.
#' @param iter How many MCMC iterations should be run? Default is 10000 but fewer can be useful to
#'     avoid longer run times when testing code or data
#' @param BurnIn The number of iterations at the begining of each chain to discard ('Burn in') to
#'     avoid biased values from starting values that do not resemble the target distribution.
#'     Default is 1000.
#' @param n.chains Number of MCMC chains to be run. Default is 4.
#' @param thin The thinning of the MCMC simulations. Deafault is 1 which means no thinning occurs.
#'     Thinning is generally only necessary for complicated models as it increases run time.
#'
#' @return A data.frame with the results of the DIC analysis for each model.
#' @export
Calculate_DIC <- function(data, Linf = NULL, Linf.se = NULL,
                          L0 = NULL, L0.se = NULL, k.max = NULL, sigma.max = NULL,  iter = 10000,BurnIn = 1000,
                          n.chains = 4, thin = iter/10){


  age_col <- grep("age", substr(tolower(names(data)),1,3))
  if(length(age_col) <1) stop("Age column heading could not be distinguished ")
  if(length(age_col) >1) stop("Multiple age columns detected. Remove unecessary variables or rename desired column to 'Age' ")
  len_col <- grep("len|tl|lt|siz", substr(tolower(names(data)),1,3))
  if(length(len_col) <1) stop("Length column heading could not be distinguished ")
  if(length(len_col) >1) stop("Multiple length columns detected. Remove unecessary variables or rename desired column to 'Length' ")

  Age <- data[,age_col]
  Length <- data[,len_col]



  VB_model <- paste0(
    "model{
    # Likelihood model for Length[i]
    for(i in 1:length(Length)) {
    Length[i] ~ dnorm(Lt[i], sigma^(-2))
    Lt[i] <- Linf-(Linf-L0)*exp(-k*Age[i])
    }

    # Prior parameters
    Linf ~ dnorm(",Linf,",",Linf.se,"^(-2))
    k ~ dunif(0,",k.max,")#,
    # L0_dist,
    L0 ~ dnorm(",L0,",",L0.se,"^(-2)) T(0,",max(Length),")
    sigma ~ dunif(0,",sigma.max,")
}")

  Gom_model <- paste0(
    "model{
    # Likelihood model for Length[i]
    for(i in 1:length(Length)) {
    Length[i] ~ dnorm(Lt[i], sigma^(-2))
    Lt[i] <-L0*exp(log(Linf/L0)*(1-exp(-k*Age[i])))
    }

    # Prior parameters
    Linf ~ dnorm(",Linf,",",Linf.se,"^(-2))
    k ~ dunif(0.0001,",k.max,")
    L0 ~ dnorm(",L0,",",L0.se,"^(-2)) T(0,",max(Length),")
    sigma ~ dunif(0,",sigma.max,")
    }")

  Log_model <- paste0(
    "model{
    # Likelihood model for Length[i]
    for(i in 1:length(Length)) {
    Length[i] ~ dnorm(Lt[i], sigma^(-2))
    Lt[i] <-(Linf*L0*exp(k*Age[i]))/(Linf+L0*(exp(k*Age[i])-1))
    }

    # Prior parameters
    Linf ~ dnorm(",Linf,",",Linf.se,"^(-2))
    k ~ dunif(0.0001,",k.max,")
    L0 ~ dnorm(",L0,",",L0.se,"^(-2)) T(0,",max(Length),")
    sigma ~ dunif(0,",sigma.max,")
    }")

  VB_jags <- rjags::jags.model(textConnection(VB_model),
                               data = list(Age = Age, Length = Length),
                               n.chains = n.chains, n.adapt = BurnIn,
                               quiet=TRUE)
  message("Running von Bertalanffy model \n")
  VB_dic <- rjags::dic.samples(VB_jags, n.iter = iter, thin = thin, type = "pD")

  Gom_jags <- rjags::jags.model(textConnection(Gom_model),
                                data = list(Age = Age, Length = Length),
                                n.chains = n.chains,n.adapt = BurnIn,
                                quiet=TRUE)
  message("Running Gompterz model \n")
  Gom_dic <- rjags::dic.samples(Gom_jags, n.iter = iter, thin = thin, type = "pD")


  Log_jags <- rjags::jags.model(textConnection(Log_model),
                                data = list(Age = Age, Length = Length),
                                n.chains = n.chains, n.adapt = BurnIn,
                                quiet=TRUE)
  message("Running Logistic model \n")
  Log_dic <- rjags::dic.samples(Log_jags, n.iter = iter, thin = thin, type = "pD")



  results <- data.frame(Model = c("VB", "Gom", "Log"),
                        `Mean.deviance` = c(sum(VB_dic$deviance), sum(Gom_dic$deviance), sum(Log_dic$deviance)),
                        penalty = c(sum(VB_dic$penalty), sum(Gom_dic$penalty), sum(Log_dic$penalty)))
  results$`Penalised.deviance` <- results$`Mean.deviance`+results$penalty
  colnames(results)<- c("Model", "Mean deviance", "pD", "Penalised deviance")
  results <- dplyr::arrange(results, `Penalised deviance`)

  return(results)

}

#' Get_MCMC_parameters
#' @description Get parameter summary statistics from the outputs of a Estimate_MCMC_Growth object
#' @param data n output from the Estimate_MCMC_Growth funtion or a similarly structure rjags output.
#'     Must be of class "mcmc.list"
#'
#' @return A tibble with 5 columns that include each parameter and the results of their posterior distributions.
#'     These include the mean, Standard deviation, lower 95th percentile and upper 95th percentile
#' @import tidyr
#' @importFrom stats sd quantile
#' @export

Get_MCMC_parameters <-function(data){
  if(class(data) != "mcmc.list") stop("`data` must be a result returned from `Estimate_MCMC_Growth()`")

  processed_data <- plyr::ldply(data)
  processed_data <- tidyr::gather(processed_data, Parameter, Value)
  processed_data <- dplyr::group_by(processed_data, Parameter)
  results <- dplyr::summarise(processed_data, mean = mean(Value),
                              SD = sd(Value),
                              `2.5%` =  quantile(Value, probs = 0.025),
                              `97.5%` =  quantile(Value, probs = 0.975))

  return(results)

}
#' Calc_Logistic_LAA
#'
#' @param Linf A single value of asymptotic length for the logistic model
#' @param k A single value of the growth completion parameter for the logistic model
#' @param L0 A single value of length-at-birth for the logistic model
#' @param Age A single value or vector of ages to convert to length based on the logistic model
#'
#' @return A vector of length-at-ages
#' @export
Calc_Logistic_LAA <- function(Linf, k, L0, Age){
  LAA <- (Linf*L0*exp(k*Age))/(Linf+L0*(exp(k*Age)-1))
  return(LAA)
}

#' Calc_VBGF_LAA
#'
#' @param Linf A single value of asymptotic length for the von Bertalanffy model
#' @param k A single value of the growth completion parameter for the von Bertalanffy model
#' @param L0 A single value of length-at-birth for the von Bertalanffy model
#' @param Age A single value or vector of ages to convert to length based on the von Bertalanffy model
#'
#' @return A vector of length-at-ages
#' @export
Calc_VBGF_LAA <- function(Linf, k, L0, Age){
  LAA <- Linf-(Linf-L0)*exp(-k*Age)
  return(LAA)
}

#' Calc_Gompertz_LAA
#'
#' @param Linf A single value of asymptotic length for the Gompertz model
#' @param k A single value of the growth completion parameter for the Gompertz model
#' @param L0 A single value of length-at-birth for the Gompertz model
#' @param Age A single value or vector of ages to convert to length based on the Gompertz model
#'
#' @return A vector of length-at-ages
#' @export
Calc_Gompertz_LAA <- function(Linf, k, L0, Age){
  LAA <-L0*exp(log(Linf/L0)*(1-exp(-k*Age)))
  return(LAA)
}



#' Calculate_MCMC_growth_curve
#' @description A mcmc.list object produced from Estimate_MCMC_Growth is converted to a dataframe
#'     and structured using requested quantiles. This function takes the list of MCMC results for multiple
#'     chains, restructures them into a dataframe and calculates quantiles around length-at-age estimates.
#'     The quantiles are produced using the tidybayes::mean_qi() function and this result is returned from the function.
#'     This can be conveniently plotted in a ggplot using the geom_lineribbon() function provided in the tidybayes
#'     package.
#'
#' @param data An output from the Estimate_MCMC_Growth funtion or a similarly structure rjags output.
#'     Must be of class "mcmc.list"
#' @param Model The model used in the Estimate_MCMC_Growth object. Either "VB", "Gom" or "Log".
#' @param max.age The max age to estimate growth up until.
#' @param probs The percentiles of the results to return. Can be a single value or a vector of
#'     values. A single quantile width is required rather than its range. For example, 50th
#'     percentiles would be width = .5 which would return a lower percentile at .25 and an upper
#'     percentile of .75.
#' @importFrom plyr ldply
#' @import tidybayes
#'
#' @return A tibble that has been formatted using  tidybayes::mean_qi(). This includes
#'     variables: Age, LAA, .lower, .upper, .width, .point and .interval.
#' @export
Calculate_MCMC_growth_curve <- function(data, Model = NULL, max.age = NULL, probs = c(0.5,0.75,0.95)){

  if(class(data) != "mcmc.list") stop("`data` must be a result returned from `Estimate_MCMC_Growth()`")
  if(!Model %in% c("VB", "Gom", "Log")) stop("'Model must be one of either 'VB', 'Gom' or 'Log")
  if(is.null(max.age)) stop("Please specify max age")

  processed_data <- plyr::ldply(data)
  processed_data <- dplyr::select(processed_data, -sigma)
  processed_data <- dplyr::mutate(processed_data, sim = row_number())
  processed_data <- dplyr::left_join(processed_data, expand.grid(sim = processed_data$sim, Age = seq(0,max.age, 0.1)), by = "sim")
  processed_data <- dplyr::mutate(processed_data, LAA = dplyr::case_when(
    Model == "VB" ~ Calc_VBGF_LAA(Linf, k, L0, Age),
    Model == "Log" ~ Calc_Logistic_LAA(Linf, k, L0, Age),
    Model == "Gom" ~ Calc_Gompertz_LAA(Linf, k, L0, Age),
    TRUE ~ NA_real_
  ))
  processed_data <- dplyr::group_by(processed_data, Age)
  results <-tidybayes::mean_qi(processed_data, LAA,.width = probs)

  return(results)
}







