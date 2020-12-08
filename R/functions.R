#' Estimate_MCMC_Growth
#' @description A wrapper function that creates a Stan MCMC model using the rstan package.
#'     The data and priors provided are combined into an rstan model that estimates a length-at-age
#'     model with a normal distribution. Three different growth models can be used: a von Bertalanffy
#'     model, Gompertz model or a logistic model. The prior on Linf and L0 are normally distributed
#'     and determined through the user providing a mean and se for each parameter. The growth completion
#'     parameter for any model (k) has a uniform prior which only requires an upper bound with the
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
#' @param BurnIn The number of iterations at the beginning of each chain to discard ('Burn in') to
#'     avoid biased values from starting values that do not resemble the target distribution.
#'     Default is 1000.
#' @param n.chains Number of MCMC chains to be run. Default is 4.
#' @param controls A named list of parameters to control the rstan models behaviour.
#' @param thin The thinning of the MCMC simulations. Default is 1 which means no thinning occurs.
#'     Thinning is generally only necessary for complicated models as it increases run time.
#' @param n_cores The number of cores to be used for parallel processing. It should be 1 core less than the
#'     maximum number available.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate output from Stan on the console,
#'     which might be helpful for model debugging.
#' @import dplyr rstan
#'
#' @return An object of class 'stanfit' from the rstan package.
#' @export
Estimate_MCMC_Growth <- function(data,  Model = NULL, Linf = NULL, Linf.se = NULL,
                                 L0 = NULL, L0.se = NULL, k.max = NULL, sigma.max = NULL,
                                 iter = 10000, BurnIn = iter*0.1, n_cores = 1, controls = NULL,
                                 n.chains = 4, thin = 1,verbose = FALSE){

  if(any(is.null(c(Linf, Linf.se, L0, L0.se, k.max, sigma.max)))) stop("At least one parameter or its error are not correctly specified")
  if(length(Model) != 1) stop("Only one growth model can be used in each function call")
  if(is.null(Model))stop("Growth model has not been specified")
  if(!Model %in% c("VB", "Gom", "Log")) stop("Model must be specified as either'VB', 'Log' or 'Gom'")


  age_col <- grep("age", substr(tolower(names(data)),1,3))
  if(length(age_col) <1) stop("Age column heading could not be distinguished ")
  if(length(age_col) >1) stop("Multiple age columns detected. Remove unecessary variables or rename desired column to 'Age' ")

  len_col <- grep("len|tl|lt|siz", substr(tolower(names(data)),1,3))
  if(length(len_col) <1) stop("Length column heading could not be distinguished ")
  if(length(len_col) >1) stop("Multiple length columns detected. Remove unecessary variables or rename desired column to 'Length' ")


  if(Linf.se == 0 | L0.se == 0) stop("L0 and Linf standard error priors cannot be zero")
  if(any(is.na(data))) stop("data contains NA's")

  if(n_cores >  parallel::detectCores()-1) {
    n_cores <- 1
    message("Not enough cores available. Reseting to 1 core")
  }

  if(is.null(controls)) controls <- list(adapt_delta = 0.9)

  Age <- data[,age_col]
  Length <- data[,len_col]

  starting_parameters <- function(chain_id){

    mean.age<-tapply(Length, round(Age), mean,na.rm = T)
    Lt1<-mean.age[2:length(mean.age)]
    Lt<-mean.age[1:length(mean.age)-1]
    model<-lm(Lt1 ~ Lt)
    k <- suppressWarnings(abs(-log(model$coef[2]))) #in case a nan occurs
    k <- ifelse(is.nan(k),0.1,k) # in case a nan occurs
    Linf<-abs(model$coef[1]/(1-model$coef[2]))

    L0<-lm(mean.age ~ poly(as.numeric(names(mean.age)), 2, raw = TRUE))$coef[1]

    return(list(Linf = Linf, L0 = L0, k = k, sigma = sigma.max/2))
  }

  if(verbose == FALSE){
    text <- 0
  }else{
    text <- iter/10
  }

  priors <- c(Linf, L0, k.max, sigma.max)
  priors_se <- c(Linf.se, L0.se)

  dat <- list(n = length(Age),
              Age = Age,
              Length = Length,
              priors = priors,
              priors_se = priors_se)

  if(Model == "VB"){
    Growth_model <- rstan::sampling(object = stanmodels$VB_stan_model,
                                    data = dat,
                                    init = starting_parameters,
                                    control = controls,
                                    warmup = BurnIn,
                                    thin = thin,
                                    verbose = verbose,
                                    iter = iter,
                                    cores = n_cores,
                                    open_progress = FALSE,
                                    refresh = text,
                                    # model_name = "Von Bertalanffy",
                                    include = TRUE,
                                    pars = c("Linf", "k","L0", "sigma"),
                                    chains=n.chains)

  } else if(Model == "Gom"){
    Growth_model <- rstan::sampling(object = stanmodels$Gompertz_stan_model,
                                    data = dat,
                                    init = starting_parameters,
                                    control = controls,
                                    warmup = BurnIn,
                                    thin = thin,
                                    verbose = verbose,
                                    iter = iter,
                                    cores = n_cores,
                                    open_progress = FALSE,
                                    refresh = text,
                                    include = TRUE,
                                    pars = c("Linf", "k","L0", "sigma"),
                                    chains=n.chains)

  } else if(Model == "Log"){
    Growth_model <- rstan::sampling(object = stanmodels$Logistic_stan_model,
                                    data = dat,
                                    init = starting_parameters,
                                    control = controls,
                                    warmup = BurnIn,
                                    thin = thin,
                                    verbose = verbose,
                                    iter = iter,
                                    open_progress = FALSE,
                                    refresh = text,
                                    cores = n_cores,
                                    include = TRUE,
                                    pars = c("Linf", "k","L0", "sigma"),
                                    chains=n.chains)


  } else{
    stop("Model must be specified as either'VB', 'Log' or 'Gom'")
  }




  return(Growth_model)

}

#' Compare_Growth_Models
#' @description Conduct growth model selection using 'Leave One Out' (LOO) cross validation analysis and
#'     Widely Applicable Information Criterion (WAIC)for three growth models:
#'     (von Bertalanffy, Gompertz and Logistic) using the same prior parameters for each.
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
#' @param BurnIn The number of iterations at the beginning of each chain to discard ('Burn in') to
#'     avoid biased values from starting values that do not resemble the target distribution.
#'     Default is 1000.
#' @param n.chains Number of MCMC chains to be run. Default is 4.
#' @param controls A named list of parameters to control the rstan models behaviour.
#' @param thin The thinning of the MCMC simulations. Default is 1 which means no thinning occurs.
#'     Thinning is generally only necessary for complicated models as it increases run time.
#' @param n_cores The number of cores to be used for parallel processing. It should be 1 core less than the
#'     maximum number available.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate output from Stan on the console,
#'     which might be helpful for model debugging.
#' @param stats Which statistics should be returned: LooIC, WAIC or both (both will return a list)
#' @import dplyr rstan loo
#'
#' @return A dataframe with the requested stats
#' @export
Compare_Growth_Models <- function(data,   Linf = NULL, Linf.se = NULL,
                                  L0 = NULL, L0.se = NULL, k.max = NULL, sigma.max = NULL,
                                  iter = 10000, BurnIn = iter*0.1, n_cores = 1,controls = NULL,
                                  n.chains = 4, thin = 1,verbose = FALSE, stats = "LooIC"){


  if(any(is.null(c(Linf, Linf.se, L0, L0.se, k.max, sigma.max)))) stop("At least one parameter or its error are not correctly specified")


  age_col <- grep("age", substr(tolower(names(data)),1,3))
  if(length(age_col) <1) stop("Age column heading could not be distinguished ")
  if(length(age_col) >1) stop("Multiple age columns detected. Remove unecessary variables or rename desired column to 'Age' ")

  len_col <- grep("len|tl|lt|siz", substr(tolower(names(data)),1,3))
  if(length(len_col) <1) stop("Length column heading could not be distinguished ")
  if(length(len_col) >1) stop("Multiple length columns detected. Remove unecessary variables or rename desired column to 'Length' ")


  if(Linf.se == 0 | L0.se == 0) stop("L0 and Linf standard error priors cannot be zero")
  if(any(is.na(data))) stop("data contains NA's")

  if(n_cores >  parallel::detectCores()-1) {
    n_cores <- 1
    message("Not enough cores available. Reseting to 1 core")
  }

  if(is.null(controls)) controls <- list(adapt_delta = 0.9)

  if(!stats %in% c("LooIC", "WAIC", "both")) stop("stats must be specified as either'LooIC', 'WAIC' or 'both'")

  Age <- data[,age_col]
  Length <- data[,len_col]

  starting_parameters <- function(chain_id){

    mean.age<-tapply(Length, round(Age), mean,na.rm = T)
    Lt1<-mean.age[2:length(mean.age)]
    Lt<-mean.age[1:length(mean.age)-1]
    model<-lm(Lt1 ~ Lt)
    k <- suppressWarnings(abs(-log(model$coef[2]))) #in case a nan occurs
    k <- ifelse(is.nan(k),0.1,k) # in case a nan occurs
    Linf<-abs(model$coef[1]/(1-model$coef[2]))

    L0<-lm(mean.age ~ poly(as.numeric(names(mean.age)), 2, raw = TRUE))$coef[1]

    return(list(Linf = Linf, L0 = L0, k = k, sigma = sigma.max/2))
  }

  if(verbose == FALSE){
    text <- 0
  }else{
    text <- iter/10
  }

  priors <- c(Linf, L0, k.max, sigma.max)
  priors_se <- c(Linf.se, L0.se)

  dat <- list(n = length(Age),
              Age = Age,
              Length = Length,
              priors = priors,
              priors_se = priors_se)

  VB_model <-
    rstan::sampling(object = stanmodels$VB_stan_model,
                    data = dat,
                    init = starting_parameters,
                    control = controls,
                    warmup = BurnIn,
                    thin = thin,
                    verbose = verbose,
                    open_progress = FALSE,
                    # # refresh = 0,
                    refresh = text,
                    iter = iter,
                    cores = n_cores,

                    chains=n.chains)



  Gom_model <- rstan::sampling(object = stanmodels$Gompertz_stan_model,
                               data = dat,
                               init = starting_parameters,
                               control = controls,
                               warmup = BurnIn,
                               thin = thin,
                               verbose = verbose,
                               open_progress = FALSE,
                               # refresh = 0,
                               refresh = text,
                               iter = iter,
                               cores = n_cores,
                               chains=n.chains)


  Logistic_model <- rstan::sampling(object = stanmodels$Logistic_stan_model,
                                    data = dat,
                                    init = starting_parameters,
                                    control = controls,
                                    warmup = BurnIn,
                                    thin = thin,
                                    verbose = verbose,
                                    open_progress = FALSE,
                                    # refresh = 0,
                                    refresh = text,
                                    iter = iter,
                                    cores = n_cores,
                                    chains=n.chains)

  # Calculate loo
  VB_loo <- suppressWarnings(loo::loo(VB_model))
  Gom_loo <- suppressWarnings(loo::loo(Gom_model))
  Logistic_loo <- suppressWarnings(loo::loo(Logistic_model))

  # Get loo comparions
  Loo_comp <- as.data.frame(loo::loo_compare(list(VB = VB_loo, Gompertz = Gom_loo, Logistic = Logistic_loo)))
  Loo_comp <- tibble::rownames_to_column(Loo_comp, "Model")

  # Get looic weights
  model_list <- list(VB = VB_model, Gompertz =Gom_model, Logistic =Logistic_model)
  log_lik_list <- lapply(model_list, loo::extract_log_lik)

  looiW <- as.data.frame(
    round(
      as.matrix(
        suppressWarnings(loo::loo_model_weights(
          log_lik_list,
          method = "stacking",
          optim_control = list(reltol=1e-10))
        ))
      ,2)
  )

  colnames(looiW) <- "looic_Weight"
  looiW <- tibble::rownames_to_column(looiW,"Model")

  # combine all Looic results
  Loo_results<- dplyr::left_join(Loo_comp, looiW,by = "Model")

  # get waics
  waic_VB <- suppressWarnings(loo::waic(loo::extract_log_lik(VB_model)))
  waic_Gom <- suppressWarnings(loo::waic(loo::extract_log_lik(Gom_model)))
  waic_Log <- suppressWarnings(loo::waic(loo::extract_log_lik(Logistic_model)))

  # get Waic Weights
  waics <- c(
    waic_VB$estimates["elpd_waic", 1],
    waic_Gom$estimates["elpd_waic", 1],
    waic_Log$estimates["elpd_waic", 1]
  )

  # Get p_waics
  p_waic <- round(c(
    waic_VB$estimates["p_waic", 1],
    waic_Gom$estimates["p_waic", 1],
    waic_Log$estimates["p_waic", 1]
  ),1)

  # Combine WAIC results
  waic_results <- data.frame(Model = c("VB", "Gompertz", "Logistic"),WAIC = waics,p_waic = p_waic,`WAIC_weight`=round(waics/sum(waics),2))

  if(stats == "both"){
    Results <-  list(LooIC = Loo_results, WAIC = waic_results)
  } else if(stats == "LooIC") {
    Results <-  Loo_results
  }else if(stats == "WAIC") {
    Results <-  waic_results
  }


  return(Results)

}




#' Get_MCMC_parameters
#' @description Get parameter summary statistics from the outputs of a Estimate_MCMC_Growth object. It is simplified set of
#'     results than is returned from summary(obj).
#' @param obj An output from the Estimate_MCMC_Growth function
#' @return A data.frame with the posterior distributions for each parameter.
#'     These include the mean, Standard error of the mean, Standard deviationof the mean, median,
#'     95th percentiles, effective sample sizes and Rhat.
#' @import tibble
#' @export
#'
Get_MCMC_parameters <- function (obj)
{
  if (class(obj) != "stanfit")
    stop("`obj` must be a result returned from `Estimate_MCMC_Growth()`")
  results <- as.data.frame(summary(obj, pars = c("Linf", "k",
                                                 "L0", "sigma"), probs = c(0.025, 0.5, 0.975))$summary)

  results <- tibble::rownames_to_column(results, var = "Parameter")
  results <- dplyr::mutate_at(results,.vars = -everything("Parameter"), .funs = ~round(.,2))
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
#' @description A 'stan.fit' object produced from Estimate_MCMC_Growth is converted to a dataframe
#'     and structured using requested quantiles. This function takes the list of MCMC results for multiple
#'     chains, restructures them into a dataframe and calculates quantiles around length-at-age estimates.
#'     The quantiles are produced using the tidybayes::mean_qi() function and this result is returned from the function.
#'     This can be conveniently plotted in a ggplot using the geom_lineribbon() function provided in the tidybayes
#'     package.
#' @param data An output from the Estimate_MCMC_Growth function
#' @param Model The model used in the Estimate_MCMC_Growth object. Either "VB", "Gom" or "Log".
#' @param max.age The max age to estimate growth up until.
#' @param probs The percentiles of the results to return. Can be a single value or a vector of
#'     values. A single quantile width is required rather than its range. For example, 50th
#'     percentiles would be width = .5 which would return a lower percentile at .25 and an upper
#'     percentile of .75.
#' @import tidybayes
#'
#' @return A tibble that has been formatted using  tidybayes::mean_qi(). This includes
#'     variables: Age, LAA, .lower, .upper, .width, .point and .interval.
#' @export
Calculate_MCMC_growth_curve <- function(obj, Model = NULL, max.age = NULL, probs = c(0.5,0.75,0.95)){

  if(!Model %in% c("VB", "Gom", "Log")) stop("'Model must be one of either 'VB', 'Gom' or 'Log")
  if(is.null(max.age)) stop("Please specify max age")

  processed_data <- Get_MCMC_parameters(obj)
  L0_sims <- rstan::extract(obj)$L0
  Linf_sims <- rstan::extract(obj)$Linf
  k_sims <- rstan::extract(obj)$k
  processed_data <- data.frame(Linf = Linf_sims,k =  k_sims,L0 = L0_sims)

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
