data {
  int<lower=1> N; // Total number of observations
  int<lower=1> J; // Number of individuals
  vector[N] age; // Age of each observation
  vector[N] length; // Length of each observation
  int<lower=1, upper=J> individual[N]; // Individual index for each observation

  //prior data
  vector[7] priors; //Linf, L0, k, sigma, and the sigma of Linf, L0 and k
  vector<lower=0>[2] priors_se; //sd of Linf, L0

}

parameters {
  real<lower=0> Linf; // Population-level asymptotic length
  real<lower=0> k; // Population-level growth rate
  real L0; // Population-level theoretical length at age  zero

  vector<lower=0>[J] L_inf_j; // Individual-level asymptotic length
  vector<lower=0>[J] K_j; // Individual-level growth rate
  vector[J] L0_j; // Individual-level theoretical length at age  zero

  real<lower=0> sigma_L_inf; // SD of individual-level L_inf
  real<lower=0> sigma_K; // SD of individual-level K
  real<lower=0> sigma_L0; // SD of individual-level t0
  real<lower=0> sigma; // SD of measurement error
}

model {
  // Priors
  Linf ~ normal(priors[1], priors_se[1]);
  k ~ uniform(0, priors[3]);
  L0 ~ normal(priors[2], priors_se[2]);

  L_inf_j ~ normal(Linf, sigma_L_inf);
  K_j ~ normal(k, sigma_K);
  L0_j ~ normal(L0, sigma_L0);

  sigma_L_inf ~ uniform(0, priors[5]);
  sigma_K ~ uniform(0, priors[7]);
  sigma_L0 ~ uniform(0, priors[6]);
  sigma ~ uniform(0, priors[4]);

  // Likelihood
  for (i in 1:N) {
    target += normal_lpdf(length[i] | L0_j[individual[i]]*exp(log(L_inf_j[individual[i]]/L0_j[individual[i]])*(1-exp(-K_j[individual[i]]*age[i]))), sigma);
  }

}

  // log likelihood for loo
  generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = normal_lpdf(length[i] | L0_j[individual[i]]*exp(log(L_inf_j[individual[i]]/L0_j[individual[i]])*(1-exp(-K_j[individual[i]]*age[i]))), sigma);
    }

}
