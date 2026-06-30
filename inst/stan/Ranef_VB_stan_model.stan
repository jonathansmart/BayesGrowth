
data {
  int<lower=1> N; // Total number of observations
  int<lower=1> J; // Number of Groups
  vector[N] Age; // Age of each observation
  vector[N] Length; // Length of each observation
  int<lower=1, upper=J> Group[N]; // Group index for each observation

  //prior data
  vector[7] priors; //Linf, L0, k, sigma, and the sigma of Linf, L0 and k
  vector<lower=0>[2] priors_se; //sd of Linf, L0

}

parameters {
  real<lower=0> Linf; // Population-level asymptotic Length
  real<lower=0> k; // Population-level growth rate
  real L0; // Population-level theoretical Length at Age  zero

  vector<lower=0>[J] Linf_j; // Group-level asymptotic Length
  vector<lower=0>[J] k_j; // Group-level growth rate
  vector[J] L0_j; // Group-level theoretical Length at Age  zero

  real<lower=0> sigma_Linf; // SD of Group-level Linf
  real<lower=0> sigma_k; // SD of Group-level k
  real<lower=0> sigma_L0; // SD of Group-level t0
  real<lower=0> sigma; // SD of measurement error
}

model {
  // Priors
  Linf ~ normal(priors[1], priors_se[1]);
  k ~ uniform(0, priors[3]);
  L0 ~ normal(priors[2], priors_se[2]);

  Linf_j ~ normal(Linf, sigma_Linf);
  k_j ~ normal(k, sigma_k);
  L0_j ~ normal(L0, sigma_L0);

 sigma_Linf ~ uniform(0, priors[5]);
  sigma_k ~ uniform(0, priors[7]);
  sigma_L0 ~ uniform(0, priors[6]);
  sigma ~ uniform(0, priors[4]);

  // Likelihood
  for (i in 1:N) {
    target += normal_lpdf(Length[i] | Linf_j[Group[i]] - (Linf_j[Group[i]] - L0_j[Group[i]])*exp(-k_j[Group[i]]*Age[i]), sigma);

  }

}

  // log likelihood for loo
  generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = normal_lpdf(Length[i] | Linf_j[Group[i]] - (Linf_j[Group[i]] - L0_j[Group[i]])*exp(-k_j[Group[i]]*Age[i]), sigma);
    }

  }
