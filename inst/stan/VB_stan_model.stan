
data {

  int<lower=1> n; //number of samples
  //data vectors
  vector<lower=0>[n] Age; //age data
  vector<lower=0>[n] Length; //length data

  //prior data
  vector[4] priors; //Linf, L0, k, sigma
  vector<lower=0>[2] priors_se; //sd of Linf, L0

}

parameters {
  //VBGM parameters
  real<lower=0> L0; //length-at-birth
  real<lower=0> Linf; //asymptotic length
  real<lower=0> k; // growth coefficient

  //Likelihood parameters
  real<lower=0> sigma; //RSE
}

model {
  //storage
  vector[n] PredL; //predicted lengths

  //VBGM priors
  Linf ~ normal(priors[1], priors_se[1]);
  L0 ~ normal(priors[2], priors_se[2]);
  k ~ uniform(0, priors[3]);

  sigma ~ uniform(0, priors[4]);


  //VBGM likelihood
  for(i in 1:n){
    PredL[i] = Linf - (Linf - L0)*exp(-k*Age[i]);
    target += normal_lpdf(Length[i]|PredL[i], sigma);//likelihood

  }

}
// Individual loglikelihoods for loo


generated quantities {
    vector[n] log_lik;
    for (i in 1:n) {
      log_lik[i] = normal_lpdf(Length[i]|(Linf - (Linf - L0)*exp(-k*Age[i])), sigma);
    }

}
