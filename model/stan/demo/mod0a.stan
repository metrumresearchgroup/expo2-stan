
data {
  int<lower=0> N ;
  int<lower=0> n_id ;
  vector[N] DV;
  vector[N] Conc;
  array[N] int ID;
  array[n_id] int start;
  array[n_id] int end;
}

parameters {
  real<lower=0, upper=100> emax;
  real tv_log_ec50;
  real tv_e0;
  real log_gamma;
  
  vector[n_id] log_ec50;
  vector[n_id] e0;
  
  real<lower=0> sigma;
  real<lower=0> omega_e0;
  real<lower=0> omega_log_ec50;
}

transformed parameters {
  real tv_ec50 = exp(tv_log_ec50);
  real gamma = exp(log_gamma);
  vector[n_id] ec50;
  vector[N] mu;
  
  ec50 = exp(log_ec50);

  for (i in 1:N) {
    mu[i] = e0[ID[i]] + (emax-e0[ID[i]])*pow(Conc[i]/100,gamma) / (pow(ec50[ID[i]]/100,gamma) + pow(Conc[i]/100,gamma));
  }
  
}

model {
  // Priors
  tv_e0 ~ normal(0, 10);
  tv_log_ec50 ~ normal(4,2);
  emax ~ uniform(0,100);
  log_gamma ~ normal(0,1);
  sigma ~ normal(0,10);
  omega_e0 ~ normal(0,2);
  omega_log_ec50 ~ normal(0,1);

  e0 ~ normal(tv_e0, omega_e0);
  log_ec50 ~ normal(tv_log_ec50, omega_log_ec50);

// Likelihhood

  DV ~ normal(mu, sigma);

}

