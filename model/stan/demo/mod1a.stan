
data {
  int<lower=0> N ;
  int<lower=0> n_id ;
  vector[N] DV;
  vector[N] Conc;
  array[N] int ID;
}

parameters {
  real<lower=0, upper=100> emax;
  real tv_log_ec50;
  real tv_e0;
  real log_gamma;
  
  matrix[2,n_id] ipars; // ipars[1,i] = log_ec50; ipars[2,i] = e0;

  real<lower=0> sigma;
  vector<lower=0>[2] omega_sds; // SDs for IIV
  corr_matrix[2] Omega_corr;
}

transformed parameters {
  real tv_ec50 = exp(tv_log_ec50);
  real gamma = exp(log_gamma);
  vector[n_id] ec50;
  vector[n_id] e0;
  
  for (i in 1:n_id) {
    ec50[i] = exp(ipars[1,i]);
    e0[i] = ipars[2,i];
  }
}

model {
  // Priors
  tv_e0 ~ normal(0, 10);
  tv_log_ec50 ~ normal(4,2);
  emax ~ uniform(0,100);
  log_gamma ~ normal(0,1);
  sigma ~ normal(0,10);
  
  Omega_corr ~ lkj_corr(2);
  omega_sds ~ cauchy(0,2.5);

  for (i in 1:n_id) {
    ipars[,i] ~ multi_normal([tv_log_ec50,tv_e0],quad_form_diag(Omega_corr, omega_sds));
  }

// Likelihhood

  vector[N] mu;
  
  for (i in 1:N) {
    mu[i] = e0[ID[i]] + (emax-e0[ID[i]])*pow(Conc[i]/100,gamma) / (pow(ec50[ID[i]]/100,gamma) + pow(Conc[i]/100,gamma));
  }
  
  DV ~ normal(mu, sigma);

}
