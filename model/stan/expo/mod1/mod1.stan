
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
  
  array[n_id] vector[2] etas; // etas[i,1] = log_ec50; etas[i,2] = e0;

  real<lower=0> sigma;
  vector<lower=0>[2] omega_sds; // SDs for IIV
  corr_matrix[2] Omega_corr;
}

transformed parameters {
  real tv_ec50 = exp(tv_log_ec50);
  real gamma = exp(log_gamma);
  vector[n_id] ec50;
  vector[n_id] e0;
  vector[N] mu;
  
  
  for (i in 1:n_id) {
    ec50[i] = exp(tv_log_ec50 + etas[i,1]);
    e0[i] = tv_e0 + etas[i,2];
  }
  
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
  
  Omega_corr ~ lkj_corr(2);
  omega_sds ~ cauchy(0,2.5);

  etas ~ multi_normal([0,0],quad_form_diag(Omega_corr, omega_sds));

// Likelihhood
  
  DV ~ normal(mu, sigma);

}


generated quantities {

  vector[N] simdv_obs;
  vector[N] simdv_new;

  // Simulated observations for observed subjects
  for (i in 1:N) {
    simdv_obs[i] = normal_rng(mu[i],sigma);
  }
  
  // Simulate observations for new subject
  {
    array[n_id] vector[2] etas_new; // etas[i,1] = log_ec50; etas[i,2] = e0;
    vector[n_id] ec50_new;
    vector[n_id] e0_new;
    vector[N] mu_new;
    
    for (i in 1:n_id) {
      etas_new[i,] = multi_normal_rng([0.0,0.0], quad_form_diag(Omega_corr, omega_sds));
        ec50_new[i] = exp(tv_log_ec50 + etas_new[i,1]);
        e0_new[i] = tv_e0 + etas_new[i,2];
    }
    
    for (i in 1:N) {
      mu_new[i] = e0_new[ID[i]] + (emax-e0_new[ID[i]])*pow(Conc[i]/100,gamma) / (pow(ec50_new[ID[i]]/100,gamma) + pow(Conc[i]/100,gamma));
      simdv_new[i] = normal_rng(mu_new[i] , sigma);
    }
  }
  
}
