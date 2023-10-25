
data {
  int<lower=0> N ;
  int<lower=0> n_id ;
  vector[N] DV;
  vector[N] Conc;
  array[N] int ID;
  array[n_id] int<lower=1, upper=N> id_start;
  array[n_id] int<lower=1, upper=N> id_end;
  
  int<lower=0,upper=1> prior_only;
}


parameters {
  real<lower=0, upper=100> emax;
  real tv_log_ec50;
  real tv_e0;
  real log_gamma;
  
  matrix[2,n_id] etas; // etas[i,1] = log_ec50; etas[i,2] = e0;

  real<lower=0> sigma;
  vector<lower=0>[2] omega_sds; // SDs for IIV
  cholesky_factor_corr[2] L_Omega_corr;
}

transformed parameters {
  real tv_ec50 = exp(tv_log_ec50);
  real gamma = exp(log_gamma);
  vector[n_id] ec50;
  vector[n_id] e0;
  matrix[2,2] Omega_corr;
  cov_matrix[2] Omega;
  vector[N] mu;
  

  {
    vector[2] tmp_eta;
    
    for (i in 1:n_id) {
      tmp_eta = diag_pre_multiply(omega_sds, L_Omega_corr)*etas[,i];
      ec50[i] = exp(tv_log_ec50 + tmp_eta[1]);
      e0[i] = tv_e0 + tmp_eta[2];
    }
  }
  
  Omega_corr = L_Omega_corr * L_Omega_corr';
  Omega = quad_form_diag(Omega_corr, omega_sds);
  
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
  
  L_Omega_corr ~ lkj_corr_cholesky(2);
  omega_sds ~ cauchy(0,2.5);

  for (i in 1:n_id) {
    etas[,i] ~ multi_normal([0.0,0.0],diag_matrix([1.0,1.0]'));
  }

// Likelihhood

  if (prior_only==0) {
    DV ~ normal(mu, sigma);
  }
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
