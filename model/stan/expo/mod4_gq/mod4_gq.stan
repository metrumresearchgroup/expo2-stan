
data {
  int<lower=0> N ; // Number of observations
  int<lower=0> K;  // Number of covariates
  int<lower=0> n_id ; // Number of individuals
  vector[N] DV;
  vector[N] Conc;
  matrix[n_id, K] X_e0;  // Matrix of covariates (one row per subject) on E0
  matrix[n_id, K] X_ec50;  // Matrix of covariates (one row per subject) on EC50
  array[N] int ID;
  array[n_id] int<lower=1, upper=N> id_start;
  array[n_id] int<lower=1, upper=N> id_end;
  
  int<lower=0,upper=1> prior_only;
  int n_sim;
}


parameters {
  real<lower=0, upper=100> emax;
  real tv_log_ec50;
  real tv_e0;
  real log_gamma;
  vector[K] beta_e0;
  vector[K] beta_log_ec50;

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
      ec50[i] = exp(tv_log_ec50 + X_ec50[i,] * beta_log_ec50 + tmp_eta[1]);
      e0[i] = tv_e0 + X_e0[i,] * beta_e0 + tmp_eta[2];
    }
  }
  
  Omega_corr = L_Omega_corr * L_Omega_corr';
  Omega = quad_form_diag(Omega_corr, omega_sds);
  
  for (i in 1:N) {
    mu[i] = e0[ID[i]] + (emax-e0[ID[i]])*pow(Conc[i]/100,gamma) / (pow(ec50[ID[i]]/100,gamma) + pow(Conc[i]/100,gamma));
  }
}


generated quantities {
  vector[n_id] log_like_approx;
  
      {
        // temporary variables
          matrix[n_id,n_sim] sim_e0;
          matrix[n_id, n_sim] sim_ec50;
          
          matrix[n_id,n_sim] sim_penalty;
          vector[n_sim] sum_ll;
          matrix[N,n_sim] mutemp;
          matrix[N,n_sim] lltemp;
          // vector[n_id] sum_ll;
          vector[2] tmp_eta;
          
          for (i in 1:n_id) {
            for (j in 1:n_sim) {
              tmp_eta = multi_normal_rng([0,0], diag_matrix([1.0,1.0]'));
              
              tmp_eta = diag_pre_multiply(omega_sds, L_Omega_corr)*tmp_eta;
              sim_ec50[i,j] = exp(tv_log_ec50 + X_ec50[i,] * beta_log_ec50 + tmp_eta[1]);
              sim_e0[i,j] = tv_e0 + X_e0[i,] * beta_e0 + tmp_eta[2];
              
              sim_penalty[i,j] = multi_normal_lpdf(tmp_eta | [0,0],Omega);
            }
          }
          
          
          for (i in 1:N) {
            for (j in 1:n_sim) {
              mutemp[i,j] = sim_e0[ID[i],j] + (emax-sim_e0[ID[i],j])*pow(Conc[i]/100,gamma) / (pow(sim_ec50[ID[i],j]/100,gamma) + pow(Conc[i]/100,gamma));
              lltemp[i,j] = normal_lpdf(DV[i] | mutemp[i,j], sigma) ;
            }
          }
          
          for (i in 1:n_id) {
            for (j in 1:n_sim) {
              sum_ll[j] = sum(lltemp[id_start[i]:id_end[i],j]) + sim_penalty[i,j];
            }
            log_like_approx[i] = log_sum_exp(sum_ll) - log(n_sim);
          }
      }
      
}
