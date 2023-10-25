
data {
  int<lower=0> N ;
  int<lower=0> n_id ;
  vector[N] DV;
  vector[N] Conc;
  array[N] int ID;
  int<lower=0,upper=1> prior_only;
  array[n_id] int<lower=1, upper=N> id_start;
  array[n_id] int<lower=1, upper=N> id_end;
}

transformed data {
  array[N] int BLQ;
  int n_sim = 20;
  
  for (i in 1:N) {
    BLQ[i] = (DV[i] > 99.9) ? 1 : 0;
  }
}

parameters {
  real<lower=0, upper=100> emax;
  real tv_log_ec50;
  real tv_e0;
  real log_gamma;
  
  array[n_id] vector[2] etas; // etas[i,1] = log_ec50; etas[i,2] = e0;

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
  
  {
    vector[2] tmp_eta;
    
    for (i in 1:n_id) {
      tmp_eta = diag_pre_multiply(omega_sds, L_Omega_corr)*etas[i];
      ec50[i] = exp(tv_log_ec50 + tmp_eta[1]);
      e0[i] = tv_e0 + tmp_eta[2];
    }
  }
  
  Omega_corr = L_Omega_corr * L_Omega_corr';
  Omega = quad_form_diag(Omega_corr, omega_sds);
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

  etas ~ multi_normal([0.0,0.0],diag_matrix([1.0,1.0]'));

// Likelihhood

  vector[N] mu;
  
  for (i in 1:N) {
    mu[i] = e0[ID[i]] + (emax-e0[ID[i]])*pow(Conc[i]/100,gamma) / (pow(ec50[ID[i]]/100,gamma) + pow(Conc[i]/100,gamma));
  }
  
  if (prior_only==0) {
    for (i in 1:N) {
      if (BLQ[i]==1) { 
        target += log1m(Phi((99.9 - mu[i])/sigma));
      } else {
        target += normal_lpdf(DV[i] | mu[i], sigma);
      }
    }
  }
}



generated quantities{
   vector[n_id] sim_ec50;
   vector[n_id] sim_e0;
   vector[N] sim_ipred;
   vector[N] yhat;
   vector[N] sim_new_dv;
   vector[N] sim_cur_dv;
   matrix[n_id,n_sim] pop_sim_e0;
   matrix[n_id, n_sim] pop_sim_ec50;
   vector[n_id] log_like_approx;
   
  
   vector[n_id] log_like;

   for (i in 1:N) {
     yhat[i] = e0[ID[i]] + (emax-e0[ID[i]])*pow(Conc[i]/100,gamma) / (pow(ec50[ID[i]]/100,gamma) + pow(Conc[i]/100,gamma));
   }
   
   for (i in 1:n_id) {
     log_like[i] = 0;
     for (j in id_start[i]:id_end[i]) {
       if (BLQ[i]==1) { 
         log_like[i] += log1m(Phi((99.9 - yhat[i])/sigma));
       } else {
         log_like[i] += normal_lpdf(DV[i] | yhat[i], sigma);
       }
     }
     log_like[i] += multi_normal_lpdf(etas[i] | [0,0],diag_matrix([1.0,1.0]') );
   }
   
   
   {
    vector[2] tmp_eta;
    
    for (i in 1:n_id) {
      tmp_eta = multi_normal_rng([0,0], diag_matrix([1.0,1.0]'));
      
      tmp_eta = diag_pre_multiply(omega_sds, L_Omega_corr)*tmp_eta;
      sim_ec50[i] = exp(tv_log_ec50 + tmp_eta[1]);
      sim_e0[i] = tv_e0 + tmp_eta[2];
    }
  }

  for (i in 1:N) {
    sim_ipred[i] = sim_e0[ID[i]] + (emax-sim_e0[ID[i]])*pow(Conc[i]/100,gamma) / (pow(sim_ec50[ID[i]]/100,gamma) + pow(Conc[i]/100,gamma));
    sim_new_dv[i] = normal_rng(sim_ipred[i], sigma);
    sim_cur_dv[i] = normal_rng(yhat[i], sigma);
    // Censor simulated data following the DGP modeled above
    sim_new_dv[i] = (sim_new_dv[i] > 99.9) ? 100.0 : sim_new_dv[i];
    sim_cur_dv[i] = (sim_cur_dv[i] > 99.9) ? 100.0 : sim_cur_dv[i];
  }

      {
        // temporary variables
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
            pop_sim_ec50[i,j] = exp(tv_log_ec50 + tmp_eta[1]);
            pop_sim_e0[i,j] = tv_e0 + tmp_eta[2];
            
            sim_penalty[i,j] = multi_normal_lpdf([log(pop_sim_ec50[i,j]),pop_sim_e0[i,j]] | [tv_log_ec50,tv_e0],Omega);
          }
        }
        
        
        for (i in 1:N) {
          for (j in 1:n_sim) {
            mutemp[i,j] = pop_sim_e0[ID[i],j] + (emax-pop_sim_e0[ID[i],j])*pow(Conc[i]/100,gamma) / (pow(pop_sim_ec50[ID[i],j]/100,gamma) + pow(Conc[i]/100,gamma));
            
            if (BLQ[i]==1) { 
              lltemp[i,j] = log1m(Phi((99.9 - mutemp[i,j])/sigma));
            } else {
              lltemp[i,j] = normal_lpdf(DV[i] | mutemp[i,j], sigma);
            }
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


