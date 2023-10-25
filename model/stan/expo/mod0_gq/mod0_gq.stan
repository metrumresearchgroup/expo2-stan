
data {
  int<lower=0> N ;
  int<lower=0> n_id ;
  vector[N] DV;
  vector[N] Conc;
  array[N] int ID;
  array[n_id] int start;
  array[n_id] int end;
  
  int n_sim; // Number of simulations used to approximate log-likelihood
}

parameters {
  real<lower=0, upper=100> emax;
  real tv_log_ec50;
  real tv_e0;
  real log_gamma;
  
  vector[n_id] eta_log_ec50;
  vector[n_id] eta_e0;
  
  real<lower=0> sigma;
  real<lower=0> omega_e0;
  real<lower=0> omega_log_ec50;
}

transformed parameters {
  real tv_ec50 = exp(tv_log_ec50);
  real gamma = exp(log_gamma);
  vector[n_id] ec50;
  vector[n_id] e0;
  vector[N] mu;
  
  ec50 = exp(tv_log_ec50 + omega_log_ec50 * eta_log_ec50);
  e0 = tv_e0 + omega_e0 * eta_e0;
  
  for (i in 1:N) {
    mu[i] = e0[ID[i]] + (emax-e0[ID[i]])*pow(Conc[i]/100,gamma) / (pow(ec50[ID[i]]/100,gamma) + pow(Conc[i]/100,gamma));
  }
  
}


generated quantities {

  vector[n_id] log_like_approx;
  vector[n_id] mean_penalty;
  vector[n_id] mean_ll;

      {
        // temporary variables
        matrix[n_id,n_sim] sim_e0;
        matrix[n_id, n_sim] sim_ec50;
        
        matrix[n_id,n_sim] sim_penalty;
        vector[n_sim] sum_ll;
        matrix[N,n_sim] mutemp;
        matrix[N,n_sim] lltemp;
        real temp_eta_e0;
        real temp_eta_log_ec50;
        
        for (i in 1:n_id) {
          for (j in 1:n_sim) {
            temp_eta_e0 = normal_rng(0,1);
            temp_eta_log_ec50 = normal_rng(0,1);
            sim_ec50[i,j] = exp(tv_log_ec50  + omega_log_ec50 * temp_eta_log_ec50);
            sim_e0[i,j] = tv_e0 + omega_e0 * temp_eta_e0;
            
           sim_penalty[i,j] = std_normal_lpdf(temp_eta_log_ec50) +
                              std_normal_lpdf(temp_eta_e0);
            // sim_penalty[i,j] = multi_normal_lpdf([temp_eta_log_ec50,temp_eta_e0] | [0,0], diag_matrix([pow(omega_log_ec50,2),pow(omega_e0,2)]'));
          }
          mean_penalty[i] = mean(sim_penalty[i,]);
        }
        
        
        for (i in 1:N) {
          for (j in 1:n_sim) {
            mutemp[i,j] = sim_e0[ID[i],j] + (emax-sim_e0[ID[i],j])*pow(Conc[i]/100,gamma) / (pow(sim_ec50[ID[i],j]/100,gamma) + pow(Conc[i]/100,gamma));
            lltemp[i,j] = normal_lpdf(DV[i] | mutemp[i,j], sigma) ;
          }
        }
        
        for (i in 1:n_id) {
          for (j in 1:n_sim) {
            sum_ll[j] = sum(lltemp[start[i]:end[i],j]) + sim_penalty[i,j];
          }
          mean_ll[i] = mean(lltemp[start[i]:end[i],]);
          log_like_approx[i] = log_sum_exp(sum_ll) - log(n_sim);
        }
      }
    
}
