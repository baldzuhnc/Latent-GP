
data {
  int<lower=2> K;                 // Number of categories
  int<lower=1> N;                 // Total number of observations
  int<lower=1> P;                 // Number of items
  int<lower=1> tmax;              // Max observed time point
  int<lower=1> I;                 // Number of participants
  real<lower=0> length_scale; // fixed globally
  array[N] int<lower=1, upper=I> participant; // Participant ID for each observation
  array[N] int<lower=1, upper=P> item;        // Item ID for each observation
  array[N] int<lower=1, upper=tmax> time;     // Observed time point for each observation
  array[N] int<lower=1, upper=K> Y;           // Response variable
}

transformed data {
  array[tmax] real time_gp = to_array_1d(linspaced_vector(tmax, 1.0, tmax * 1.0));
  real delta = 1e-5; // jiggle to main diagonal
}
 
parameters {

  // Non-centered parameters for GP hyperparameters, per participant
  vector[I] z_log_gp_amplitude;                   

  // Population-level standard deviation for GP hyperparameters (on log scale for positive params)
  // real mu_log_gp_amplitude;  // REMOVED FOR IDENTIFICATION (fixes mean log-amp to 0)
  real<lower=0> sigma_log_gp_amplitude;
  
  // Non-centered GP latent variables
  matrix[I, tmax] z_theta; 

  // IRT parameters
  //vector[P] item_disc;
  array[P] ordered[K - 1] item_thresh;           // Category thresholds per item
  // real mean_thresh;                           // REMOVED FOR IDENTIFICATION (fixes mean thresh to 0)
  real<lower=0> sd_thresh;                       // Prior SD for category thresholds
}

transformed parameters {

  // mu_log_gp_amplitude removed from this calculation
  vector[I] log_gp_amplitude_ind = sigma_log_gp_amplitude * z_log_gp_amplitude; 
  vector[I] gp_amplitude_ind = exp(log_gp_amplitude_ind);                                     

  matrix[I, tmax] theta;

  for (i in 1:I) {
    matrix[tmax, tmax] K_gp = gp_exp_quad_cov(time_gp, gp_amplitude_ind[i], length_scale);
    
    for (t in 1:tmax)
      K_gp[t, t] += delta;

    matrix[tmax, tmax] L = cholesky_decompose(K_gp);
    theta[i] = (L * z_theta[i]')';
  }
}

model {
  // IRT priors
  //item_disc ~ normal(1, 0.5);
  // mean_thresh ~ normal(0, 5); // REMOVED
  sd_thresh ~ normal(0, 3);
  for (j in 1:P) {
    item_thresh[j] ~ normal(0, 1); // mean_thresh changed to 0
  }

  // mu_log_gp_amplitude ~ normal(-0.25, 0.5); // REMOVED
  sigma_log_gp_amplitude ~ normal(-0.25, 0.5); // halfnoral(0,1), amplitude centered around 0 due to log(0)

  // Non-centered hierarchical priors
  z_log_gp_amplitude ~ normal(0, 1);           
  to_vector(z_theta) ~ normal(0, 1); 

  // Likelihood
  for (n in 1:N) {
    int pid = participant[n];
    int it = item[n];
    int t = time[n];
    //real eta = item_disc[it] * theta[pid, t];
    real eta = theta[pid, t];
    Y[n] ~ ordered_logistic(eta, item_thresh[it]);
  }
}

generated quantities {
  vector[N] log_lik;                 // Log-likelihood for each observation
  array[I, P, tmax] int<lower=1, upper=K> Y_rep; // Posterior predictive responses for observed time points

  // Compute log-likelihood
  for (n in 1:N) {
    int pid = participant[n];
    int it = item[n];
    int t = time[n];
    //real eta = item_disc[it] * theta[pid, t];
    real eta = theta[pid, t];
    log_lik[n] = ordered_logistic_lpmf(Y[n] | eta, item_thresh[it]);
  }

  // Generate posterior predictive responses for the observed time points
  for (i in 1:I) {
    for (j in 1:P) {
      for (t in 1:tmax) {
        //real eta = item_disc[j] * theta[i, t];
        real eta = theta[i, t];
        Y_rep[i, j, t] = ordered_logistic_rng(eta, item_thresh[j]);
      }
    }
  }
}
