data {
  // Romiguier et al subset
  int<lower=1> N_subset;
  real<lower=0> body_length_subset[N_subset];
  real<lower=0> body_mass_subset[N_subset];

  // Damuth et al data
  int<lower=1> N_damuth;
  real<lower=0> density_damuth[N_damuth];
  real<lower=0> body_mass_damuth[N_damuth];

  // main dataset
  int<lower=1> N;
  real<lower=0> body_lengths[N];
  real<lower=0> ranges[N];
}

transformed data {
  real log10_body_lengths_subset[N_subset];
  real log10_body_mass_subset[N_subset];
  real log10_body_lengths[N];
  real log10_density_damuth[N_damuth];
  real log10_mass_damuth[N_damuth];

  // Romiguier et al subset
  log10_body_lengths_subset = log10(body_length_subset);
  log10_body_mass_subset = log10(body_mass_subset);

  // main dataset
  log10_body_lengths = log10(body_lengths);

  // Damuth data
  log10_density_damuth = log10(density_damuth);
  log10_mass_damuth = log10(body_mass_damuth);


}

parameters {
  // body size / mass 
  real beta_0;
  real beta_size;
  real<lower=0> sigma_mass;

  // Damuth
  real gamma; // intercept
  real theta; // slope
  real<lower=0> sigma_mass_damuth;
}

transformed parameters {
  vector[N_subset] mu_mass;
  vector[N_damuth] mu_density_damuth;
  for (i in 1:N_subset) {
    mu_mass[i] = beta_0 + beta_size*log10_body_lengths_subset[i];
  }
  for (i in 1:N_damuth) {
    mu_density_damuth[i] = gamma + theta*log10_mass_damuth[i];
  }
}

// normal version
model {
  // model for body length -->  mass  
  beta_0 ~ student_t(3, 0, 3);
  beta_size ~ student_t(3, 0, 3);
  sigma_mass ~ student_t(3, 0, 3);
  body_mass_subset ~ lognormal(mu_mass, sigma_mass);
 
  // Damuth model
  gamma ~ student_t(3, 0, 3); //
  theta ~ student_t(3, 0, 3);
  sigma_mass_damuth ~ student_t(3, 0, 3);
  density_damuth ~ lognormal(mu_density_damuth, sigma_mass_damuth);

} 

generated quantities {
  vector[N] log10_pred_popsize; 
  vector[N] log10_pred_body_mass; 
  for (i in 1:N)  {
    log10_pred_body_mass[i] = normal_rng((beta_0 + beta_size*log10_body_lengths[i])/log(10), 
                                         sigma_mass/log(10));
    log10_pred_popsize[i] = log10(ranges[i]*pow(10, 
                                  normal_rng((gamma + theta*log10_pred_body_mass[i])/log(10), 
                                              sigma_mass_damuth/log(10))));
  }
}

/* model { */
/*   // model for body length -->  mass */  
/*   beta_0 ~ student_t(3, 0, 3); */
/*   beta_size ~ student_t(3, 0, 3); */
/*   sigma_mass ~ student_t(3, 0, 3); */
/*   body_mass_subset ~ lognormal(mu_mass, sigma_mass); */
 
/*   // Damuth model */
/*   gamma ~ student_t(3, 0, 3); */
/*   theta ~ student_t(3, 0, 3); */
/*   sigma_mass_damuth ~ student_t(3, 0, 3); */
/*   density_damuth ~ lognormal(mu_density_damuth, sigma_mass_damuth); */

/* } */ 

/* generated quantities { */
/*   vector[N] pred_popsize; */ 
/*   vector[N] pred_body_mass; */ 
/*   for (i in 1:N)  { */
/*     pred_body_mass[i] = lognormal_rng(beta_0 + beta_size*log10_body_lengths[i], */ 
/*                                       sigma_mass); */
/*     pred_popsize[i] = ranges[i]*lognormal_rng(gamma + theta*log10(pred_body_mass[i]), */ 
/*                                               sigma_mass_damuth); */
/*   } */
/* } */
