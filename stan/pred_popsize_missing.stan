data {
  // Romiguier et al subset 
  // Complete data on body length -> body mass relationship
  int<lower=1> N_r;
  real<lower=0> body_length_r[N_r];
  real<lower=0> body_mass_r[N_r];

  // Damuth et al data
  // Complete data on body mass -> pop density relationship
  int<lower=1> N_d;
  real<lower=0> density_d[N_d];
  real<lower=0> body_mass_d[N_d];

  // main dataset
  // Incomplete data on body length -> pop size (=range * density)
  int<lower=1> N;
  real<lower=0> body_length[N];
  real<lower=0> range_meas[N];
}

transformed data {
  real log10_body_length_r[N_r];
  real log10_body_mass_d[N_d];
  real log10_body_length[N];
  real log10_range_meas[N];

  // Romiguier et al data
  log10_body_length_r = log10(body_length_r);

  // Damuth data
  log10_body_mass_d = log10(body_mass_d);

  // main dataset
  log10_body_length = log10(body_length);
  log10_range_meas = log10(range_meas);
}

parameters {
  // body size / mass 
  real alpha;
  real beta;
  real<lower=0> sigma_mass;

  // Damuth
  real gamma; // intercept
  real theta; // slope
  real<lower=0> sigma_density;

  // Mising data
  real<lower=0> body_mass[N];
  //real<lower=0> log10_body_mass[N];
  real<lower=0> density[N];
  //real<lower=0> log10_density[N];

  // Ranges -- see note in model section
  //real<lower=0> log10_range[N];
  //real<lower=0> sigma_range;
  //real<lower=0> tau;
  //real lambda; // intercept in range ~ body_mass
  //real eta; // slope in range ~ body_mass
}

transformed parameters {
  vector[N_r] mu_mass_r;
  vector[N_d] mu_density_d;
  vector[N] mu_unknown_density;
  vector[N] mu_unknown_mass;
  vector[N] log10_mu_range;

  // Transformed parameters for known data
  for (i in 1:N_r) {
    mu_mass_r[i] = alpha + beta*log10_body_length_r[i];
  }
  for (i in 1:N_d) {
    mu_density_d[i] = gamma + theta*log10_body_mass_d[i];
  }


  // Transformed parameters for missing data
  for (i in 1:N) {
    mu_unknown_mass[i] = alpha + beta*log10(body_length[i]);
    mu_unknown_density[i] = gamma + theta*log10(body_mass[i]);
    //log10_mu_range[i] = lambda + eta*log10(body_mass[i]);
  }

}

// normal version
model {
  // model for body length -->  mass based on Romiguier data
  beta ~ normal(0, 3);
  alpha ~ normal(0, 3);
  sigma_mass ~ normal(0, 3);
  body_mass_r ~ lognormal(mu_mass_r, sigma_mass);
 
  // model for body mass --> density based on Damuth data
  gamma ~ normal(0, 3); //
  theta ~ normal(0, 3);
  sigma_density ~ normal(0, 3);
  density_d ~ lognormal(mu_density_d, sigma_density);

  // range model
  // This regularization approach doesn't match the error in the 
  // range inference close enough, and over regularlizes, leading
  // to smaller ranged organisms having inflated ranges that are unrealistic
  // params for range ~ size model
  //lambda ~ normal(6, 1);  // based on scale/location of parameter
  //tau ~ normal(0, 1); // measurement error; uninformative prior
  //sigma_range ~ normal(0, 3); // variance in range ~ body_mass
  //eta ~ normal(0, 3);
  //log10_range ~ normal(log10_mu_range, sigma_range);
  //log10_range_meas ~ normal(log10_range, tau);

  // full dataset to be augmented/imputed
  body_mass ~ lognormal(mu_unknown_mass, sigma_mass);
  density ~ lognormal(mu_unknown_density, sigma_density);

} 

generated quantities {
  vector[N] log10_body_mass;
  vector[N] log10_density;
  vector[N] log10_popsize;
  for (i in 1:N) {
    log10_body_mass[i] = log10(body_mass[i]);
    log10_density[i] = log10(density[i]);
    log10_popsize[i] = log10(pow(10., log10(range_meas)[i]) * density[i]);
  }
}
