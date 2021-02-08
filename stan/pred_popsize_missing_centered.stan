functions {
  vector scale_center(vector x) {
    return (x - mean(x))/sd(x);
  }

  real unscale_slope(real b, vector x, vector y) {
    return b * sd(y) / sd(x);
  }
  real unscale_intercept(real a, real b, vector x, vector y) {
    return sd(y) * (a - b * mean(x) / sd(x)) + mean(y);
  }
}

data {
  // Romiguier et al subset 
  // Complete data on body length -> body mass relationship
  int<lower=1> N_r;
  vector[N_r] body_length_r;
  vector[N_r] body_mass_r;

  // Damuth et al data
  // Complete data on body mass -> pop density relationship
  int<lower=1> N_d;
  vector[N_d] density_d;
  vector[N_d] body_mass_d;

  // main dataset
  // Incomplete data on body length -> pop size (=range * density)
  int<lower=1> N;
  vector[N] body_length;
  vector[N] range_meas;
}

transformed data {
  // The data transformations here are used to put the data
  // back into a linear scale. The data is provided in log10
  // centered data, so prior choice is clearer.
  vector[N_r] sclog10_body_length_r;
  vector[N_d] sclog10_body_mass_d;
  vector[N_r] sclog10_body_mass_r;
  vector[N] sclog10_body_length;
  vector[N_d] sclog10_density_d;

  sclog10_body_length_r = scale_center(log10(body_length_r));
  sclog10_body_mass_d = scale_center(log10(body_mass_d));
  sclog10_body_length = scale_center(log10(body_length));
  sclog10_body_mass_r = scale_center(log10(body_mass_r));
  sclog10_density_d = scale_center(log10(density_d));

}

parameters {
  // body size / mass 
  real alpha_sc;
  real beta_sc;
  real<lower=0> sigma_mass_sc;

  // Damuth
  real gamma_sc; // intercept
  real theta_sc; // slope
  real<lower=0> sigma_density_sc;

  // Mising data
  vector[N] sclog10_body_mass;
  //real<lower=0> log10_body_mass[N];
  vector[N] sclog10_density;
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
  //vector[N] log10_mu_range;

  // Transformed parameters for known data
  for (i in 1:N_r) {
    mu_mass_r[i] = alpha_sc + beta_sc*sclog10_body_length_r[i];
  }
  for (i in 1:N_d) {
    mu_density_d[i] = gamma_sc + theta_sc*sclog10_body_mass_d[i];
  }


  // Transformed parameters for missing data
  for (i in 1:N) {
    mu_unknown_mass[i] = alpha_sc + beta_sc*sclog10_body_length[i];
    mu_unknown_density[i] = gamma_sc + theta_sc*sclog10_body_mass[i];
    //log10_mu_range[i] = lambda + eta*log10(body_mass[i]);
  }

}

// normal version
model {
  // model for body length -->  mass based on Romiguier data
  alpha_sc ~ normal(0, 2);
  beta_sc ~ normal(0, 2);
  sigma_mass_sc ~ cauchy(0, 5);
  sclog10_body_mass_r ~ normal(mu_mass_r, sigma_mass_sc);
 
  // model for body mass --> density based on Damuth data
  gamma_sc ~ normal(0, 2); //
  theta_sc ~ normal(0, 2);
  sigma_density_sc ~ cauchy(0, 5);

  sclog10_density_d ~ normal(mu_density_d, sigma_density_sc);

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
  sclog10_body_mass ~ normal(mu_unknown_mass, sigma_mass_sc);
  sclog10_density ~ normal(mu_unknown_density, sigma_density_sc);
} 

generated quantities {
  // parameters on original scales
  real alpha;
  real beta;
  real gamma;
  real theta;
  real<lower=0> sigma_mass;
  real<lower=0> sigma_density;
  vector[N] log10_body_mass;
  vector[N] log10_density;
  vector[N] log10_body_mass_pred;
  vector[N] log10_density_pred;
  vector[N] log10_popsize;
  real log10_popsize_range;

  alpha = unscale_intercept(alpha_sc, beta_sc, log10(body_length_r), log10(body_mass_r));
  beta = unscale_slope(beta_sc, log10(body_length_r), log10(body_mass_r));
  sigma_mass = sigma_mass_sc * sd(log10(body_mass_r));
  gamma = unscale_intercept(gamma_sc, theta_sc, log10(body_mass_d), log10(density_d));
  theta = unscale_slope(theta_sc, log10(body_mass_d), log10(density_d));
  sigma_density = sigma_density_sc * sd(log10(density_d));

  // missing data transformed to original scales
  log10_body_mass = alpha + beta*log10(body_length);
  log10_density = gamma + theta*log10_body_mass;

  log10_popsize = log10_density + log10(range_meas);
  log10_popsize_range = max(log10_popsize) - min(log10_popsize);

  for (i in 1:N)  {
    log10_body_mass_pred[i] = normal_rng(alpha + beta*log10(body_length[i]), sigma_mass);
    log10_density_pred[i] = normal_rng(gamma + theta*log10_body_mass_pred[i], sigma_density);
  }

}
