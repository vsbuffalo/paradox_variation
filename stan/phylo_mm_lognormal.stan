data {
  int<lower=0> N;
  vector[N] map_length;
  vector[N] log10_popsize;
  corr_matrix[N] cor;  // tree's correlation matrix
  vector[100] x_new;
}

transformed data {
  real shift = 1000;
  vector[N] shift_log10_popsize;
  matrix[N, N] L;
  shift_log10_popsize = log10(shift) + log10_popsize;
  L = cholesky_decompose(cor);
}

parameters {
  real phi;               // intercept: warning, this is shifted!!
  real beta;
  vector[N] log10_N;      // measurement noise
  real<lower=0> sigma_p;  // phylogenetic RE scaling SD
  real<lower=0> sigma_r;  // residual error SD
  vector[N] gamma;    // phylogenetic random effects
  /* vector[N] eta;          // this used to reparameterize the MVN */ 
}

transformed parameters {
  vector[N] l_mu;
  vector[N] kappa;
  vector[N] mu;
  for (i in 1:N) l_mu[i] = phi + beta*log10_N[i];
  kappa = sigma_p * (L*gamma);
  mu = l_mu + kappa;
}

model {
  phi ~ normal(0, 5);
  beta ~ normal(0, 5);
  sigma_r ~ cauchy(0, 5);
  sigma_p ~ cauchy(0, 5);
  log10_N ~ normal(log10_popsize, 1.1);
  gamma ~ std_normal(); // implies: beta ~ multi_normal(mu, Sigma)

  map_length ~ lognormal(mu, sigma_r);
}

generated quantities {
  vector[100] y_rep; 
  vector[100] mu_rep; 
  for (i in 1:100) {
    mu_rep[i] = exp(phi + beta*x_new[i]) - log10(shift);
    y_rep[i] = lognormal_rng(phi + beta*x_new[i], sigma_r) - log10(shift);
  }
}


