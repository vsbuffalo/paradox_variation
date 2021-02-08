data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  corr_matrix[N] cor;  // tree's correlation matrix
  vector[100] x_new;
}

transformed data {
  matrix[N, N] L;
  L = cholesky_decompose(cor);
}

parameters {
  real phi;
  real beta;
  real<lower=0> sigma_p;  // phylogenetic RE scaling SD
  real<lower=0> sigma_r;  // residual error SD
  vector[N] gamma;    // phylogenetic random effects
}

transformed parameters {
  vector[N] l_mu;
  vector[N] kappa;
  vector[N] mu;
  for (i in 1:N) l_mu[i] = phi + beta*x[i];
  kappa = sigma_p * (L*gamma);
  mu = l_mu + kappa;
}

model {
  phi ~ normal(0, 2);
  beta ~ normal(0, 2);
  sigma_r ~ normal(0, 2);
  sigma_p ~ normal(0, 2);
  gamma ~ std_normal(); // implies: beta ~ multi_normal(mu, Sigma)
  y ~ normal(mu, sigma_r);
}

generated quantities {
  vector[100] y_rep; 
  vector[100] mu_rep; 
  for (i in 1:100) {
    mu_rep[i] = phi + beta*x_new[i];
    y_rep[i] = normal_rng(phi + beta*x_new[i], sigma_r);
  }
}


