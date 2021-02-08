data {
  int<lower=2> N;
  vector[N] Y;
  vector[N] X;
  cov_matrix[N] cov;  // tree's covariance matrix
  int prior_only;  // should the likelihood be ignored?
}

transformed data {
  matrix[N, N] L;
  L = cholesky_decompose(cov);
}


parameters {
  real alpha;             // intercept
  real beta;              // slope 
  real power;             // power on x
  real<lower=0> sigma_p;  // phylogenetic RE scaling SD
  real<lower=0> sigma_r;  // residual error SD
  vector[N] gamma;    // phylogenetic random effects
  vector[N] eta;          // this used to reparameterize the MVN 
}

transformed parameters {
  vector[N] nl_mu;
  vector[N] kappa;
  vector[N] mu;
  for (i in 1:N) nl_mu[i] = alpha + beta*pow(X[i], power);
  kappa = sigma_p * (L*gamma);
  mu = nl_mu + kappa;
}

model {
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  power ~ normal(0, 2);
  sigma_r ~ cauchy(0, 1);
  sigma_p ~ cauchy(0, 1);
  gamma ~ normal(0, 5);
  eta ~ std_normal(); // implies: beta ~ multi_normal(mu, Sigma)
  Y ~ normal(mu, sigma_r);
}

generated quantities {
  real lambda; 
  lambda = pow(sigma_p, 2.) / (pow(sigma_r, 2.) + pow(sigma_p, 2.));
}
