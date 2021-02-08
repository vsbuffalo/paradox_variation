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
  vector[N] gamma;    // phylogenetic random effects
}

transformed parameters {
  vector[N] mu;
  for (i in 1:N) mu[i] = phi + beta*x[i];
}

model {
  matrix[N, N] L_Sigma;
  phi ~ normal(0, 5);
  beta ~ normal(0, 5);
  sigma_p ~ normal(0, 2);

  L_Sigma = sigma_p * L;
  y ~ multi_normal_cholesky(mu, L_Sigma);
}

generated quantities {
  vector[100] y_rep; 
  vector[100] mu_rep; 
  for (i in 1:100) {
    mu_rep[i] = phi + beta*x_new[i];
  }
}


