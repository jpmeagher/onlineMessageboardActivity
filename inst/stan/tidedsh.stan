//
// This Stan Model defines the Time-Dependant Doubly Stochastic
// Hawkes process Model we use to Analyse Online-Message-board Activity.
//

data {
  int<lower=0> N; // total arrivals
  int N_immigrant; // arrivals of generation 0
  int N_offspring; // offspring arrivals
  real a; // observed interval lower bound
  real<lower=a> b; // observed interval upper boiund
  vector<lower=a, upper=b>[N] t; // arrival times
  int<lower=0, upper=N> v[N]; // branching structure
  int<lower=0, upper=1> is_ds; // are there latent reproduction numbers
  // int<lower=0, upper=1> is_td; // is there time-dependance
  int K; // number of sinusoidal basis
  real<lower=0> f[K]; // basis frequencies
}
transformed data {
  int ind_immigrant; // generation 0 index
  vector[N_immigrant] t_immigrant; // generation 0 arrival times
  int ind_offspring; // offspring index
  vector[N_offspring] t_offspring; // offspring arrival times
  int v_offspring[N_offspring]; // parent indices
  vector<lower=0, upper=N>[N] z; // offspring
  int N_latent; // offspring arrivals
  vector<lower=0>[K] w; // angular frequency
  matrix[N, 2*K] sinusoidal_basis; // sinusoidal basis
  matrix[N_immigrant, 2*K] sinusoidal_basis_immigrant;
  matrix[N_offspring, 2*K] sinusoidal_basis_offspring;
  vector[2*K] sinusoidal_basis_integral; // integral basis

  // separate out immigrants and offspring
  ind_immigrant = 0;
  ind_offspring = 0;
  for (i in 1:N) {
    if (v[i] == 0) {
      ind_immigrant += 1;
      t_immigrant[ind_immigrant] = t[i];
    } else {
      ind_offspring += 1;
      t_offspring[ind_offspring] = t[i];
      v_offspring[ind_offspring] = v[i];
    }
  }
  for (i in 1:N) {
    z[i] = 0;
    for(j in (i+1):N) {
      z[i] = z[i] + (v[j] == i);
    }
  }
  if (is_ds) {
    N_latent = N;
  } else {
    N_latent = 0;
  }
  // // set the sinusoidal basis
  for (k in 1:K) {
    w[k] = 2 * pi() * f[k];
    sinusoidal_basis[, (2*k) - 1] = sin(w[k] * t);
    sinusoidal_basis[, 2*k] = cos(w[k] * t);
    sinusoidal_basis_immigrant[, (2*k) - 1] = sin(w[k] * t_immigrant);
    sinusoidal_basis_immigrant[, 2*k] = cos(w[k] * t_immigrant);
    sinusoidal_basis_offspring[, (2*k) - 1] = sin(w[k] * t_offspring);
    sinusoidal_basis_offspring[, 2*k] = cos(w[k] * t_offspring);
    sinusoidal_basis_integral[(2*k) - 1] = (cos(a * w[k]) - cos(b * w[k])) / w[k];
    sinusoidal_basis_integral[2*k] = (sin(b * w[k]) - sin(a * w[k])) / w[k];
  }
}

parameters {
  real<lower=0> mu; // immigrant intensity
  real log_alpha[is_ds];
  // real<lower=0> sigma_log_alpha[is_ds];
  real<lower=0> R;
  // real log_R;
  // real<lower=0> sigma_log_R;
  vector<lower=0>[N_latent] nu; // stochastic branching numbers
  real<lower=0> eta; // exponential rate
  vector<lower=-1, upper=1>[2*K] beta;
}
transformed parameters {
  // real<lower=0> R;
  real<lower=0> alpha[is_ds];
  vector<lower=0>[N] rho; // relative activity levels
  vector<lower=0>[N_immigrant] rho_immigrant; // relative activity levels
  vector<lower=0>[N_offspring] rho_offspring; // relative activity levels
  matrix[N, 2*K] exp_sinusoidal_basis_integral; // compute the integral for the exponent times the basis
  real immigrant_log_intensity; // immigrant arrival intensity
  real<lower=0> immigrant_scalar; // immigrant likelihood scalar
  vector[N_offspring] tmp_offspring_log_intensity_vector; // temporary vector for offspring intensity
  vector<lower=0>[N] tmp_offspring_scalar_vector; // temporary vector for offspring likelihood scalar
  real offspring_log_intensity; // offspring arrival intensity
  real<lower=0> offspring_scalar; // offspring likelihood scalar

  // R = exp(log_R);
  if (is_ds) {
    alpha = exp(log_alpha);
  }
  rho = rep_vector(1, N);
  rho_immigrant = rep_vector(1, N_immigrant);
  rho_offspring = rep_vector(1, N_offspring);
  if (K != 0) {
    rho += sinusoidal_basis * beta;
    rho_immigrant += sinusoidal_basis_immigrant * beta;
    rho_offspring += sinusoidal_basis_offspring * beta;
  }
  for (k in 1:K) {
    exp_sinusoidal_basis_integral[, 2*k - 1] = ( eta / (eta^2 + w[k]^2) ) * (
      eta * (sin(w[k] * t) - (sin(w[k] * b) * exp(- eta * (b - t)))) +
      w[k] * (cos(w[k] * t) - (cos(w[k] * b) * exp(- eta * (b - t))))
    );
    exp_sinusoidal_basis_integral[, 2*k] = ( eta / (eta^2 + w[k]^2) ) * (
      w[k] * ((sin(w[k] * b) * exp(- eta * (b - t)) - sin(w[k] * t))) -
      eta * ((cos(w[k] * b) * exp(- eta * (b - t)) - cos(w[k] * t)))
    );
  }
  // immigrant contribution to the likelihood
  // intensity
  immigrant_log_intensity = (N_immigrant * log(mu)) + sum(log(rho_immigrant));
  // scalar
  immigrant_scalar = (b - a);
  if (K != 0){
    immigrant_scalar += sinusoidal_basis_integral' * beta;
  }
  immigrant_scalar = mu * immigrant_scalar;
  //  offspring contribution
  // intensity
  // number of offspring
  if (is_ds) {
    offspring_log_intensity = sum(z .* log(nu));
  } else {
    offspring_log_intensity = N_offspring * log(R);
  }
  // // relative activity
  offspring_log_intensity += sum(log(rho_offspring));
  // generation interval
  for (i in 1:N_offspring) {
    tmp_offspring_log_intensity_vector[i] = exponential_lpdf(t_offspring[i] - t[v_offspring[i]] | eta);
  }
  offspring_log_intensity += sum(tmp_offspring_log_intensity_vector);
  // scalar
  for (i in 1:N) {
    tmp_offspring_scalar_vector[i] = exponential_cdf(b - t[i], eta);
  }
  if (K != 0) {
    tmp_offspring_scalar_vector += exp_sinusoidal_basis_integral * beta;
  }
   if (is_ds) {
    offspring_scalar = sum(nu .* tmp_offspring_scalar_vector);
  } else {
    offspring_scalar = R * sum(tmp_offspring_scalar_vector);
  }
  // tmp_offspring_scalar_vector = rep_vector(1, N);
  // if (K != 0) {
  //   tmp_offspring_scalar_vector = 1 + exp_sinusoidal_basis_integral * beta;
  // }
  //
  // offspring_scalar = eta * sum(tmp_offspring_scalar_vector);
}

model {
  target += immigrant_log_intensity;
  target += - immigrant_scalar;
  target += offspring_log_intensity;
  target += - offspring_scalar;

  // target += normal_lpdf(log_R | 0, sigma_log_R);
  // target += cauchy_lpdf(sigma_log_R | 0, 1);

  if (is_ds) {
    target += gamma_lpdf(nu | alpha[is_ds], alpha[is_ds] / R);
    target += normal_lpdf(log_alpha | 0, 1);
    // target += cauchy_lpdf(sigma_log_alpha | 0, 1);
  }
  // if (is_ds) {
  //   target += gamma_lpdf(nu | alpha[is_ds], alpha[is_ds] / R);
  //   target += normal_lpdf(log_alpha | 0, sigma_log_alpha);
  //   target += cauchy_lpdf(sigma_log_alpha | 0, 1);
  // }

  target += uniform_lpdf(beta | -1, 1);

  // if (is_td) {
  //   // target += normal_lpdf(beta | 0, delta_beta[is_td]);
  //   // target += inv_gamma_lpdf(delta_beta | 2, 1);
  //   target += std_normal_lpdf(beta_raw);
  //   target += cauchy_lpdf(sigma_beta | 0, 1);
  // }


  // for(k in 1:K) {
  //   target += std_normal_lpdf(beta_raw[, k]);
  //   target += cauchy_lpdf(sigma_beta[k] | 0, 1);
  // }
}
generated quantities {
  vector[K] A;
  vector[K] phase;

  // for (k in 1:K) {
  //   A[k] = sqrt(dot_self(beta[k]));
  //   phase[k] = atan2(beta[k][2], beta[k][1]);
  // }
  for (k in 1:K) {
    A[k] = sqrt(beta[(2*k) - 1]^2 + beta[2*k]^2);
    phase[k] = atan2(beta[2*k], beta[(2*k) - 1]);
  }
}
