//
// A Poisson-mixture cluster process
//

data {
  int<lower=0> N; // total events
  int<lower=0, upper=N> N_1; // offspring events
  real a; // observed interval lower bound
  real<lower=a> b; // observed interval upper bound
  vector<lower=a, upper=b>[N] t; // arrival times
  int<lower=0, upper=N> bs[N]; // branching structure
  int<lower=1> N_type; // types of event
  int<lower=1, upper=N_type> type[N]; // branching structure
  int<lower=0, upper=1> is_hetero; // is there heterogeneous reproduction
  int K; // number of sinusoidal basis
  real<lower=0> f[K]; // basis frequencies
  real<lower=0> sigma_R; // log R prior standard deviation
  real<lower=0> sigma_phi; // log phi prior standard deviation
  real<lower=0> sigma_eta; // log eta prior standard deviation
}

transformed data {
  int ind_1; // offspring index
  vector[N_1] t_1; // offspring event times
  int bs_1[N_1]; // offspring branching structure
  vector<lower=0, upper=N>[N] z; // offspring count
  int N_latent; // latent reproduction
  vector<lower=0>[K] w; // angular frequency
  matrix[N_1, 2*K] sinusoidal_basis_1; // sinusoidal basis
  // separate out offspring
  ind_1 = 0;
  for (i in 1:N) {
    if (bs[i] != 0) {
      ind_1 += 1;
      t_1[ind_1] = t[i];
      bs_1[ind_1] = bs[i];
    }
    z[i] = 0;
    for(j in (i+1):N) {
      z[i] = z[i] + (bs[j] == i);
    }
  }
  // include latent variables
  if (is_hetero) {
    N_latent = N;
  } else {
    N_latent = 0;
  }
  // set the sinusoidal basis
  for (k in 1:K) {
    w[k] = 2 * pi() * f[k];
    sinusoidal_basis_1[, (2*k) - 1] = sin(w[k] * t_1);
    sinusoidal_basis_1[, 2*k] = cos(w[k] * t_1);
  }
}

parameters {
  real log_R[N_type]; // log reproduction number
  real log_phi[is_hetero * N_type]; // log dispersion parameter
  real log_eta[N_type]; // log gi rate
  vector<lower=0>[N_latent] nu; // heterogeneous reproduction numbers
  vector<lower=-1, upper=1>[2*K] beta; // sinusoid coefficients
}

transformed parameters {
  real<lower=0> R[N_type]; // reproduction number
  real<lower=0> phi[is_hetero * N_type]; // dispersion parameter
  real<lower=0> eta[N_type]; // gi rate
  vector<lower=0>[N_1] alpha_1; // relative activity levels
  matrix[N, 2*K] exp_sinusoidal_basis_integral; // compute the integral for the exponent times the basis
  vector[N_1] tmp_offspring_log_intensity_vector; // temporary vector for offspring intensity
  vector<lower=0>[N] tmp_offspring_scalar_vector; // temporary vector for offspring likelihood scalar
  real offspring_log_intensity; // offspring arrival intensity
  real<lower=0> offspring_scalar; // offspring likelihood scalar

  R = exp(log_R);
  if (is_hetero) {
    phi = exp(log_phi);
  }
  eta = exp(log_eta);
  // sinusoids
  alpha_1 = rep_vector(1, N_1);
  if (K != 0) {
    alpha_1 += sinusoidal_basis_1 * beta;
  }
  for (i in 1:N) {
    for (k in 1:K) {
      exp_sinusoidal_basis_integral[i, 2*k - 1] = ( eta[type[i]] / (eta[type[i]]^2 + w[k]^2) ) * (
        eta[type[i]] * (sin(w[k] * t[i]) - (sin(w[k] * b) * exp(- eta[type[i]] * (b - t[i])))) +
        w[k] * (cos(w[k] * t[i]) - (cos(w[k] * b) * exp(- eta[type[i]] * (b - t[i]))))
      );
      exp_sinusoidal_basis_integral[i, 2*k] = ( eta[type[i]] / (eta[type[i]]^2 + w[k]^2) ) * (
        w[k] * ((sin(w[k] * b) * exp(- eta[type[i]] * (b - t[i])) - sin(w[k] * t[i]))) -
        eta[type[i]] * ((cos(w[k] * b) * exp(- eta[type[i]] * (b - t[i])) - cos(w[k] * t[i])))
      );
    }
  }
  // offspring intensity
  // time-dependent contribution
  offspring_log_intensity = sum(log(alpha_1));
  // number of offspring contribution
  if (is_hetero) {
    offspring_log_intensity += sum(z .* log(nu));
  } else {
    for (i in 1:N) {
      offspring_log_intensity += z[i] * log_R[type[i]];
    }
  }
  // generation interval contribution
  for (i in 1:N_1) {
    tmp_offspring_log_intensity_vector[i] = exponential_lpdf(t_1[i] - t[bs_1[i]] | eta[type[bs_1[i]]]);
  }
  offspring_log_intensity += sum(tmp_offspring_log_intensity_vector);
  // scalar
  // generation interval contribution
  for (i in 1:N) {
    tmp_offspring_scalar_vector[i] = exponential_cdf(b - t[i], eta[type[i]]);
  }
  // time-dependant modulation
  if (K != 0) {
    tmp_offspring_scalar_vector += exp_sinusoidal_basis_integral * beta;
  }
  // reproduction number contribution
  if (is_hetero) {
    for (i in 1:N) {
      tmp_offspring_scalar_vector[i] = nu[i] * tmp_offspring_scalar_vector[i];
    }
  } else {
    for (i in 1:N) {
      tmp_offspring_scalar_vector[i] = R[type[i]] * tmp_offspring_scalar_vector[i];
    }
  }
  offspring_scalar = sum(tmp_offspring_scalar_vector);
}

model {
  // likelihood
  target += offspring_log_intensity;
  target += - offspring_scalar;
  if (is_hetero) {
    for (i in 1:N) {
      target += gamma_lpdf(nu[i] | phi[type[i]], phi[type[i]] / R[type[i]]);
    }
  }
  // priors
  target += normal_lpdf(log_R | 0, sigma_R);
  if (is_hetero) {
    target += normal_lpdf(log_phi | 0, sigma_phi);
  }
  target += normal_lpdf(log_eta | 0, sigma_eta);
  // target += uniform_lpdf(beta | -1, 1);
  target += normal_lpdf(beta | 0, 1 / (K * sqrt(2)));
}

generated quantities {
  vector[K] A; // sinusoidal amplitude
  vector[K] phase; // sinusoidal phase

  for (k in 1:K) {
    A[k] = sqrt(beta[(2*k) - 1]^2 + beta[2*k]^2);
    phase[k] = atan2(beta[2*k], beta[(2*k) - 1]);
  }
}

