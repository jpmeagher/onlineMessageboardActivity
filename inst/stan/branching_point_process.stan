//
// This Stan program fits a branching point process model to data.
//

// Functions for computing the cluster likelihood
functions {
  // count the number of offspring from each point
  int[] offspring_count(
    int N, int[] beta
  ){
    int z[N];

    for (i in 1:N) {
      z[i] = 0;
      for (j in 1:N) {
        z[i] += beta[j] == i;
      }
    }
    return z;
  }

  // Evaluates the sinusoidal basis functions specified by
  // angular frequencies omega at time points t.
  matrix sinusoidal_basis(
    int K, int N, vector t, vector omega
  ){
    // INPUTS:
    //   K:     The number of distinct frequencies in the sinusoidal basis
    //   N:     The number of time points at which the basis is evaluated.
    //   t:     Time points at which sinusoidal basis functions are evaluated.
    //   omega: The vector of frequencies defining the sinusoidal basis functions.
    matrix[N, 2*K] X;

    for (j in 1:K) {
      X[, (2*j) - 1] = sin(omega[j] * t);
      X[, 2*j] = cos(omega[j] * t);
    }

    return X;
  }

  // Evaluates the sinusoidal additive model specified by
  // coefficients beta and angular frequencies omega at time points t.
  vector sinusoidal_additive_model(
    int K, int N, vector alpha, vector t, vector omega
  ){
    // INPUTS:
    //   K:     The number of distinct frequencies in the sinusoidal basis
    //   N:     The number of time points at which the basis is evaluated.
    //   alpha: Sinusoidal basis coefficients.
    //   t:     Time points at which sinusoidal basis functions are evaluated.
    //   omega: The vector of frequencies defining the sinusoidal basis functions.
    matrix[N, 2*K] X = sinusoidal_basis(K, N, t, omega);
    vector[N] f = X * alpha;

    return(f);
  }

  // Evaluate the integral of each sinusoidal basis function times an
  // exponential probability density function to provide a basis for
  // integrating over the sinusoidal additive model.
  matrix sam_exp_prod_integral_basis(
    int K, int N,
    vector t, real a,
    vector omega, vector xi
  ){
    // INPUTS:
    //   K:     The number of distinct frequencies in the sinusoidal basis
    //   N:     The number of time points at which the basis is evaluated.
    //   t:     Time points at which sinusoidal basis functions are evaluated.
    //   a:     The upper limit of integration.
    //   omega: The vector of frequencies defining the sinusoidal basis functions.
    //   xi:    The rate parameter of the exponential pdf.
    matrix[N, 2*K] W;

    for (i in 1:N) {
      for (j in 1:K){
        W[i, (2*j) - 1] = xi[i] / (xi[i]^2 + omega[j]^2) * (
          xi[i] * (sin(omega[j] * t[i]) - (sin(omega[j] * a) * exp( - xi[i] * (a - t[i])))) +
          omega[j] * (cos(omega[j] * t[i]) - (cos(omega[j] * a) * exp( - xi[i] * (a - t[i]))))
        );
        W[i, (2*j)] = xi[i] / (xi[i]^2 + omega[j]^2) * (
          omega[j] * (sin(omega[j] * a) * exp( - xi[i] * (a - t[i])) - (sin(omega[j] * t[i]))) -
          xi[i] * (cos(omega[j] * a) * exp( - xi[i] * (a - t[i])) - (cos(omega[j] * t[i])))
        );
      }
    }

    return W;
  }

  // Evaluate the integral of a sinusoidal additive model times an exponential
  // probability density function originating at time t up to time b.
  vector sam_exp_prod_integral(
    int K, int N,
    vector alpha,
    vector t, real a,
    vector omega, vector xi
  ){
    // INPUTS:
    //   K:     The number of distinct frequencies in the sinusoidal basis
    //   N:     The number of time points at which the basis is evaluated.
    //   alpha: Sinusoidal basis coefficients.
    //   t:     Time points at which sinusoidal basis functions are evaluated.
    //   a:     The upper limit of integration.
    //   omega: The vector of frequencies defining the sinusoidal basis functions.
    //   xi:    The rate parameter of the exponential pdf.
    matrix[N, 2*K] W = sam_exp_prod_integral_basis(
        K, N, t, a, omega, xi
      );
    vector[N] F = W * alpha;

    return(F);
  }

  // Evaluate the likelihood of a cluster under Hawkes model given the point
  // process of event times, the underlying branching structure, the observation
  // interval, individual reproduction numbers, individual memory decay rates,
  // along with angular frequencies and coefficients allowing for exogenous
  // periodicities in the offspring intentity function.
  real cluster_lpdf(
    vector t, int[] beta, int[] z,
    real a,
    vector mu,
    vector psi,
    vector eta,
    int[] inf_idx,
    vector omega, vector alpha
  ){
    // INPUTS:
    //   t:       Point process event times.
    //   beta:    Underlying branching structure.
    //   z:       Number of offspring from each point.
    //   a:       Upper limit on the observation interval
    //   mu:      Reproduction numbers.
    //   psi:     Dispersion parameters.
    //   eta:     The rate parameters of the exponential pdf.
    //   inf_idx: Index of infinite dispersion parameters.
    //   omega:   The vector of frequencies defining the sinusoidal basis functions.
    //   alpha:   Sinusoidal basis coefficients.

    int K = num_elements(omega);
    int N = num_elements(t);

    vector[N-1] alpha_fun;
    vector[N-1] l_alpha_fun;
    vector[N-1] l_rho;

    vector[N] C_rate;
    vector[N] l_norm;
    vector[N] l_off;
    vector[N] l_disp;
    real ll;
    // Point intensity
    alpha_fun = rep_vector(1, N-1);
    if ((K > 0) && (N > 1)) {
    alpha_fun += sinusoidal_additive_model(
    K, N-1, alpha, t[2:N], omega
    );
    }
    l_alpha_fun = log(alpha_fun);
    for (i in 2:N) l_rho[i-1] = exponential_lpdf(t[i] - t[beta[i]] | eta[beta[i]]);
    // Rate scalar
    for (i in 1:N) C_rate[i] = exponential_cdf(a - t[i], eta[i]);
    if (K > 0) {
      C_rate += sam_exp_prod_integral(
        K, N, alpha, t, a, omega, eta
      );
    }
    // neg-bin likelihood
    for (i in 1:N){

      l_off[i] = z[i] * log(mu[i]);

      if (inf_idx[i] == 0) {
        l_norm[i] = lgamma(z[i] + psi[i]) - lgamma(psi[i]) - z[i] * log(mu[i] * C_rate[i] + psi[i]);
        l_disp[i] = - psi[i] * log1p(mu[i] * C_rate[i] / psi[i]);
      } else {
        l_norm[i] = 0;
        l_disp[i] = - mu[i] * C_rate[i];
      }
    }

    ll = sum(l_alpha_fun + l_rho) + sum(l_norm + l_off + l_disp);

    return ll;
  }
}

// The input data.
data {
  int<lower=0> N;                        // # observations
  int<lower=0> S;                        // # clusters
  int<lower=0> K;                        // # sinusoidal basis frequencies
  int<lower=1> M;                        // # distinct types of offspring process
  vector<lower=0>[S] a;                  // cluster-wise observation horizon
  vector[K] omega;                       // angular frequencies (2 * pi * f)
  vector[N] t;                           // event times
  int beta[N];                           // within-cluster branching structure
  int n_i[S];                            // cluster sizes
  int<lower=1, upper=M> type[N];         // type of event
  int<lower=0, upper=1> is_hetero[M];    // is there heterogeneous reproduction
  real<lower=0> alpha_eta;               // gamma prior shape for eta
  real<lower=0> beta_eta;                // gamma prior rate for eta
  real<lower=0> alpha_mu;                // gamma prior shape for mu
  real<lower=0> beta_mu;                 // gamma prior rate for mu
}

// Data transformations for easy computation
transformed data {
  int M_hetero = sum(is_hetero); // Number of heterogeneous types
  int j = 1;

  int pos[S];
  int<lower=0> z[N];
  int hetero_idx[M_hetero];
  int<lower=0, upper=1> inf_idx[N];
  // Starting index for each cluster
  pos[1] = 1;
  for (i in 2:S) {
    pos[i] = pos[i-1] + n_i[i-1];
  }
  // how many offspring
  for (i in 1:S) {
    z[pos[i]:(pos[i] + n_i[i] - 1)] = offspring_count(n_i[i], segment(beta, pos[i], n_i[i]));
  }
  // Index the heterogeneous type
  for (i in 1:M) {
    if (is_hetero[i] == 1) {
      hetero_idx[j] = i;
      j += 1;
    }
  }
  // index points belonging to homogeneous types
  for (i in 1:N) {
    inf_idx[i] = is_hetero[type[i]] == 0;
  }
}

// The parameters accepted by the model.
parameters {
  vector<lower=0>[M] mu;     // Reproduction numbers
  vector[M_hetero] log_psi;  // Log dispersion parameters
  vector<lower=0>[M] eta;    // Memory decay rates
  vector[2*K] alpha;         // Sinusoidal basis coefficients
}

// Transform parameters for easy computation
transformed parameters {
  vector<lower=0>[M] psi;

  for (i in 1:M) {
    psi[i] = positive_infinity();
  }
  for (i in 1:M_hetero) {
    psi[hetero_idx[i]] =  exp(log_psi[i]);
  }
}

// The model to be estimated.
model {
  // Likelihood
  for (i in 1:S) {
    target += cluster_lpdf(
      segment(t, pos[i], n_i[i]) |
      segment(beta, pos[i], n_i[i]),
      segment(z, pos[i], n_i[i]),
      a[i],
      mu[segment(type, pos[i], n_i[i])],
      psi[segment(type, pos[i], n_i[i])],
      eta[segment(type, pos[i], n_i[i])],
      segment(inf_idx, pos[i], n_i[i]),
      omega, alpha
    );
  }
  // Priors
  target += gamma_lpdf(mu | alpha_mu, beta_mu);
  target += normal_lpdf(log_psi | 0, 1);
  target += gamma_lpdf(eta | alpha_eta, beta_eta);
  target += normal_lpdf(alpha | 0, 1 / sqrt(2*K));
}

// Additional quantities of interest
generated quantities {
  vector[K] A;
  vector[K] phase;

  for (i in 1:K) {
    A[i] = sqrt(alpha[(2*i) - 1]^2 + alpha[2*i]^2);
    phase[i] = atan2(alpha[2*i], alpha[(2*i) - 1]);
  }
}
