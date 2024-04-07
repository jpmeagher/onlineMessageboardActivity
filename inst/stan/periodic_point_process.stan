//
// This Stan program fits a periodic Poisson point process model to data.
//

// Functions for computing the process likelihood
functions {

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

  // Evaluate the integral of each sinusoidal basis function to provide a basis for
  // integrating over a sinusoidal additive model.
  matrix sam_integral_basis(
    int K, int N,
    real a, real b,
    vector omega
  ){
    // INPUTS:
    //   K:     The number of distinct frequencies in the sinusoidal basis
    //   N:     The number of time points at which the basis is evaluated.
    //   t:     Time points at which sinusoidal basis functions are evaluated.
    //   a:     The lower limit of integration.
    //   b:     The upper limit of integration.
    //   omega: The vector of frequencies defining the sinusoidal basis functions.
    matrix[N, 2*K] S;

    for (i in 1:N) {
      for (j in 1:K){
        S[i, (2*j) - 1] = (cos(omega[j] * a) - cos(omega[j] * b)) / omega[j];
        S[i, (2*j)] = (sin(omega[j] * b) - sin(omega[j] * a)) / omega[j];
      }
    }

    return S;
  }

  // Evaluate the integral of a sinusoidal additive model times an exponential
  // probability density function originating at time t up to time b.
  real sam_integral(
    int K,
    vector omega,
    vector alpha,
    real a,
    real b
  ){
    // INPUTS:
    //   K:     The number of distinct frequencies in the sinusoidal basis
    //   alpha: Sinusoidal basis coefficients.
    //   t:     Time points at which sinusoidal basis functions are evaluated.
    //   a:     The upper limit of integration.
    //   omega: The vector of frequencies defining the sinusoidal basis functions.
    //   xi:    The rate parameter of the exponential pdf.
    matrix[1, 2*K] S = sam_integral_basis(
        K, 1, a, b, omega
      );
    vector[1] G = S * alpha;
    real G1 = G[1];

    return(G1);
  }

  // Evaluate the likelihood of a cluster under Hawkes model given the point
  // process of event times, the underlying branching structure, the observation
  // interval, individual reproduction numbers, individual memory decay rates,
  // along with angular frequencies and coefficients allowing for exogenous
  // periodicities in the offspring intentity function.
  real periodic_poisson_process_lpdf(
    vector t,
    real lambda_0, real a, real b,
    vector omega, vector alpha
  ){
    // INPUTS:
    //   t:       Point process event times.
    //   a:       Lower limit on the observation interval
    //   b:       Upper limit on the observation interval
    //   mu:      Reproduction numbers.
    //   omega:   The vector of frequencies defining the sinusoidal basis functions.
    //   alpha:   Sinusoidal basis coefficients.

    int K = num_elements(omega);
    int N = num_elements(t);

    vector[N] periodic_intensity;
    vector[N] log_periodic_intensity;
    real compensator;
    real loglik;

    // intensity
    periodic_intensity = rep_vector(1, N);
    if (K > 0) {
    periodic_intensity +=
      sinusoidal_additive_model(K, N, alpha, t, omega);
    }
    log_periodic_intensity = rep_vector(log(lambda_0), N) + log(periodic_intensity);

    // compensator
    compensator = (b - a);
    if (K > 0) {
      compensator += sam_integral(K, omega, alpha, a,  b);
    }
    compensator = lambda_0 * compensator;

    // likelihood
    loglik = sum(log_periodic_intensity) - compensator;

    return loglik;
  }
}

// The input data.
data {
  int<lower=0> N;                        // # observations
  int<lower=0> K;                        // # sinusoidal basis frequencies
  vector[K] omega;                       // angular frequencies (2 * pi * f)
  vector[N] t;                           // event times
  real a_0;                              // immigrant arrival interval
}

// The parameters accepted by the model.
parameters {
  real<lower=0> lambda_0;    // average intensity
  vector[2*K] alpha;         // Sinusoidal basis coefficients
}

// The model to be estimated.
model {
  // Likelihood
  target += periodic_poisson_process_lpdf(
    t | lambda_0, 0, a_0, omega, alpha
  );
  // Priors
  target += normal_lpdf(alpha | 0, 1 / sqrt(2*K));
}

// Additional quantities of interest
generated quantities {
  vector[K] A;
  vector[K] phase;
  real loglik;

  for (i in 1:K) {
    A[i] = sqrt(alpha[(2*i) - 1]^2 + alpha[2*i]^2);
    phase[i] = atan2(alpha[2*i], alpha[(2*i) - 1]);
  }

  loglik = periodic_poisson_process_lpdf(
    t | lambda_0, 0, a_0, omega, alpha
  );
}
