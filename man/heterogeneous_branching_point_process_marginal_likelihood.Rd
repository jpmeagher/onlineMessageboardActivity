% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal_likelihood.R
\name{heterogeneous_branching_point_process_marginal_likelihood}
\alias{heterogeneous_branching_point_process_marginal_likelihood}
\title{Heterogeneous Branching Point Process Marginal Likelihood}
\usage{
heterogeneous_branching_point_process_marginal_likelihood(
  t,
  beta,
  a,
  mu_vec,
  psi_vec,
  eta_vec,
  omega = NULL,
  alpha = NULL,
  z = NULL,
  log = TRUE,
  perform_checks = TRUE
)
}
\arguments{
\item{t}{A real-valued vector of length N. The times of the point process.}

\item{beta}{An integer vector of length N. The branching structure of the
point process.}

\item{a}{A real valued scalar. The upper limit of the observation interval.}

\item{mu_vec}{A positive real-valued vector of length N. The reproduction
number for each point in the cluster.}

\item{psi_vec}{A positive real-valued vector of length N. The dispersion
parameter for each point in the cluster.}

\item{eta_vec}{A positive real-valued vector of length N. The memory decay
rate for each point in the cluster.}

\item{omega}{A positive, increasing real-valued vector of length K.
Frequencies defining the sinusoidal basis functions for exogeneous
excitation.}

\item{alpha}{A real-valued vector of length 2K. Coefficients of the
sinusoidal basis function for exogenous excitation.}

\item{z}{A non-negative integer-valued vector of length N. The number of
offspring from each point.}

\item{log}{A logical scalar. If `TRUE` then the log-likelihood is returned.}

\item{perform_checks}{A logical scalar. Should inputs be checked before
executing code.}
}
\value{
The (log)-likelihood for the branching point process.
}
\description{
Computes the heterogeneous branching point process marginal likelihood for
the class of branching processes described in Meagher & Friel.
}
