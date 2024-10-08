% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-lik.R
\name{log_lik_beta_binom}
\alias{log_lik_beta_binom}
\title{Beta-Binomial Log-Likelihood}
\usage{
log_lik_beta_binom(x, size = 1, prob = 0.5, theta = 0)
}
\arguments{
\item{x}{A non-negative whole numeric vector of values.}

\item{size}{A non-negative whole numeric vector of the number of trials.}

\item{prob}{A numeric vector of values between 0 and 1 of the probability of
success.}

\item{theta}{A non-negative numeric vector of the dispersion for the mixture
models (student, gamma-Poisson and beta-binomial).}
}
\value{
An numeric vector of the corresponding log-likelihoods.
}
\description{
This parameterization of the beta-binomial distribution uses an expected
probability parameter, \code{prob}, and a dispersion parameter, \code{theta}. The
parameters of the underlying beta mixture are \code{alpha = (2 * prob) / theta}
and \code{beta = (2 * (1 - prob)) / theta}. This parameterization of \code{theta} is
unconventional, but has useful properties when modelling. When \code{theta = 0},
the beta-binomial reverts to the binomial distribution. When \code{theta = 1} and
\code{prob = 0.5}, the parameters of the beta distribution become \code{alpha = 1} and
\code{beta = 1}, which correspond to a uniform distribution for the beta-binomial
probability parameter.
}
\examples{
log_lik_beta_binom(c(0, 1, 2), 3, 0.5, 0)
}
\seealso{
Other log_lik_dist: 
\code{\link{log_lik_bern}()},
\code{\link{log_lik_binom}()},
\code{\link{log_lik_gamma}()},
\code{\link{log_lik_gamma_pois}()},
\code{\link{log_lik_gamma_pois_zi}()},
\code{\link{log_lik_lnorm}()},
\code{\link{log_lik_neg_binom}()},
\code{\link{log_lik_norm}()},
\code{\link{log_lik_pois}()},
\code{\link{log_lik_pois_zi}()},
\code{\link{log_lik_skewnorm}()},
\code{\link{log_lik_student}()}
}
\concept{log_lik_dist}
