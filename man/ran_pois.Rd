% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ran.R
\name{ran_pois}
\alias{ran_pois}
\title{Poisson Random Samples}
\usage{
ran_pois(n = 1, lambda = 1)
}
\arguments{
\item{n}{A non-negative whole number of the number of random samples to
generate.}

\item{lambda}{A non-negative numeric vector of means.}
}
\value{
A numeric vector of the random samples.
}
\description{
Poisson Random Samples
}
\examples{
ran_pois(10)
}
\seealso{
Other ran_dist: 
\code{\link{ran_bern}()},
\code{\link{ran_beta_binom}()},
\code{\link{ran_binom}()},
\code{\link{ran_gamma}()},
\code{\link{ran_gamma_pois}()},
\code{\link{ran_gamma_pois_zi}()},
\code{\link{ran_lnorm}()},
\code{\link{ran_neg_binom}()},
\code{\link{ran_norm}()},
\code{\link{ran_pois_zi}()},
\code{\link{ran_skewnorm}()},
\code{\link{ran_student}()}
}
\concept{ran_dist}
