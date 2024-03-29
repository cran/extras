% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pzeros.R
\name{pzeros}
\alias{pzeros}
\title{Proportion of Zeros}
\usage{
pzeros(x, na_rm = FALSE)
}
\arguments{
\item{x}{A numeric object of MCMC values.}

\item{na_rm}{A flag specifying whether to remove missing values.}
}
\value{
A number between 0 and 1.
}
\description{
The proportion of zeros in an numeric object.
}
\examples{
pzeros(c(0:2))
}
\seealso{
Other summary: 
\code{\link{kurtosis}()},
\code{\link{lower}()},
\code{\link{pvalue}()},
\code{\link{skewness}()},
\code{\link{svalue}()},
\code{\link{upper}()},
\code{\link{variance}()},
\code{\link{xtr_mean}()},
\code{\link{xtr_median}()},
\code{\link{xtr_sd}()},
\code{\link{zeros}()},
\code{\link{zscore}()}
}
\concept{summary}
