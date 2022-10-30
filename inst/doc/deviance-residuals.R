## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(extras)
mu <- 2
sigma <- 0.5
y <- 1

(y - mu) / sigma
dev_norm(y, mu, sigma, res = TRUE)
sign(y - mu) * sqrt(dev_norm(y, mu, sigma))
sign(y - mu) * sqrt(2 * (log(dnorm(y, y, sigma)) - log(dnorm(y, mu, sigma))))

