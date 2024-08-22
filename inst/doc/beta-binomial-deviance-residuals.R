## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(extras)
library(tidyr)
library(ggplot2)
library(viridis)
library(scales)

## ----echo = FALSE, fig.width = 7, fig.height = 4, fig.cap = "Fig. 1: Beta-binomial likelihood profile for $x = 1$ and $n = 5$, for different values of $\\theta$. $\\theta = 0$ corresponds to the binomial case."----
n_samp <- 1000
x <- 1
size <- 5
prob <- 0.3
theta <- c(0, 0.1, 0.5, 1)

prob_seq <- seq(0, 1, length.out = n_samp)

lik <- data.frame(
  prob_seq = prob_seq,
  lik_0.0 = exp(log_lik_beta_binom(x, size, prob_seq, theta[1])),
  lik_0.1 = exp(log_lik_beta_binom(x, size, prob_seq, theta[2])),
  lik_0.5 = exp(log_lik_beta_binom(x, size, prob_seq, theta[3])),
  lik_1.0 = exp(log_lik_beta_binom(x, size, prob_seq, theta[4]))
)

tb <- as_tibble(lik)
tb <- pivot_longer(tb,
  cols = c(lik_0.0, lik_0.1, lik_0.5, lik_1.0),
  names_to = "theta", names_prefix = "lik_",
  values_to = "lik"
)

ggplot(tb, aes(x = prob_seq, y = lik, colour = theta)) +
  geom_line() +
  scale_colour_viridis(discrete = TRUE) +
  xlab("p") +
  ylab("likelihood")

## ----echo = FALSE, fig.width = 7, fig.height = 4, fig.cap = "Fig. 2: Beta-binomial likelihood profile for $x = 1$, $n = 5$, and $\\theta = 0.5$. The dashed vertical line shows the likelihood at the expected value of the beta-binomial distribution ($n \\cdot p$) where $p = \\frac{x}{n} = \\frac{1}{5} = 0.2$. As you can see, this is not the $p$ for which the likelihood is maximized. The solid vertical line shows the likelihood at its maximum point. In this case, $p^* = 0.26$."----
n_samp <- 1000
x <- 1
size <- 5
prob <- 0.3
theta <- 0.5

opt_beta_binom <- function(prob, x, size, theta) {
  -log_lik_beta_binom(x = x, size = size, prob = prob, theta = theta)
}

opt_p <- stats::optimize(opt_beta_binom,
  interval = c(0, 1), x = x,
  size = size, theta = theta
)$minimum

prob_seq <- seq(0, 1, length.out = n_samp)

lik <- data.frame(
  prob_seq = prob_seq,
  lik = exp(log_lik_beta_binom(x, size, prob_seq, theta))
)
p_df <- data.frame(
  val = c(x / size, opt_p),
  probability = c("Expected value (x / n) = 0.20", "Optimized value (p*) = 0.26")
)

ggplot(lik, aes(x = prob_seq, y = lik)) +
  geom_line(colour = viridis(1, begin = 0.66)) +
  geom_vline(data = p_df, mapping = aes(linetype = probability, xintercept = val)) +
  xlab("p") +
  ylab("likelihood")

## ----echo = FALSE, fig.height = 8, fig.width = 7, fig.cap = "Fig. 4: Histograms of deviance residuals for 10,000 beta-binomial data points simulated with $n = 50$, $p = \\{0.3, 0.5, 0.9 \\}$, and $\\theta = \\{0.0, 0.1, 0.5, 1.0 \\}$. The percentages within each panel (i.e., each combination of $p$ and $\\theta$) sum to 100%."----
n_samp <- 10000
size <- 50
prob <- c(0.3, 0.5, 0.9)
theta <- c(0, 0.1, 0.5, 1)

res_mat <- matrix(NA, nrow = n_samp * length(theta), ncol = length(prob))

set.seed(101)
for (i in seq_along(prob)) {
  res <- NULL
  for (j in seq_along(theta)) {
    x <- ran_beta_binom(n_samp, size, prob[i], theta[j])
    res <- c(res, res_beta_binom(x, size, prob[i], theta[j]))
  }
  res_mat[, i] <- res
}

df <- as.data.frame(res_mat)
colnames(df) <- c("prob_0.3", "prob_0.5", "prob_0.9")
df$theta <- rep(theta, each = n_samp)

df <- pivot_longer(df,
  cols = c(prob_0.3, prob_0.5, prob_0.9), names_to = "prob",
  names_prefix = "prob_", values_to = "res"
)

ggplot(df, aes(x = res)) +
  geom_histogram(aes(y = after_stat(density) * 0.1), binwidth = 0.1) +
  scale_y_continuous(labels = percent) +
  facet_grid(~prob ~ theta, labeller = label_bquote(cols = theta:.(theta), rows = p:.(prob))) +
  xlab("Deviance residuals") +
  ylab("Percent (%)")

