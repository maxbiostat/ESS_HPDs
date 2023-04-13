## First calculation: everything is known
### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
sym <- FALSE
if (sym) {
  Mu <- 4.01
  Sigma <- 0.0352
} else{
  Mu <- 0 # 2.58
  Sigma <- 1# 0.254
}

true_F <- function(x)
  plnorm(q = x, Mu, Sigma)
true_Finv <- function(p)
  qlnorm(p, Mu, Sigma)
true_f <- function(x)
  dlnorm(x = x, Mu, Sigma)

source("aux.r")
(HPD <- lognormal_hpd(alpha = 0.95, lmean = Mu, lsd = Sigma))
(qs <- true_F(HPD))

curve(true_f, HPD[1], HPD[2], lwd = 3)

q <- qs[1]
true.csi <- true_Finv(q) ## true quantile
ci.level <- .95
delta <- .05 ## relative precision
Z <- qnorm(p = (1 - ci.level) / 2)
R <- (Z ^ 2 * q * (1 - q)) / (delta * true.csi * true_f(true.csi)) ^ 2
S <- floor(R) + 1
S
true.csi
delta * true.csi