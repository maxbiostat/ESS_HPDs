## First calculation: everything is known
### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
sym <- FALSE
logN <- FALSE
if(sym){
  Mu <- 4.01 
  Sigma <- 0.0352 
}else{
  Mu <- 2.58
  Sigma <- 0.254 
}
if(logN){
  true_F <- function(x) plnorm(q = x, Mu, Sigma)
  true_Finv <- function(p) qlnorm(p, Mu, Sigma)
  true_f <- function(x) dlnorm(x = x, Mu, Sigma)
}else{
  true_F <- function(x) pnorm(q = x, Mu, Sigma)
  true_Finv <- function(p) qnorm(p, Mu, Sigma)
  true_f <- function(x) dnorm(x = x, Mu, Sigma)
}

q <- .95
true.csi <- true_Finv(q) ## true quantile 
ci.level <- .95
delta <- 0.005
Z <- qnorm(p = (1-ci.level)/2)
R <- (Z^2 * q *(1-q))/(delta * true.csi * true_f(true.csi))^2
S <- floor(R) + 1
S 