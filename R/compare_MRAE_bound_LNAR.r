source("aux_LNAR1.r")
source("aux_fake_MCMC.r")

###############
MRAE_bound_LNAR <- function(n, m, v, rho){
  
  Mu <- exp(m + v/2)
  Sigma_sq <- (exp(v) - 1) * exp(2 * m + v)
  
  C_n <- 2 * sum(autocov_expAR1(1:(n - 1), m, v, rho))
  bound <- sqrt(Sigma_sq + C_n)/(sqrt(n) * Mu)
  
  return(bound)
}
MRAE_bound_LNAR <- Vectorize(MRAE_bound_LNAR)

###############

m <- 0 #2.58
v <-  1 # .254^2
eff <- .025
Phi <-  find_rho(eff, m, v)  ## if Phi = 0, independent samples
M <- 1e4

Mu <- exp(m + v/2)
Sigma_sq <- (exp(v) - 1) * exp(2*m + v)

NRep <- 500
simus <- parallel::mclapply(1:NRep, function(i){
  generate_fake_MCMC(
    N = M,
    mu = m,
    sigma = sqrt(v),
    phi = Phi,
    LN = TRUE
  )
}, mc.cores = 8)

Xbars <-  unlist( lapply(simus, mean) )

hist(Xbars, probability = TRUE)

RAEs <- abs(Xbars - Mu)/Mu

MRAE.brute.force <- mean(RAEs)
MRAE.brute.force
( bound <- MRAE_bound_LNAR(n = M, m = m, v = v, rho = Phi) )
sqrt(Sigma_sq/M)/Mu ## IID bound


ns <- c(100, 200, 500, 625, 1000, 2000)

plot(ns, MRAE_bound_LNAR(n = ns,
                         m = m, v = v, rho = Phi),
     ylab  = "Bound on MRAE", xlab = "Sample size (n) ")
