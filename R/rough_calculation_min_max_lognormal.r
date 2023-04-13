source("aux_prob_range.r")
source("aux.r")
###########
### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
sym.Mu <- 4.01
sym.Sigmasq <- 0.0352 ^ 2
asym.Mu <- 2.58
asym.Sigmasq <- 0.254 ^ 2
Alpha <- .95
symm_lcdf <- function(x)
  plnorm(
    q = x,
    meanlog =  sym.Mu,
    sdlog =  sqrt(sym.Sigmasq),
    log.p = TRUE
  )
asymm_lcdf <- function(x)
  plnorm(
    q = x,
    meanlog = asym.Mu,
    sdlog =  sqrt(asym.Sigmasq),
    log.p = TRUE
  )

sym.HPD <- lognormal_hpd(alpha = Alpha,
                         lmean = sym.Mu,
                         lsd = sqrt(sym.Sigmasq))

asym.HPD <- lognormal_hpd(alpha = Alpha,
                          lmean = asym.Mu,
                          lsd = sqrt(asym.Sigmasq))

pS <- function(k)
  prob_of_interest(logcdf = symm_lcdf,
                   n = k,
                   a = sym.HPD[1],
                   b = sym.HPD[2])
pS <- Vectorize(pS)

pAS <- function(k)
  prob_of_interest(logcdf = asymm_lcdf,
                   n = k,
                   a = asym.HPD[1],
                   b = asym.HPD[2])
pAS <- Vectorize(pAS)

IncProb <- .80
curve(
  pS,
  1,
  500,
  xlab = "Sample size (n)",
  ylab = "Probability HPD is within range",
  main = "",
  lwd = 4,
  col = "blue"
)
curve(pAS,
      1,
      500,
      col = "gold",
      lwd = 4,
      add = TRUE)
abline(h = IncProb, lwd = 2, lty = 2)
legend(
  x = "bottomright",
  legend = c("Symmetric", "Asymmetric"),
  col = c("blue", "gold"),
  lwd = 3,
  bty = 'n'
)

optimise(function(n)
  (pS(n) - IncProb) ^ 2,
  c(1, 500))

optimise(function(n)
  (pAS(n) - IncProb) ^ 2,
  c(1, 5E4))
