source("aux_prob_range.r")
source("aux_lognormal.r")
source("aux_lognormal_targets.r")
###########
sym.Mu <-  subset(target.info, target == "symmetric")$m
sym.Sigmasq <-  subset(target.info, target == "symmetric")$v
asym.Mu <-  subset(target.info, target == "moderate")$m
asym.Sigmasq <-  subset(target.info, target == "moderate")$v
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

sym.BCI <- lognormal_bci(alpha = Alpha,
                         lmean = sym.Mu,
                         lsd = sqrt(sym.Sigmasq))

asym.BCI <- lognormal_bci(alpha = Alpha,
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


pBCI <- function(k)
  prob_of_interest(logcdf = asymm_lcdf,
                   n = k,
                   a = asym.BCI[1],
                   b = asym.BCI[2])
## we could implement this directly, but this way is just convenient
pBCI <- Vectorize(pBCI)

#####################################
## The actual plot

IncProb <- .80
colours <- c("#2297E6", "#DF536B", "black")


pdf(file = "../figures/min_sample_size_cover_HPD_lognormal.pdf",   
    width = 11.7, 
    height = 8.3)

curve(
  pS,
  1,
  500,
  xlab = "Sample size (n)",
  ylab = "Probability HPD is within range",
  ylim = c(0, 1),
  main = "",
  lwd = 4,
  col = colours[1]
)
curve(pAS,
      1,
      500,
      col = colours[2],
      lwd = 4,
      add = TRUE)
curve(
  pBCI,
  1,
  500,
  col = colours[3],
  lwd = 4,
  lty = 2,
  add = TRUE
)
abline(h = IncProb, lwd = 2, lty = 3)
legend(
  x = "bottomright",
  legend = c("Symmetric", "Asymmetric"),
  col = colours,
  lwd = 3,
  bty = 'n'
)

dev.off()

optimise(function(n)
  (pS(n) - IncProb) ^ 2,
  c(1, 1500))

optimise(function(n)
  (pAS(n) - IncProb) ^ 2,
  c(1, 2000))
