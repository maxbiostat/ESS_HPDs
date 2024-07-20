source("aux_LNAR1.r")
source("aux_fake_MCMC.r")

library(cmdstanr)

m <- -2
v <- 1^2 #0.254
eff <- .02
Phi <- find_rho(eff)
M <- 1e4

Samples <- generate_fake_MCMC(
  N = M,
  m = m,
  sigma = sqrt(v),
  phi = Phi,
  LN = TRUE
)

hist(Samples, probability = TRUE)
curve(dlnorm(x, m, sqrt(v)), min(Samples), max(Samples),
      lwd = 3, add = TRUE)

LNAR_mod <- cmdstan_model("../stan/LNAR_exp.stan")

LNAR.MAP <- LNAR_mod$optimize(data = list(N = M, Y = Samples))

# LNAR.MCMC <- LNAR_mod$sample(data = list(N = M, Y = Samples),
#                            chains = 4,
#                            parallel_chains = 4)

source("aux_LNAR1.r")

m.hat <- LNAR.MAP$mle("m")
s.hat <- LNAR.MAP$mle("s")
Rho.hat <- LNAR.MAP$mle("rho")

LNAR1_eff(m, v, rho = Phi)
LNAR1_eff(m.hat, s.hat^2, rho = Rho.hat)

mean(Samples)
exp(m + v/2)
exp(m.hat + s.hat^2/2)
## Conclusion: model is not identifiable.
## so we can just do this

eff.hat <- posterior::ess_basic(Samples)/M
Rho.hat.ESS <- (1 - eff.hat) / (eff.hat + 1)


Phi
Rho.hat
Rho.hat.ESS

ks <- 0:50

rs.true <- autocorr_expAR1(ks, m, v, Phi)
rs.hat <- autocorr_expAR1(ks, m.hat, s.hat^2, Rho.hat)

plot(rs.hat ~ rs.true)

mean(log(Samples))
m
sd(log(Samples))
sqrt(v)

