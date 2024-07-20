library(tidyverse)
library(ggplot2)
library(ggh4x)

source("aux_lognormal.r")

kl_kernel <- function(x, m, s, lmu, lsd){
  dnorm(x = x, mean = m, sd = s, log = TRUE) -
    dlnorm(x = x, meanlog = lmu, sdlog = lsd, log = TRUE)
}
kl_kernel <- Vectorize(kl_kernel)

fit_both <- function(X, M = 1E6){
  ## Takes a set of draws X, fits normal (P) and lognormal (Q)
  ### computes the KL(P||Q) distance taking Q as reference
  #### returns parameter estimates and KL estimated via Monte Carlo
  #### we cheat a little bit by truncating the draws, so beware
  require(fitdistrplus)
  require(truncnorm)
  
  normal.fit <- fitdistrplus::fitdist(data = X, distr = "norm")
  lognormal.fit <- fitdistrplus::fitdist(data = X, distr = "lnorm")
  q0 <- pnorm(q = 0,
              mean = normal.fit$estimate["mean"],
              sd = normal.fit$estimate["sd"])
  us <- runif(n = M, q0, 1)
  normal.draws <- qnorm(p = us, mean = normal.fit$estimate["mean"],
                        sd = normal.fit$estimate["sd"])
  
  kls <- kl_kernel(x = normal.draws,
                   m = normal.fit$estimate["mean"],
                   s = normal.fit$estimate["sd"],
                   lmu = lognormal.fit$estimate["meanlog"],
                   lsd = lognormal.fit$estimate["sdlog"])
  
  out <- tibble::tibble(
    normal_mean = normal.fit$estimate["mean"],
    normal_sd = normal.fit$estimate["sd"],
    lognormal_mean = lognormal.fit$estimate["meanlog"],
    lognormal_sd = lognormal.fit$estimate["sdlog"],
    lognormal_skew = lognormal_skewness(lognormal.fit$estimate["sdlog"]),
    KL_MC = mean(kls, na.omit = TRUE)
  )
  return(out)
}

########################################

clade.out <- list(
  lookup = read.csv("../data/clade_lookup_Yang_golden.csv"),
  clade_log = read.csv("../data/combined_golden_Yang.log")
)

clog <- split(clade.out$clade_log, clade.out$clade_log$chain)
## removing first two obs of each chain as "burn-in"
bnin.clog <- do.call(rbind,
                      lapply(clog, function(x){
                        x[-c(1:2), ]
                      }))


only.clades <- bnin.clog[, grep("clade_", names(bnin.clog))]

Ncores <- 10
fits <- do.call(rbind,
                parallel::mclapply(1:ncol(only.clades),
                                   function(i){
                                     fit_both(X = only.clades[, i])
                                   }, mc.cores = Ncores))

fits <- data.frame(clade = names(only.clades), fits)

res <- fits[, -1]

write.csv(res,
          file = "../saved_data/distribution_fits_Yang.csv",
          row.names = FALSE)

similar <- data.frame(fits[order(fits$KL_MC), ],
                      selection = "similar")

dissimilar <- data.frame(
  fits[order(fits$KL_MC, decreasing = TRUE), ],
  selection = "dissimilar"
)

selected <- rbind(similar[1:10, ], dissimilar[1:10, ])
selected

write.csv(selected,
          file = "../saved_data/clades_of_interest_Yang.csv",
          row.names = FALSE)
