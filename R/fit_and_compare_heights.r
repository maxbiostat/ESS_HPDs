library(tidyverse)
library(ggplot2)
library(ggh4x)
source("aux.r")

kl_kernel <- function(x, m, s, lmu, lsd){
  dnorm(x = x, mean = m, sd = s, log = TRUE) -
    dlnorm(x = x, meanlog = lmu, sdlog = lsd, log = TRUE)
}
kl_kernel <- Vectorize(kl_kernel)

fit_both <- function(X, M = 1E6){
  ## Takes a data set X, fits normal (P) and lognormal (Q)
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
    KL_MC = mean(kls, na.omit = TRUE)
  )
  return(out)
}
make_all_fits <- function(fname, Ncores = 6){
  
  raw.log <- read.table(fname, header = TRUE)
  only.clades <- raw.log[, -1]
  cladenames <- colnames(only.clades)
  cladenames <- gsub("yang.", "", cladenames)
  newNames <- paste0("clade_", match(cladenames, sort(cladenames)))
  lookup <- tibble::tibble(
    original = cladenames,
    new = newNames
  )
  
  raw.fits <- parallel::mclapply(1:ncol(only.clades),
                                 function(i){
                                   fit_both(X = only.clades[, i])
                                 }, mc.cores = Ncores)
  
  ans <- tibble::tibble(
    lookup,
    do.call(rbind, raw.fits)
  )
  return(
    ans
  )
}
#########

the.file <- "../data/yang-mcc-clade-heights.log"
fits <- make_all_fits(fname = the.file, Ncores = 6)
similar <- fits[order(fits$KL_MC), -1]
dissimilar <- fits[order(fits$KL_MC, decreasing = TRUE), -1]

selected <- rbind(similar[1:10, ], dissimilar[1:10, ])
selected

write.csv(selected,
          file = "../saved_data/clades_of_interest_Yang.csv",
          row.names = FALSE)