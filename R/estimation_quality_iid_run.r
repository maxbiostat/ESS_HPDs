library(tidyverse)
source("aux.r")
##############################
############# Functions
generate_samples <- function(N, mu, sigma, LN = TRUE) {
  if (LN) {
    samples <- rlnorm(n = N,
                      meanlog = mu,
                      sdlog = sigma)
  } else{
    samples <- rnorm(n = N, mean = mu, sd = sigma)
  }
  return(samples)
}
run_once <- function(i, Mu, Sigma, Alpha, M, logn) {
  Samples <- generate_samples(
    N = M,
    mu = Mu,
    sigma = Sigma,
    LN = logn
  )
  #
  HPD <- HDInterval::hdi(Samples, credMass = Alpha)
  BCI <- quantile(Samples, probs = c(1 - Alpha, 1 + Alpha) / 2)
  xbar <- mean(Samples)
  md <- median(Samples)
  #
  out <- tibble::tibble(
    min = min(Samples),
    max = max(Samples),
    point = c(HPD, BCI, xbar, md),
    quantity = c(
      "HPD_L",
      "HPD_U",
      "BCI_L",
      "BCI_U",
      "sample_mean",
      "sample_median"
    ),
    type = "IID",
    M = M,
    replicate = i
  )
  return(out)
}

do_for_M <- function(M) {
  simus <- do.call(rbind,
                   parallel::mclapply(seq_len(Nrep),
                                      function(i) {
                                        run_once(
                                          i,
                                          M = M,
                                          Mu = Mu,
                                          Sigma = Sigma,
                                          Alpha = Alpha,
                                          logn = logNormal
                                        )
                                      }, mc.cores = 6))
  return(simus)
}

### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
## should the target be log-normal (TRUE) or normal (FALSE)
logNormal <- TRUE

Mu <-  0 # 2.58
Sigma <- 1 # 0.254

Alpha <- 0.95
M <- 4E3
Nrep <- 500

if (logNormal) {
  true.Mean <- exp(Mu + Sigma ^ 2 / 2)
  true.Median <- exp(Mu)
  true.BCI <- qlnorm(p = c(1 - Alpha, 1 + Alpha) / 2,
                     meanlog = Mu,
                     sdlog = Sigma)
  true.HPD <- lognormal_hpd(alpha = Alpha,
                            lmean = Mu,
                            lsd = Sigma)
} else{
  true.Mean <- Mu
  true.Median <- Mu
  true.HPD <- qnorm(p = c(1 - Alpha, 1 + Alpha) / 2,
                    mean = Mu,
                    sd = Sigma)
  true.BCI <- true.HPD
}

true.quants <- tibble::tibble(
  true = c(true.HPD, true.BCI, true.Mean, true.Median),
  quantity = c(
    "HPD_L",
    "HPD_U",
    "BCI_L",
    "BCI_U",
    "sample_mean",
    "sample_median"
  ),
)

Ms <- c(100, 200, 625, 1000)

results <- do.call(rbind, lapply(Ms, do_for_M))
results.df <- merge(results, true.quants, by = "quantity")

write.csv(results.df, 
          file = paste0("../saved_data/iid_lognormal_m=", Mu,
          "_s=", Sigma, "_results.csv"), 
          row.names = FALSE)