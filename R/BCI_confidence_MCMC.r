library(tidyverse)
source("aux.r")
source("aux_fake_MCMC.r")
source("aux_estimation_quality.r")
source("aux_linear_mixture.r")
source("aux_quantile_CI.r")

##############################
############# Functions
source("aux_BCI_estimators.r")
##
run_once <- function(i, Mu, Sigma, Phi, Alpha, Omega, M, logn) {
  Samples <- generate_fake_MCMC(
    N = M,
    mu = Mu,
    sigma = Sigma,
    phi = Phi,
    LN = logn
  )
  
  Meeker.ests <- get_bci_intervals_Meeker(X = Samples,
                                   b_level = Alpha,
                                   ci_level = Omega)
  
  Stan.ests <- get_bci_intervals_Stan(X = Samples,
                                          b_level = Alpha,
                                          ci_level = Omega)
  
  Doss.ests <- get_bci_intervals_Doss(X = Samples,
                                      b_level = Alpha,
                                      ci_level = Omega)
  
  Batches.ests <- get_bci_intervals_batches(X = Samples,
                                            b_level = Alpha,
                                            ci_level = Omega)
  
  concat <- do.call(
    rbind,
    list(
      Meeker.ests,
      Doss.ests,
      Stan.ests,
      Batches.ests
    )
  )
  
  out <- tibble::tibble(
    min = min(Samples),
    max = max(Samples),
    concat,
    replicate = i
  )
  return(out)
}
##########################################
### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
## should the target be log-normal (TRUE) or normal (FALSE)?
logNormal <- TRUE

Mu <- 0
Sigma <- 1
Type <- "high"
eff <- .05
Phi <- (1 - eff) / (eff + 1) ## if Phi = 0, independent samples
exper <- 9

Alpha <- 0.95
TargetCoverage <- .95
M <- 4E3
Nrep <- 1e4

if (logNormal) {
  true.BCI <- qlnorm(
    p = c(1 - Alpha, 1 + Alpha) / 2,
    meanlog = Mu,
    sdlog =  Sigma
  )
} else{
  true.BCI <- qnorm(
    p = c(1 - Alpha, 1 + Alpha) / 2,
    mean = Mu,
    sd =  Sigma
  )
}

simus <- do.call(rbind,
                 parallel::mclapply(seq_len(Nrep),
                                    function(i) {
                                      run_once(
                                        i,
                                        M = M,
                                        Mu = Mu,
                                        Sigma = Sigma,
                                        Phi = Phi,
                                        Alpha = Alpha,
                                        Omega = TargetCoverage,
                                        logn = logNormal
                                      )
                                    }, mc.cores = 6))

true.vals <- tibble::tibble(true_value = true.BCI,
                            quantity = c("BCI_L", "BCI_U"))
simus.df <- merge(simus, true.vals, by = "quantity")
simus.df <-
  simus.df %>% mutate(covers = is_in(true_value, lwr, upr))
simus.df <-
  simus.df %>% mutate(lwr = ifelse(is.infinite(lwr), 0, lwr))

coverage <-
  aggregate(covers ~ quantity + method, simus.df, mean) ## coverage
ci.width <-
  aggregate((upr - lwr) ~ quantity + method, simus.df, mean) ## average CI width
bias <-
  aggregate((point - true_value) ~ quantity + method, simus.df, mean)
variance <- aggregate(point ~ quantity + method, simus.df, var)
mse <-
  aggregate((point - true_value) ^ 2 ~ quantity + method, simus.df, mean)
ans <-
  Reduce(function(...)
    merge(..., by = c("quantity", "method")),
    list(coverage, ci.width,
         bias,
         variance, mse))
names(ans) <- c("quantity",
                "method",
                "est_coverage",
                "ci_width",
                "bias",
                "variance",
                "MSE")

ans$mu <- Mu
ans$sigma <- Sigma
ans$efficiency <- eff
ans$alpha <- Alpha
ans$ci_target_coverage <- TargetCoverage
ans$n_draws <- M
ans$asymmetry <- Type

ans

write.csv(ans,
          file = paste0("BCI_experiment_", exper, ".csv"),
          row.names = FALSE)

pp <- ggplot(simus.df,
             aes(x = point, fill = method)) +
  geom_density(alpha = 0.4) +
  facet_wrap(. ~ quantity, scales = "free") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(data = true.vals,
             aes(xintercept = true_value),
             linetype = "longdash") +
  theme_bw()

pp

ggsave(
  plot = pp,
  filename =
    paste0("../figures/BCI_experiment_",
           exper,
           ".pdf"),
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)

# ggplot(simus.df,
#        aes(x = point, fill = method)) +
#   geom_boxplot() +
#   facet_wrap(. ~ quantity, scales = "free") +
#   geom_vline(data = true.vals,
#              aes(xintercept = true_value),
#              linetype = "longdash") +
#   theme_bw()
