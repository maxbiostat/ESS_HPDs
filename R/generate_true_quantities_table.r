library(tibble)
source("aux_true_quantities.r")
source("aux_lognormal_targets.r")

Alpha <- 0.95

K <- nrow(target.info)

BCIs <- vector(length = K, mode = "list")
HPDs <- vector(length = K, mode = "list")
BCI.ps <- vector(length = K, mode = "list")
HPD.ps <- vector(length = K, mode = "list")

means <-  vector(length = K, mode = "numeric")
medians <- vector(length = K, mode = "numeric")
vars <- vector(length = K, mode = "numeric")

for(k in 1:K){
  
  BCIs[[k]] <- lognormal_bci(alpha = Alpha,
                             lmean = target.info$m[k],
                             lsd = target.info$sigma[k])
  
  HPDs[[k]] <- lognormal_hpd(alpha = Alpha,
                             lmean = target.info$m[k],
                             lsd = target.info$sigma[k])
  
  BCI.ps[[k]] <- plnorm(q = BCIs[[k]],
                        meanlog = target.info$m[k],
                        sdlog = target.info$sigma[k]) ## could be skipped since it's just (1-alpha, 1+alpha)/2
  
  HPD.ps[[k]] <- plnorm(q = HPDs[[k]],
                        meanlog = target.info$m[k],
                        sdlog = target.info$sigma[k])
  
  means[k] <- lognormal_mean(lmean = target.info$m[k],
                             lsd = target.info$sigma[k])
  
  medians[k] <- lognormal_median(lmean = target.info$m[k],
                                 lsd = target.info$sigma[k])
  
  vars[k] <- lognormal_variance(lmean = target.info$m[k],
                                lsd = target.info$sigma[k])
}

bci.df <- data.frame(do.call(rbind, BCIs))
hpd.df <- data.frame(do.call(rbind, HPDs))

colnames(bci.df) <- c("bci_lwr", "bci_upr")
colnames(hpd.df) <- c("hpd_lwr", "hpd_upr")

bci.p.df <- data.frame(do.call(rbind, BCI.ps))
hpd.p.df <- data.frame(do.call(rbind, HPD.ps))

colnames(bci.p.df) <- c("bci_lwr_p", "bci_upr_p")
colnames(hpd.p.df) <- c("hpd_lwr_p", "hpd_upr_p")

complete <- data.frame(target.info, 
                       mean = means,
                       median = medians,
                       var = vars,
                       cv = sqrt(vars)/means,
                       bci.df, hpd.df,
                       bci.p.df, hpd.p.df)

complete

write.csv(complete,
          file = "../saved_data/lognormal_targets.csv",
          row.names = FALSE)