library(tidyverse)
library(ggplot2)

target <- "high"
Mu <- 0
Sigma <- 1

MRAE.iid <- read_csv(file = paste0("../saved_data/MRAE_iid_", target, ".csv"))
RMSE.iid <- read_csv(file = paste0("../saved_data/RMSE_iid_", target, ".csv"))
MRAE.MCMC <- read_csv(file = paste0("../saved_data/MRAE_MCMC_", target, ".csv"))
RMSE.MCMC <- read_csv(file = paste0("../saved_data/RMSE_MCMC_", target, ".csv"))  

RMSE.iid <- RMSE.iid %>% rename(target_ESS = M)
RMSE.iid$sampler <- "iid"

RMSE.MCMC$sampler <- "MCMC"

RMSE <- rbind(RMSE.iid, RMSE.MCMC)

MRAE.iid <- MRAE.iid %>% rename(target_ESS = M)
MRAE.iid$sampler <- "iid"

MRAE.MCMC$sampler <- "MCMC"

MRAE <- rbind(MRAE.iid, MRAE.MCMC)


############

ggplot(data = RMSE,
       aes(x = target_ESS, y = rmse, colour = sampler)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous("Effective sample size") + 
  scale_y_continuous("Root mean squared error", expand = c(0, 0)) +
  facet_grid(quantity~., scales = "free_y") + 
  theme_bw(base_size = 16)

source("aux_MRAE.r")
tmrae <- function(target_ESS){
  mm <- exp(Mu + Sigma^2/2)
  vv <- mm^2 *(exp(Sigma^2)-1)  
  ans <- theo_mrae(n = target_ESS,
                   mu = mm,
                   sigma = sqrt(vv))
  return(ans)
}
tmrae <- Vectorize(tmrae)

ggplot(data = MRAE,
       aes(x = target_ESS, y = rae, colour = sampler)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous("Effective sample size") + 
  scale_y_continuous("Mean absolute relative error", expand = c(0, 0)) +
  geom_function(fun = tmrae, inherit.aes = FALSE, linetype = "dotted",
                linewidth = 1) + 
  geom_hline(yintercept = .01, linetype = "longdash") +
  facet_grid(quantity~., scales = "free_y") + 
  theme_bw(base_size = 16)


