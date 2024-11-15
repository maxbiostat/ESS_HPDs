source("aux_LNAR1.r")
source("aux_lognormal.r")

calculate_quantities <- function(m, v, p, xi_p, eff){
  
  Phi <- find_rho(eff = eff, mean = m, variance = v)
  
  return(
    data.frame(
      rho = Phi,
      v_p = compute_LT_variance(p = p,
                                x = xi_p,
                                m = m, v = v,
                                rho = Phi)
    )
  )
  
}

#######################################

target.info.table <- read.csv("../saved_data/lognormal_targets.csv")

Gamma <- 0.95
M <- 1e4
target.esss <- c(50, 100, 200, 500, 625, 1000, 5000)
target.effs <- target.esss/M

J <- length(target.effs)

for.computing  <- reshape2::melt(target.info.table,
                id.vars = c("target", "m", "v",
                            "sigma", "mean", "median", "var", 
                            "bci_lwr", "bci_upr",
                            "hpd_lwr", "hpd_upr"),
                variable.name = "percentile",
                value.name = "p")

the.grid <- do.call(rbind,
                    lapply(target.effs, function(x){
                      data.frame(for.computing,
                                 eff = x)
                    }))

compute.time <- system.time(
  result <- do.call(rbind,
                    parallel::mclapply(1:nrow(the.grid), function(k){
                      ddf <- the.grid[k, ]
                      cname <- gsub("_p", "", ddf$percentile)
                      out <- calculate_quantities(m = ddf$m,
                                                  v = ddf$v,
                                                  p = ddf$p,
                                                  xi_p = ddf[, cname],
                                                  eff = ddf$eff)
                      return(out)
                    }, mc.cores = 10))  
)


compute.time
head(result)

final <- cbind(the.grid, result)

write.csv(final,
          file = "../saved_data/efficiency_table_lognormal.csv",
          row.names = FALSE)
