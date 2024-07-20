## Import table with the specs (m, v, rho, L, U, v_q) for each scenario
### and compute the number of necessary samples to attain a given
### precision, be it absolute or relative.
source("aux_min_ESS.r")

specs <- read.csv(
  "../saved_data/efficiency_table_lognormal.csv"
)

Delta <- 0.05

bounds <- unlist(
  parallel::mclapply(1:nrow(specs),
                     function(k){
                       row <- specs[k, ]
                       ESS_bound_Liu_known_lognormal(m = row$m, v = row$v,
                                                     p = row$p,
                                                     asymp_var = row$v_p.asymp_var,
                                                     delta = Delta)
                     }, mc.cores = 10)
  
)

bounds.iid <- unlist(
  parallel::mclapply(1:nrow(specs),
                     function(k){
                       row <- specs[k, ]
                       ESS_bound_Liu_known_lognormal(m = row$m, v = row$v,
                                                     p = row$p,
                                                     asymp_var = row$v_p.iid_var,
                                                     delta = Delta)
                     }, mc.cores = 10)
  
)


with.bounds <- data.frame(
  specs,
  Neff_bound = bounds,
  Neff_bound_iid =  bounds.iid
)

head(with.bounds)

write.csv(with.bounds, 
          file = "../saved_data/ESS_bounds_quantiles_lognormal.csv",
          row.names = FALSE)
