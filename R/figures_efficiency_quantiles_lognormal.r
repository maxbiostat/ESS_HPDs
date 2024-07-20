library(ggplot2)

with.bounds <- read.csv( "../saved_data/ESS_bounds_quantiles_lognormal.csv")

bounds.plot <- 
  ggplot(data = with.bounds,
         aes(x = eff, y = Neff_bound)) +
  geom_line() +
  scale_y_log10("Minimum number of samples") +
  scale_x_continuous("MCMC efficiency") +
  geom_hline(yintercept = 200, linetype = "dotted") +
  geom_hline(yintercept = 625, linetype = "longdash") +
  facet_grid(percentile ~ target) +
  theme_bw(base_size = 16)

bounds.plot
