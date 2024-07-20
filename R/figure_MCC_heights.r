library(tidyverse)
library(ggplot2)
library(ggh4x)
source("aux_estimation.r")

compute_height_hpds <- function(obj, level = 0.95, ncores = 4) {
  lookup <- obj$lookup
  only.clades <- obj$clade_log
  K <- ncol(only.clades)
  ints.list <- parallel::mclapply(1:(K-2), function(i)
    get_intervals(X = only.clades[, i],
                  alpha = level), mc.cores = ncores)
  intervals <- do.call(rbind, ints.list)
  out <- tibble::tibble(lookup,
                        intervals)
  return(out)
}

organise_heights <- function(obj) {
  
  lookup <- obj$lookup
  only.clades <- obj$clade_log
  
  rawout <- reshape2::melt(only.clades)
  
  out <- tibble::tibble(clade = rawout$variable,
                        height = rawout$value)
  
  return(list(clade_lookup = lookup,
              melted_heights = out))
}
#########

clade.out <- list(
  lookup = read.csv("../data/clade_lookup_Yang_golden.csv"),
  clade_log = read.csv("../data/combined_golden_Yang.log")
)

clog <- split(clade.out$clade_log, clade.out$clade_log$chain)
## removing first two obs of each chain as "burn-in"
bnin.clog <-  do.call(rbind,
                      lapply(clog, function(x){
                        x[-c(1:2), ]
                      }))

clade.out$clade_log <- bnin.clog


summaries <- compute_height_hpds(obj = clade.out,
                                 level = .95)

heights.df <- organise_heights(obj = clade.out)$melted_heights

heights.df$iteration <- NULL


methods <- rep(c("BCI", "HPD", "Normal_approx"),
               each = nrow(summaries))

interval.properties <- tibble::tibble(
  method = methods,
  lwr = c(summaries$bci_lwr,
          summaries$hpd_lwr,
          summaries$napp_lwr),
  upr = c(summaries$bci_upr,
          summaries$hpd_upr,
          summaries$napp_upr),
  width = c(
    summaries$bci_width,
    summaries$hpd_width,
    summaries$napp_width
  )
)

prop.df <- reshape2::melt(interval.properties)


endpoint.plot <- ggplot(data = prop.df,
                        aes(x = method,
                            y = value, fill = method)) +
  geom_boxplot(alpha = .4) +
  facet_wrap(. ~ variable, scales = "free") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

ggsave(
  plot = endpoint.plot,
  filename = "../figures/Yang_endpoints.pdf",
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)


heights.df$clade <-
  factor(heights.df$clade,
         levels = unique(heights.df$clade) )

colours <- c("normal" = "blue",
             "log-normal" = "gold")


heights.plot <- ggplot(data = heights.df,
                       aes(x = height)) +
  geom_histogram(alpha = .4, aes(y = ..density..)) +
  stat_theodensity(distri = "norm",
                   linewidth = .7,  aes(color = "blue")) +
  stat_theodensity(distri = "lnorm",
                   linewidth = .7,  aes(color = "gold")) +
  facet_wrap(. ~ clade, scales = "free") +
  theme_bw(base_size = 16) +
  scale_x_continuous("Node height", expand = c(0, 0)) + 
  scale_color_identity(
    name = "Distribution",
    breaks = c("blue", "gold"),
    labels = c("normal", "log-normal"),
    guide = "legend"
  ) +
  theme(legend.position = "bottom")


ggsave(
  plot = heights.plot,
  filename = "../figures/Yang_heights.pdf",
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)

interesting.clades <-
  read.csv("../saved_data/clades_of_interest_Yang.csv")

heights.plot.selected <-
  ggplot(data = subset(heights.df,
                       clade %in%
                         interesting.clades$clade),
         aes(x = height)) +
  geom_histogram(alpha = .4, aes(y = ..density..)) +
  stat_theodensity(distri = "norm",
                   linewidth = .7,  aes(color = "blue")) +
  stat_theodensity(distri = "lnorm",
                   linewidth = .7,  aes(color = "gold")) +
  scale_x_continuous("Node height", expand = c(0, 0)) + 
  facet_wrap(. ~ clade, scales = "free") +
  theme_bw(base_size = 16) +
  scale_color_identity(
    name = "Distribution",
    breaks = c("blue", "gold"),
    labels = c("normal", "log-normal"),
    guide = "legend"
  ) +
  theme(legend.position = "bottom")

ggsave(
  plot = heights.plot.selected,
  filename = "../figures/Yang_heights_selected.pdf",
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)
