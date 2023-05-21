library(tidyverse)
library(ggplot2)
library(ggh4x)
source("aux.r")

compute_height_hpds <- function(fname, level = 0.95) {
  raw.log <- read.table(fname, header = TRUE)
  only.clades <- raw.log[,-1]
  cladenames <- colnames(only.clades)
  cladenames <- gsub("yang.", "", cladenames)
  newNames <- paste0("clade_", match(cladenames, sort(cladenames)))
  lookup <- tibble::tibble(original = cladenames,
                           new = newNames)
  K <- ncol(only.clades)
  ints.list <- parallel::mclapply(1:K, function(i)
    get_intervals(X = only.clades[, i],
                  alpha = level), mc.cores = 4)
  intervals <- do.call(rbind, ints.list)
  out <- tibble::tibble(lookup,
                        intervals)
  return(out)
}

organise_heights <- function(fname) {
  raw.log <- read.table(fname, header = TRUE)
  only.clades <- raw.log[,-1]
  cladenames <- colnames(only.clades)
  cladenames <- gsub("yang.", "", cladenames)
  newNames <- paste0("clade_", match(cladenames, sort(cladenames)))
  lookup <- tibble::tibble(original = cladenames,
                           new = newNames)
  K <- ncol(only.clades)
  rawout <- reshape2::melt(only.clades)
  
  out <- tibble::tibble(clade = lookup$new[match(rawout$variable, sort(unique(rawout$variable)))],
                        height = rawout$value)
  return(list(clade_lookup = lookup,
              melted_heights = out))
}
#########

the.file <- "../data/yang-mcc-clade-heights.log"
summaries <- compute_height_hpds(the.file, level = .99)
heights.df <- organise_heights(the.file)


methods <- rep(c("BCI", "HPD", "Normal_approx"),
               each = nrow(summaries))

interval.properties <- tibble::tibble(
  method = methods,
  lwr = c(summaries$bci_lwr, summaries$hpd_lwr, summaries$napp_lwr),
  upr = c(summaries$bci_upr, summaries$hpd_upr, summaries$napp_upr),
  width = c(
    summaries$bci_width,
    summaries$hpd_width,
    summaries$napp_width
  )
)

prop.df <- reshape2::melt(interval.properties)


endpoint.plot <- ggplot(data = prop.df,
                        aes(x = method, y = value, fill = method)) +
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


heights.df$melted_heights$clade <-
  factor(heights.df$melted_heights$clade,
         levels = paste0("clade_", 1:nrow(heights.df$clade_lookup)))

colours <- c("normal" = "blue",
             "log-normal" = "gold")


heights.plot <- ggplot(data = heights.df$melted_heights,
                       aes(x = height)) +
  geom_histogram(alpha = .4, aes(y = ..density..)) +
  stat_theodensity(distri = "norm",
                   size = .7,  aes(color = "blue")) +
  stat_theodensity(distri = "lnorm",
                   size = .7,  aes(color = "gold")) +
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
  ggplot(data = subset(heights.df$melted_heights,
                       clade %in%
                         interesting.clades$new),
         aes(x = height)) +
  geom_histogram(alpha = .4, aes(y = ..density..)) +
  stat_theodensity(distri = "norm",
                   size = .7,  aes(color = "blue")) +
  stat_theodensity(distri = "lnorm",
                   size = .7,  aes(color = "gold")) +
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