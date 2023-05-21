source("aux.r")
### Asymmetric : (2.58, 0.254)
Mu <- 2.58
Sigmasq <- 0.254 ^ 2
Alpha <- .95
true.HPD <- lognormal_hpd(alpha = Alpha,
                          lmean = Mu,
                          lsd = sqrt(Sigmasq))

dlnorm_limit <- function(x) {
  y <- dlnorm(x, meanlog = Mu, sdlog = sqrt(Sigmasq))
  y[x < true.HPD[1]  |  x > true.HPD[2]] <- NA
  return(y)
}

library(ggplot2)

p <- ggplot(data.frame(x = c(0, true.HPD[2] * 2)),
            aes(x = x))

q <- p +
  stat_function(
    fun = dlnorm_limit,
    geom = "area",
    fill = "blue",
    alpha = 0.2
  ) +
  stat_function(fun = function(x) dlnorm(x, Mu, sqrt(Sigmasq))) +
  scale_x_continuous(expression(x), expand = c(0, 0)) +
  scale_y_continuous(expression(pi(x)), expand = c(0, 0)) +
  theme_bw(base_size = 16) 

q

ggsave(
  plot = q,
  filename = "../figures/Asymmetric_lognormal_HPD_cartoon.pdf",
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)
