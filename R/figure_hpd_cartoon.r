source("aux_true_quantities.r")
source("aux_lognormal_targets.r")
library(ggplot2)

make_cartoon <- function(Mu, Sigmasq, title = "", Alpha = 0.95){
  
  true.BCI <- lognormal_bci(alpha = Alpha,
                            lmean = Mu,
                            lsd = sqrt(Sigmasq))
  true.HPD <- lognormal_hpd(alpha = Alpha,
                            lmean = Mu,
                            lsd = sqrt(Sigmasq))
  
  dlnorm_limit_HPD <- function(x) {
    y <- dlnorm(x, meanlog = Mu, sdlog = sqrt(Sigmasq))
    y[x < true.HPD[1]  |  x > true.HPD[2]] <- NA
    return(y)
  }
  
  dlnorm_limit_BCI <- function(x) {
    y <- dlnorm(x, meanlog = Mu, sdlog = sqrt(Sigmasq))
    y[x < true.BCI[1]  |  x > true.BCI[2]] <- NA
    return(y)
  }
  
  p <- ggplot(data.frame(x = c(0, true.HPD[2] * 2)),
              aes(x = x))
  
  maxf <- dlnorm(x = exp(Mu - Sigmasq), meanlog = Mu, sdlog = sqrt(Sigmasq))
  delta <- maxf/50
  lred <- 1e-02
  lblue <- lred + 2*delta
  
  lims <- range(c(true.BCI, true.HPD))
  mx <- max(lims[1] - .2 * lims[1], 0)
  Mx <- lims[2] + .2 * lims[2]
  
  q <- p +
    # stat_function(
    #   fun = dlnorm_limit_BCI,
    #   geom = "area",
    #   fill = "blue",
    #   alpha = 0.2
    # ) +
    # stat_function(
    #   fun = dlnorm_limit_HPD,
    #   geom = "area",
    #   fill = "red",
    #   alpha = 0.2
    # ) +
    geom_segment(aes(x = true.HPD[1], y = lred,
                     xend = true.HPD[2], yend = lred), 
                 colour = "red") + 
    geom_segment(aes(x = true.HPD[1], y = lred + delta,
                     xend = true.HPD[1], yend = lred - delta), 
                 colour = "red") + 
    geom_segment(aes(x = true.HPD[2], y = lred + delta,
                     xend = true.HPD[2], yend = lred - delta), 
                 colour = "red") + 
    geom_segment(aes(x = true.BCI[1], y = lblue,
                     xend = true.BCI[2], yend = lblue),
                 colour = "blue") + 
    geom_segment(aes(x = true.BCI[1], y = lblue + delta,
                     xend = true.BCI[1], yend = lblue - delta),
                 colour = "blue") + 
    geom_segment(aes(x = true.BCI[2], y = lblue + delta,
                     xend = true.BCI[2], yend = lblue - delta),
                 colour = "blue") + 
    stat_function(fun = function(x) dlnorm(x, Mu, sqrt(Sigmasq)),
                  xlim = c(mx, Mx)) +
    # scale_x_continuous(expression(x), expand = c(0, 0)) +
    # scale_y_continuous(expression(f(x)), expand = c(0, 0)) +
    scale_x_continuous("", expand = c(0, 0)) +
    scale_y_continuous("", expand = c(0, 0)) +
    ggtitle(title) + 
    theme_bw(base_size = 16) 
  
  return(q)
}

S.plot <- make_cartoon(Mu = subset(target.info, target == "symmetric")$m,
                       Sigmasq = subset(target.info, target == "symmetric")$v,
                       title = "Symmetric")

M.plot <- make_cartoon(Mu = subset(target.info, target == "moderate")$m,
                       Sigmasq = subset(target.info, target == "moderate")$v,
                       title = "Moderate")

A.plot <- make_cartoon(Mu = subset(target.info, target == "asymmetric")$m,
                       Sigmasq = subset(target.info, target == "asymmetric")$v,
                       title = "Asymmetric")


library(gridExtra)

big.plot <- grid.arrange(S.plot, M.plot, A.plot,
                         nrow = 1)


ggsave(
  plot = big.plot,
  filename = "../figures/All_lognormal_HPD_cartoon.pdf",
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)
