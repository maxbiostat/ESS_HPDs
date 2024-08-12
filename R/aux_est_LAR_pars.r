est_LAR_pars <- function(X, level = 0.95) {
  M <- length(X)
  Y <- log(X)
  fit <- arima(Y, order = c(1, 0, 0))
  coefs <- c(coef(fit), sd(Y))
  ints <- confint(fit)
  x2qs <- qchisq(p = c(1 - level, 1 + level)/2, df = M - 1)
  S2 <- var(Y)
  out <- data.frame(
    parameter = c("phi", "mu", "sigma"),
    point = coefs,
    lwr = c(ints[, 1], sqrt((M - 1)/x2qs[2] * S2)),
    upr = c(ints[,2 ], sqrt((M - 1)/x2qs[1] * S2))
  )
  rownames(out) <- NULL
  return(out)
}
