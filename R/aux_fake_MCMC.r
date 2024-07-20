generate_fake_MCMC <- function(N, m, v, phi, LN = FALSE) {
  sigma <- sqrt(v)
  sigma_e <- sigma * sqrt(1 - phi ^ 2)
  rawSimu <- arima.sim(n = N, list(ar = c(phi)), sd = sigma_e) + m
  if (LN) {
    exampleSimu <- exp(rawSimu)
  } else{
    exampleSimu <- rawSimu
  }
  return(exampleSimu)
}