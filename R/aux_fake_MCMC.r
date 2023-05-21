generate_fake_MCMC <- function(N, mu, sigma, phi, LN = FALSE) {
  sigma_e <- sigma * sqrt(1 - phi ^ 2)
  rawSimu <- arima.sim(n = N, list(ar = c(phi)), sd = sigma_e) + mu
  if (LN) {
    exampleSimu <- exp(rawSimu)
  } else{
    exampleSimu <- rawSimu
  }
  return(exampleSimu)
}