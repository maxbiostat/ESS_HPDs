theo_mrae <- function(n, mu, sigma){
  a_n <- sqrt(2/(pi*n))
  w <- sigma/abs(mu) ## coefficient of variation
  return(
    a_n * w
  )
}