theo_mrae <- function(n, mu, sigma){
  a_n <- sqrt(2/(pi*n))
  w <- sigma/abs(mu) ## coefficient of variation
  return(
    a_n * w
  )
}

target.ess <- 200
allowable.mrae <- 0.05

f1 <- function(x) theo_mrae(n = target.ess, mu = 1, sigma = x)

curve(f1, 0, 10, lwd = 3,
      ylab = "MRAE", xlab = "Coefficient of variation")
abline(h = allowable.mrae, lwd = 2, lty  = 2)

ESSs <- c(2500, 625, 200, 156.25, 25)

sapply(ESSs, 
       function(n) theo_mrae(n = n, mu = 1, sigma = .25))

max_cv <- function(n, err){
  err/sqrt(2/(pi*n))
} 

sapply(ESSs, 
       function(n) max_cv(n = n, err = 0.05))

theo_mrae(n = ESSs[4], mu = 1, sigma = .78)