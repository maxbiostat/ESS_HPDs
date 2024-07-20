### Here we are interested in the question:
### If I draw n iid variates X = (X_1, ..., X_n) from a distribution with cdf F
### What is the probability that min(X) <= a and max(X) > b, for a < b two known
### constants?
### Turns out that the answer is
### Pr(min(X) > a and max(X) <= b) =  1 - {F(b)^n - (F(b)-F(a))^n}
Rcpp::cppFunction(
  'double logDiffExp(double x, double y)
{return x > y ?
Rf_logspace_sub(x, y) :
Rf_logspace_sub(y, x);}'
)
prob_of_interest <- function(logcdf, n, a, b) {
  if (b <= a)
    stop("b needs to be greater than a")
  lFa <- logcdf(a)
  lFb <- logcdf(b)
  if (lFb < lFa)
    stop("logcdf needs to be non-decreasing")
  lPA <- logDiffExp(0, n * logDiffExp(0, lFa))
  lPA_and_B <- logDiffExp(n * lFb, n * logDiffExp(lFb, lFa))
  lprob <- logDiffExp(lPA, lPA_and_B)
  return(exp(lprob))
}