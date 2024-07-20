target.info <- tibble(
  target = c("symmetric", "moderate", "asymmetric"),
  m = c(4.01, 0.8314, 0),
  v = c(0.0345, 0.30696, 1)^2
)
target.info$sigma <- sqrt(target.info$v)