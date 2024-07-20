compute_all_ess <- function(x, alpha = .95){
  out <- tibble::tibble(
    
    ess = c(coda::effectiveSize(x),
            posterior::ess_basic(x),
            convenience:::essTracerC(x),
            mcmcse::ess(x)
            ),
    
    method = c("coda",
               "Stan",
               "Tracer",
               "mcmcse")
  )
  return(out)
}