### set parameters controlling the model fit
coxaalenic.control <- function(eps = 1e-7, eps.norm = c("max", "grad"),
                               iter.max = 5000, armijo = 1/3, coef.typ = 1,
                               coef.max = 10, trace = FALSE, thread.max = 1)
{
  eps.norm <- match.arg(eps.norm)
  if (eps <= .Machine$double.eps)
    stop("Invalid epsilon. Choose a small value > ", .Machine$double.eps, ".")
  if (!is.element(eps.norm, c("max", "grad")))
    stop(paste("Unknown stopping rule norm", eps.norm))
  if (iter.max < 0)
    stop("Invalid maximum iterations. Choose a large positive integer.")
  if (coef.typ < eps)
    stop("Invalid coefficient magnitude. Choose a positive value.")
  if (coef.max <= coef.typ)
    stop("Invalid maximum coefficient size. Choose a value > ", coef.typ, ".")
  if (armijo < eps | armijo >= 1/2)
    stop("Invalid scale for Armijo's rule. Choose a value in (0, 1/2).")
  if (thread.max < 0)
    stop(paste("Invalid maximum threads. Choose integer in 0, ..., ",
               detectCores(logical = TRUE), ".", sep = ""))
  if (thread.max > detectCores(logical = TRUE)) {
    thread.max <- detectCores(logical = TRUE)
    warning(paste("Invalid maximum threads. Setting 'thread.max' to ",
                  detectCores(logical = TRUE), ".", sep = ""))
  }
  list(eps = eps, eps.norm = eps.norm, iter.max = iter.max, armijo = armijo,
       coef.typ = coef.typ, coef.max = coef.max, trace = trace,
       thread.max = thread.max)
}
