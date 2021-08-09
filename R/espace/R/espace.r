
# SPACE incorporating information of potential hub nodes

espace <- function(X,hub_indx, alpha, lambda,maxit_in=1000,maxit_out=5,tol=1e-6)
{
  n = nrow(X)
  p = ncol(X)
  nh = length(hub_indx)
  rho = matrix(0,p,p)
  rsd = matrix(0,n,p)
  sigma = rep(0,p)
  
  
  out <- .C('espace',n = as.integer(n), p = as.integer(p), nh = as.integer(nh),
             X = as.double(X), hub_indx = as.integer(hub_indx), alpha = as.double(alpha),
             lam = as.double(lambda), niter_in = as.integer(maxit_in), niter_out=as.integer(maxit_out),
             tol = as.double(tol), rho = as.double(rho), residual = as.double(rsd), sigma=as.double(sigma),PACKAGE='espace')
    
  out$rho <- matrix(out$rho,p,p)

  output <- list(rho=out$rho, alpha=alpha, lambda=lambda, residual=matrix(out$residual,n,p), w_d=out$sigma)
  
  return(output)
}                     