\name{espace}
\alias{espace}
\title{
Enhanced construction of gene regulatory networks
using hub gene information}
\usage{
espace(X,hub_indx, alpha, lambda, maxit_in=1000,maxit_out=5,tol=1e-6)
}
\arguments{
  \item{X}{Observed data matrix with n rows (samples) and
  p columns (variables)}

  \item{hub_indx}{A vector of indices corresponding to the potential hub nodes}

  \item{alpha}{A weight to reduce the penalty on the edges connected to
  the potential hub nodes (0<alpha<=1)}

  \item{lambda}{A tuning parameter imposed on the partial correlations}

  \item{maxit_in}{Maximum number of inner iterations for the partial correlations}
  
  \item{maxit_out}{Maximum number of outer iterations for the diagonal
  elements of the precision matrix}
  

  \item{tol}{Convergence tolerance}

}

\value{
  \item{rho}{Estimated Partial correlation matrix}
  \item{residual}{Residual matrix where the i-th column denotes
  residuals for the multiple regression model with a response X[,i]
  and covariates X[,-i]}
  \item{w_d}{Estimated diagonal elements of the inverse covariance matrix}
} 
\description{
Extended Sparse PArtial Correlation Estimation method}

\details{
Implemented package of the ESPACE method in Yu et al. (2016)
}

\references{
Yu, D., Lim, J., Wang, X., Liang, F. and Xiao, G. (2016).
Enhanced Construction of Gene regulatory Networks using Hub gene
Information. Preprint.
}

\examples{
library(MASS)

## True precision matrix
omega <- matrix(0,5,5)
omega[1,2] <- omega[1,3] <- omega[1,4] <- 1
omega[2,3] <- omega[3,4] <- 1.5
omega <- t(omega) + omega
diag(omega) <- 3

Sig = solve(omega)

X = mvrnorm(50,rep(0,5),Sig)

hub_indx = 1
alpha = 0.7
lambda = 5;
maxit_in = 1000;
maxit_out = 5;
tol = 1e-6;

res <- espace(X,hub_indx, alpha, lambda,maxit_in,maxit_out,tol)

}

