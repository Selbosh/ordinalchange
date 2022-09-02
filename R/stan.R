#' Fit an autoregressive ordinal probit model using Stan
#' 
#' This is an alternative to the sequential estimation procedure proposed by Li and Mo.
#' See \url{https://mc-stan.org/docs/stan-users-guide/change-point.html} to better understand the implementation.
#' 
#' Internally we are using \code{optimizing}, i.e. fitting the model using L-BFGS optimization, rather than \code{sampling},
#' which performs Hamiltonian Monte Carlo (automatic differentation). This is because the model is otherwise extremely slow to fit.
#' However, it may be possible to optimize the Stan code such that MCMC becomes more feasible in future.
#' 
#' @param Y Integer vector of ordered categorical responses.
#' @param K Total number of response categories (assumed to go \eqn{1,\dots,K}).
#' @param period The period of the seasonal component of the latent process.
#' @param ... Further arguments passed to [rstan::optimizing()].
#' 
#' @export
ordinalchange_stan <- function(Y, K, period, ...) {
    standata <- list(Y = Y, K = K, Period = Period, T = length(Y))
    out <- rstan::optimizing(stanmodels$aop, data = standata, ...)
    return(out)
}
