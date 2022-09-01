#' Estimate the AOP(1) model from Li and Lu (2022)
#' 
#' The parameter vector is \eqn{\theta = (c_2, \dots, c_{K-1}, \alpha_0, \alpha_1, A, D, \Delta)}.
#' Starting values for \eqn{\alpha_1, A, D, \Delta} are set to 0.
#' Starting values for \eqn{\alpha_0} and thresholds \eqn{c_2,\dots,c_{K-1}} are \eqn{\alpha_0 = -\Phi^s that \eqn{Y_t} lies in category \eqn{j}.
#' {-1}(p_1)} and \eqn{c_k=-\Phi^{-1}(\sum_{j=1}^k p_j) + \alpha_0}, where \eqn{p_j} is the sample proportion of time
#' @param Y Integer vector.
#' @param K Integer number of response categories, denoted \eqn{k = 1, \dots, K}.
#' @param iterations Number of iterations for Newton's method.
#' 
#' @source \url{https://doi.org/10.1002/env.2752}
#' 
#' @examples 
#' Y <- sample(1:3, 50, replace = T)
#' seq_estimation(Y)
#' 
#' @export
seq_estimation <- function(Y, K = max(Y), tau = 10, L = 7, ...) {
    # Input checking.
    stopifnot(all(Y > 0))
    stopifnot(K > 1)
    stopifnot(tau > 1 && tau < length(Y))
    stopifnot(L > 0)
    
    # Initialization.
    cutpoint <- vector(length = K - 2)
    alpha0 <- -qnorm(mean(Y == 1))
    for (k in 2:(K - 1))
      cutpoint[k - 1] <- qnorm(mean(Y <= k)) - alpha0
    alpha1 <- A <- D <- Delta <- 0
    init_theta <- c(cutpoint, alpha0, alpha1, A, D, Delta)


    # Newton's method.
    result <- nlm(Q, init_theta, Y = Y, tau = tau, L = L)

    # Method of moments.
    theta_hat <- result$estimate
    mu <- mu(alpha0 = theta_hat[K-1],
             alpha1 = theta_hat[K],
             A = theta_hat[K + 1],
             D = theta_hat[K + 2],
             Delta = theta_hat[K + 3],
             Y, tau, L)
    rho <- estimate_rho(theta_hat, Y, K, mu)
    return(list(theta = theta_hat, rho = rho))
}

mu <- function(alpha0, alpha1, A, D, Delta, Y, tau, L) {
    t <- seq_along(Y)
    n <- length(Y)
    s <- A * cos(2 * pi * t / L) + D * sin(2 * pi * t / L)
    alpha0 + alpha1 * t / n + s + Delta * (t >= tau)
}

#' @import mvtnorm
estimate_rho <- function(theta, Y, K = max(Y), mu) {
    init <- acf(Y, lag.max = 1, plot = FALSE)$acf[2]
    m <- function(rho) {
        E_YY <- function(t, k, rho) { # E[Y_t,k Y_{t+1},k]
            threshold <- c(-Inf, 0, theta[1:(K - 2)], Inf)
            lower <- threshold[k] # indexed from c_0
            upper <- threshold[k + 1]
            mvtnorm::pmvnorm(lower = rep(lower, 2), upper = rep(upper, 2),
                             mean = mu[t:(t + 1)],
                             sigma = matrix(c(1, rho, rho, 1), 2))
        }
        m <- 0
        n <- length(Y)
        for (t in 1:(n - 1))
          for (k in 1:K)
            m <- m + E_YY(t, k, rho) # careful!
        return(m)
    }

    Y_b <- outer(Y, 1:K, '==') # n x K
    m_hat <- sum(head(Y_b, -1) * tail(Y_b, -1))

    optim(init, function(x) abs(m_hat - m(x)),
          method = 'L-BFGS-B', lower = -abs(init), upper = abs(init))
}

#' Sum of squared discrepancies between observations and their expected values
#' 
#' @param theta A vector of parameters.
#' @param Y A vector of ordered categorical responses.
#' @param tau A candidate changepoint.
#' @param L The known period of the seasonal component.
#' 
#' @note The symbol \code{T} has been replaced with \code{L} to avoid confusion with \code{T == TRUE}.
Q <- function(theta, Y, tau, L) {
    # Unpack the parameters
    K <- length(theta) - 5 + 2
    cutpoints <- c(0, head(theta, K - 2), +Inf)
    alpha0 <- theta[K - 1]
    alpha1 <- theta[K]
    A <- theta[K + 1]
    D <- theta[K + 2]
    Delta <- theta[K + 3]

    # Compute level probabilities and expected values.
    n <- length(Y)
    mu <- mu(alpha0, alpha1, A, D, Delta, Y, tau, L)
    E_Y <- xi <- omega <- matrix(nrow = n, ncol = K)
    E_Y[, 1] <- pnorm(0 - mu) - 0 # c0 = -Inf
    xi[, 1] <- dnorm(0 - mu) - 0  # c0 = -Inf
    omega[, 1] <- (Y == 1) - E_Y[, 1]
    for (k in 2:K) {
      E_Y[, k] <- pnorm(cutpoints[k] - mu) - pnorm(cutpoints[k - 1] - mu)
      xi[, k] <- dnorm(cutpoints[k] - mu) - dnorm(cutpoints[k - 1] - mu)
      omega[, k] <- (Y == k) - E_Y[, k]
    }

    # Compute the sum of squared differences.
    result <- sum((outer(Y, 1:K, '==') - E_Y)^2)
    # Compute the gradient.
    attr(result, 'gradient') <- gradient_theta(theta, omega, xi, mu, tau, L)
    return(result)
}

gradient_theta <- function(theta, omega, xi, mu, tau, L) {
    # Unpack the parameters
    K <- length(theta) - 5 + 2
    cutpoints <- c(0, head(theta, K - 2))
    alpha0 <- theta[K - 1]
    alpha0 <- theta[K]
    A <- theta[K + 1]
    D <- theta[K + 2]
    Delta <- theta[K + 3]
    gradient(cutpoints, alpha0, alpha1, A, D, Delta, omega, xi, mu, tau, L)
}

gradient <- function(cutpoints, alpha0, alpha1, A, D, Delta, omega, xi, mu, tau, L) {
    K <- length(cutpoints) + 1
    d_cutpoints <- vector(length = K - 1) # d_c1 is redundant and left NA
    for (k in 2:(K - 1)) {
        d_cutpoints[k] <- 2 * sum((omega[, k + 1] - omega[, k]) * dnorm(cutpoints[k] - mu))
    }

    t <- seq_along(Y)
    n <- length(Y)
    d_alpha0 <- 2 * sum(omega * xi * 1)
    d_alpha1 <- 2 * sum(omega * xi * t / n)
    d_A <- 2 * sum(omega * xi * cos(2 * pi / L * t))
    d_D <- 2 * sum(omega * xi * sin(2 * pi / L * t))
    d_Delta <- 2 * sum (omega * xi * (t >= tau))

    c(d_cutpoints[-1], d_alpha0, d_alpha1, d_A, d_D, d_Delta)
}

# hessian <- function(theta) {
# }