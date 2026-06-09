#' Fixed weighted generalized covariance measure test
#'
#' Tests against the null hypothesis that X is independent of Y given Z. Allows
#' for some nonlinearity in the dependence relationship, which is accommodated
#' by taking a weighted average of the test statistic across the distribution of
#' Z. This version allows for users to fit their own approximating regressions
#' outside of the wGCM test function. See Scheidegger, Hoerrmann and Buehlmann
#' (2022) "The Weighted Generalised Covariance Measure" for more details, and
#' the weightedGCM package for a fuller implementation that includes model
#' fitting.
#'
#' @param resid.XonZ vector of raw residuals from a regression of X ~ Z
#' @param resid.YonZ vector of raw residuals from a regression of Y ~ Z
#' @param Z vector or matrix of condition variables
#' @param k0 number of quantiles of Z over which to calculate test statistics
#' @param nsim number of simulations used to construct the null distribution
#'
#' @returns A list containing the test statistic (`test.statistic`) and p-value
#' (`p.value`) of the wGCM test.
#' @export
#'
#' @examples
#' #Simulate X --> Z --> Y
#' X <- runif(100, -1, 1)
#' Z <- X^2 + rnorm(100, 0, 0.1)
#' Y <- exp(Z) + rnorm(100, 0, 0.1)
#'
#' #X and Y should be independent conditional on Z
#' residXZ <- residuals(mgcv::gam(X ~ Z))
#' residYZ <- residuals(mgcv::gam(Y ~ Z))
#'
#' test_result <- WGCM_fix(resid.XonZ = residXZ,
#'   resid.YonZ = residYZ,
#'   Z = Z,
#'   k0 = 10,
#'   nsim = 499)
#'
#' test_result
#'
WGCM_fix <- function(resid.XonZ, resid.YonZ, Z, k0, nsim = 499){

  ###
  # Get dimensions of Z
  ###
  Z <- as.matrix(Z)
  n <- nrow(Z)
  dZ <- ncol(Z)

  ###
  # Construct a matrix of weight function values
  ###

  #Starting with w(z) = 1
  W <- rep(1, n)

  #Quantiles given k0 value
  quants <- seq(1, k0)/(k0 + 1)

  #Calculate weights for each quantile for each dimension of Z
  for(d in seq(1, dZ)){

    #Get z and a values
    Z_d <- Z[,d]
    a_values <- stats::quantile(Z_d, quants, names = FALSE)

    #Return value of sign function
    W_d <- outer(Z_d, a_values, function(x, a){
      return(sign(x - a))
    })

    #Add to W
    W <- cbind(W, W_d)
  }

  #Calculate weighted residual product
  R <- resid.XonZ * resid.YonZ * W
  R <- t(R)

  #Normalize
  R_norm <- R / sqrt(rowMeans(R^2) - rowMeans(R)^2)

  #Calculate test stat
  test.statistic <- sqrt(n) * max(abs(rowMeans(R_norm)))

  #Null
  test.stat.null <- apply(abs(R_norm %*% matrix(rnorm(n*nsim), n, nsim)), 2, max) / sqrt(n)

  #Calculate p
  p.value <- (sum(test.stat.null >= test.statistic) + 1)/(nsim + 1)

  #Return
  return(list(test.statistic = test.statistic, p.value = p.value))

}
