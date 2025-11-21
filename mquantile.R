#' M-Quantile Regression via Iteratively Reweighted Least Squares (IRLS)
#'
#' Fits M-quantile regression models using the robust framework of
#' iteratively reweighted least squares (IRLS).  
#' The approach generalizes classical quantile regression by incorporating
#' robustness weights derived from an influence function (e.g., Huber or
#' bisquare), enabling estimation of conditional M-quantiles for heavy-tailed
#' or contaminated data.
#'
#' @details
#' For each quantile level \eqn{q}, the algorithm performs a fixed-point iteration:
#' \deqn{
#'   \hat{\beta}^{(k+1)} =
#'     \arg\min_\beta \sum_{t=1}^n w_t^{(k)} (y_t - x_t^\top \beta)^2,
#' }
#' where the weights \eqn{w_t^{(k)}} depend on an influence function
#' \eqn{\psi}, a specified scale estimator, and the sign of the residuals to
#' encode quantile asymmetry.
#'
#' Supported features include:
#' * Huber or bisquare influence functions  
#' * MAD, Huber, or Proposal–2 scale estimators  
#' * Case weights and variance weights  
#' * M- and MM-estimation  
#' * Multiple quantile levels  
#'
#'
#' @param x A numeric design matrix of predictors.
#' @param y Numeric response vector.
#' @param q A numeric vector of quantile levels in \eqn{(0,1)} (e.g., `0.5` for median).
#' @param psi Influence function used in robustification.
#'   Default: `psi.huber`. Alternative: `psi.bisquare`.
#' @param scale.est Scale estimator to use. One of:
#'   * `"MAD"`  
#'   * `"Huber"`  
#'   * `"proposal 2"`  
#' @param case.weights Observation-level weights. Defaults to equal weights.
#' @param var.weights Variance weights. Defaults to equal weights.
#' @param maxit Maximum number of IRLS iterations (default: `20`).
#' @param acc Convergence tolerance for IRLS updates (default: `1e-4`).
#' @param init Initial estimator: `"ls"` (default), `"lts"`, or a vector/list of coefficients.
#' @param method Estimation method. `"M"` (default) or `"MM"`.
#' @param k2 Tuning constant for Huber-type scale estimators (default: `1.345`).
#' @param test.vec Convergence diagnostic type (`"resid"`, `"coef"`, `"w"`, or `"NULL"`).
#' @param ... Additional arguments passed to the influence function.
#'
#' @return A list with the following elements:
#' \item{fitted.values}{Matrix of fitted values for each quantile level.}
#' \item{residuals}{Matrix of residuals for each quantile level.}
#' \item{q.values}{Quantile levels used in the estimation.}
#' \item{q.weights}{Matrix of final IRLS weights for each quantile.}
#' \item{coefficients}{Matrix of estimated regression coefficients for each quantile level.}
#'
#' @author Patrick Ferreira Patrocinio
#'
#' @import MASS
#'
#' @references
#' Breckling, J. & Chambers, R. (1988). *M-quantiles*. Biometrika, 75, 761–771.
#'
#' @seealso
#' \code{\link[MASS]{rlm}},  
#' \code{\link[MASS]{lqs}},  
#' \code{\link{quantreg}} (quantile regression).
library(MASS)

mqlm <- function(x, y,
                 case.weights = rep(1, nrow(x)),
                 var.weights  = rep(1, nrow(x)),
                 ...,
                 w = rep(1, nrow(x)),
                 init = "ls",
                 psi = psi.huber,
                 scale.est = c("MAD", "Huber", "proposal 2"),
                 k2 = 1.345,
                 method = c("M", "MM"),
                 maxit = 20,
                 acc = 1e-4,
                 test.vec = "resid",
                 q = 0.5) {

  #---------------------------------------------------------------------------
  # Helper functions
  #---------------------------------------------------------------------------
  irls.delta <- function(old, new) {
    sqrt(sum((old - new)^2) / max(1e-20, sum(old^2)))
  }

  irls.rrxwr <- function(x, w, r) {
    w_sqrt <- sqrt(w)
    num <- abs((matrix(r * w_sqrt, 1) %*% x) /
               sqrt(matrix(w_sqrt, 1) %*% (x^2)))
    max(num) / sqrt(sum(w * r^2))
  }

  #---------------------------------------------------------------------------
  # Preprocessing
  #---------------------------------------------------------------------------
  method <- match.arg(method)
  scale.est <- match.arg(scale.est)

  # Ensure matrix form
  x <- as.matrix(x)
  y <- as.numeric(y)

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("X", seq_len(ncol(x)))
  }

  # Rank check
  if (qr(x)$rank < ncol(x)) {
    stop("x is singular: singular fits are not implemented in mqlm()")
  }

  # Validate test.vec
  if (!(test.vec %in% c("resid", "coef", "w", "NULL") || is.null(test.vec))) {
    stop("invalid test.vec")
  }

  # Validate weights
  n <- nrow(x)
  if (length(var.weights) != n) stop("'var.weights' must have length n")
  if (length(case.weights) != n) stop("'case.weights' must have length n")
  if (any(var.weights < 0)) stop("Negative values in var.weights")

  # Combined weights (unchanged formula)
  w <- (w * case.weights) / var.weights

  #---------------------------------------------------------------------------
  # Initial estimation
  #---------------------------------------------------------------------------
  if (method == "M") {

    # Ensure psi is a function
    if (!is.function(psi)) {
      psi <- get(psi, mode = "function")
    }

    # Pass "..." to psi()
    extra_args <- list(...)
    if (length(extra_args)) {
      pm <- pmatch(names(extra_args), names(formals(psi)), nomatch = 0)
      if (any(pm == 0)) warning("some arguments in ... do not match psi() formals")
      args_matched <- names(extra_args)[pm > 0]
      formals(psi)[args_matched] <- unlist(extra_args[args_matched])
    }

    # Initialize
    if (is.character(init)) {
      if (init == "ls") {
        temp <- lm.wfit(x, y, w)
      } else if (init == "lts") {
        temp <- lqs(x, y, intercept = TRUE, nsamp = 200)
      } else stop("unknown init method")

      coef <- temp$coef
      resid <- temp$resid
    } else {
      if (is.list(init)) coef <- init$coef else coef <- init
      resid <- y - x %*% coef
    }

  } else { # MM-estimator -----------------------------------------------------

    temp <- lqs(x, y, intercept = TRUE, method = "S", k0 = 1.548)
    coef <- temp$coef
    resid <- temp$resid
    psi <- psi.bisquare

    # pass c argument
    extra_args <- list(...)
    if (length(extra_args) && "c" %in% names(extra_args)) {
      c0 <- extra_args$c
      if (c0 > 1.548) {
        psi$c <- c0
      } else {
        warning("c must be >= 1.548 and was ignored")
      }
    }

    scale <- temp$scale
  }

  #---------------------------------------------------------------------------
  # IRLS setup
  #---------------------------------------------------------------------------
  done <- FALSE
  converg <- NULL
  n_params <- ncol(x) + 1

  if (scale.est != "MM") {
    scale <- mad(resid / sqrt(var.weights), 0)
  }

  # constants
  theta <- 2 * pnorm(k2) - 1
  gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)

  # Allocate outputs (same structure as original)
  Q <- length(q)
  coef_mat <- matrix(0, n_params, Q)
  w_mat    <- matrix(0, n, Q)
  fit_mat  <- matrix(0, n, Q)
  res_mat  <- matrix(0, n, Q)

  #---------------------------------------------------------------------------
  # Main loop over quantiles
  #---------------------------------------------------------------------------
  for (qi in seq_len(Q)) {

    for (it in seq_len(maxit)) {

      test_prev <- if (!is.null(test.vec)) get(test.vec) else NULL

      # scale update
      if (scale.est != "MM") {
        if (scale.est == "MAD") {
          scale <- median(abs(resid / sqrt(var.weights))) / 0.6745
        } else {
          scale <- sqrt(sum(pmin(resid^2 / var.weights, (k2 * scale)^2)) /
                         ((n - ncol(x)) * gamma))
        }
        if (scale == 0) {
          done <- TRUE
          break
        }
      }

      # compute psi weights
      w <- psi(resid / (scale * sqrt(var.weights))) * case.weights

      # asymmetric quantile weighting
      adj <- ifelse(resid > 0, 2 * q[qi], 2 * (1 - q[qi]))
      w <- w * adj

      # weighted LS step
      temp <- lm.wfit(cbind(1, x), y, w)
      coef <- temp$coef
      resid <- temp$residuals

      # convergence check
      if (!is.null(test.vec)) {
        conv_i <- irls.delta(test_prev, get(test.vec))
      } else {
        conv_i <- irls.rrxwr(x, w, resid)
      }

      converg <- c(converg, conv_i)
      done <- (conv_i <= acc)
      if (done) break
    }

    if (!done) {
      warning("mqlm did not converge at q = ", q[qi])
    }

    coef_mat[, qi] <- coef
    w_mat[, qi]    <- w
    fit_mat[, qi]  <- temp$fitted.values
    res_mat[, qi]  <- resid
  }

  list(
    fitted.values = fit_mat,
    residuals     = res_mat,
    q.values      = q,
    q.weights     = w_mat,
    coefficients  = coef_mat
  )
}

