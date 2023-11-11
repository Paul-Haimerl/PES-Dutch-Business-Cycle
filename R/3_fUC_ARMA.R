#' @description Generalized (including ARMA cycle) and comprehensive (including both CSS and ML) fUC function
#' @param theta Parameter vector
#' @param y Series to be decomposed into trend and cycle
#' @param nulim Interval for nu
#' @param quiet If FALSE, the value of the objective function is printed out
#' @param approx If TRUE, the steady state Kalman filter will be used once the variance of
#' the prediction error has sufficiently converged
#' @param corr Either TRUE or FALSE. If TRUE, trend and cycle are allowed to be correlated,
#' and nu must specify the 2x2 covariance matrix. If FALSE, trend and cycle are uncorrelated,
#' and nu can be a scalar.
#' @param deterministics Either FALSE, or "const" for constant, "trend" for linear trend,
#' or "frac" for trend of order d
#' @param START First observable y that enters the objective function
#' @param diffuse If TRUE, the Kalman filter will be initialized diffusely
#' @param return.det If true, estimates for the deterministic components are returned
#' @param nu.opt If TRUE, optimization will be conducted setting the variance of
#' the long run innovations to unity. Highly recommended to set this to FALSE, unless
#' the QML estimator is compared to the CSS estimator in a simulation
#' @param d.int Interval for integration order#
#' @param eta ?
#' @param Q.trans Type of transformation for covariance matrix in the estimation. Use
#' of "mlv" is highly recommended
#' @param pq lag order of the ARMA cycle
#' @param m ?
#' @param penalty.corr If TRUE, correlation is bounded away from +/- 1
#' @param return.det.resid If TRUE, estimates for the residuals after controlling for the deterministic
#' component are returned
#' @param return.filter.output If TRUE, the filtered trend (+ deterministics), cycle and pred
#' errors are returned
#' @param d If FALSE, d is treated as a parameter to be estimated. Else, the trend integration order
#' @param rho If FALSE, the correlation is treated as a parameter to be estimated. Else, the
#' correlation parameter
#' @param GLS If TRUE, the GLS correction for the inclusion of deterministics in the ME is used

fUC_opt_ML_ARMA <- function(theta, y, nulim = c(0.05, 30), quiet = TRUE, approx = FALSE,
                            corr = FALSE, deterministics = FALSE, START = 1, diffuse = TRUE,
                            return.det = FALSE, nu.opt = TRUE, d.int = c(0, 2.5),
                            eta = NULL, Q.trans = "mlv", flip = F, neg = F, pq = c(2, 0), m = 30, penalty.corr = TRUE,
                            return.det.resid = FALSE, return.filter.output = FALSE, d = FALSE, rho = FALSE, GLS = FALSE) {
  n <- length(y)
  
  #----------------------------#
  #### Set d                ####
  #----------------------------#
  
  if (is.logical(d)) {
    d <- theta[1]
    theta <- theta[-1]
  }
  
  if (d <= d.int[1] | d >= d.int[2]) {
    return(.Machine$integer.max)
  }
  
  #----------------------------#
  #### Set correlation      ####
  #----------------------------#
  
  if (!corr) {
    # No corr CSS
    if (nu.opt) {
      # Var ratio of trend and cycle innovations
      nu <- exp(theta[1])
      Q <- diag(c(1, nu))
      # Checks
      if (nu[1] > nulim[2] | nu[1] < nulim[1]) {
        return(.Machine$integer.max)
      }
      theta <- theta[-1]
    } else {
      # No corr ML
      sigma_par <- exp(theta[1:2])
      Q <- diag(sigma_par)
      if (sigma_par[1] > nulim[2] | sigma_par[1] < nulim[1]) {
        return(.Machine$integer.max)
      }
      if (sigma_par[2] > nulim[2] | sigma_par[2] < nulim[1]) {
        return(.Machine$integer.max)
      }
      nu <- sigma_par[2] / sigma_par[1]
      theta <- theta[-(1:2)]
    }
    penalty <- 0
  } else {
    # Corr CSS
    if (nu.opt) {
      # Note: Calculation of x and c is invariant to dividing Q by a constant.
      # Thus, Q[1,1] is fixed to unity
      nu_pars <- c(exp(theta[1]), theta[2])
      # check if sigma_eps is bounded
      if (nu_pars[1] < nulim[1] | nu_pars[1] > nulim[2]) {
        return(.Machine$integer.max)
      }
      Q <- matrix(c(1, nu_pars[2], nu_pars[2], nu_pars[1]), 2, 2)
      # In case the correlation is not estimated. Workaround, where the the correlation is still
      # part of the theta vector but simply neglected
      if (!is.logical(rho)) {
        Q[1, 2] <- Q[2, 1] <- rho * prod(sqrt(diag(Q)))
        if (any(eigen(Q)$values <= 0)) {
          return(.Machine$integer.max)
        }
      }
      # check if variance is pd
      corr <- nu_pars[2] / sqrt(nu_pars[1])
      if (corr > 1 | corr < -1) {
        return(.Machine$integer.max)
      }
      nu <- c(Q[1, 1], Q[2, 2], Q[1, 2])
      
      # Check something
      if (!is.null(eta)) {
        ev <- tryCatch(
          {
            eigen(cov2cor(Q))$values
          },
          error = function(e) {
            return(NA)
          }
        )
        if (any(is.na(ev)) | any(ev < 2 * eta)) {
          return(.Machine$integer.max)
        }
      }
      theta <- theta[-(1:2)]
    } else {
      # Corr ML
      Q <- pd_cov(theta[1:3], Q.trans, flip, neg)
      if (any(is.nan(Q))) {
        return(.Machine$integer.max)
      }
      
      # In case the correlation is not estimated. Workaround, where the the correlation is still
      # part of the theta vector but simply neglected
      if (!is.logical(rho)) {
        Q[1, 2] <- Q[2, 1] <- rho * prod(sqrt(diag(Q)))
        if (any(eigen(Q)$values <= 0)) {
          return(.Machine$integer.max)
        }
      }
      nu <- c(Q[1, 1], Q[2, 2], Q[1, 2])
      if (nu[1] > nulim[2] | nu[1] < nulim[1] | nu[2] > nulim[2] | nu[2] < nulim[1]) {
        return(.Machine$integer.max)
      }
      theta <- theta[-(1:3)]
    }
    
    # penalize correlation if too close to +-1
    if (penalty.corr) {
      corr <- cov2cor(Q)[2, 1]
      if (abs(corr) > .99) {
        suppressWarnings(penalty <- -log((1 - abs(corr)) / .01) * 100)
      } else {
        penalty <- 0
      }
    } else {
      penalty <- 0
    }
    if (is.nan(penalty)) penalty <- .Machine$integer.max
  }
  
  
  #----------------------------#
  #### ARMA to AR           ####
  #----------------------------#
  
  
  if (pq[1] > 0) {
    ar <- theta[1:pq[1]]
    theta <- theta[-(1:pq[1])]
  } else {
    ar <- NULL
  }
  if (pq[2] > 0) {
    ma <- theta[1:pq[2]]
    theta <- theta[-(1:pq[2])]
  } else {
    ma <- NULL
  }
  
  # Set the cycle
  if (!is.null(ar)) {
    if (!toComp(-ar)$stable) {
      return(.Machine$integer.max)
    }
    if (!is.null(eta)) {
      if (any(abs(toComp(-ar)$eigv) > (1 - eta))) {
        return(.Machine$integer.max)
      }
    }
  }
  
  if (!is.null(ma)) {
    maInvert <- function(ma) {
      q <- length(ma)
      q0 <- max(which(c(1, ma) != 0)) - 1L
      if (!q0) {
        return(ma)
      }
      roots <- polyroot(c(1, ma[1L:q0]))
      ind <- Mod(roots) < 1
      if (all(!ind)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    if (!maInvert(ma)) {
      return(.Machine$integer.max)
    }
    b <- ARMAtoMA(ar = -ma, ma = ar, lag.max = n - 1)
  } else {
    b <- ARMAtoMA(ma = ar, lag.max = length(ar))
  }
  
  # Trend-cycle part of the UC model (without determ.)
  if (!corr) {
    v <- fUC_comp(y, d, nu, b, corr = FALSE)$v
  } else {
    v <- fUC_comp(y, d, Q, b, corr = TRUE)$v
  }
  resid.tc <- v
  
  
  #----------------------------#
  #### Deterministics       ####
  #----------------------------#
  
  
  # Set the determ.
  if (is.logical(deterministics)) {
    if (deterministics) {
      det.type <- "frac"
    } else {
      det.type <- "none"
      mu <- NA
      v2 <- rep(NA, length(v))
      if (return.det.resid) v <- v2
    }
  } else {
    det.type <- deterministics
    deterministics <- TRUE
  }
  
  # Produce the second vector of determ. residuals and "marry" with first ones via simple OLS reg (GLS theory)
  if (deterministics) {
    if (det.type == "frac") {
      # Omit estimation of the constant. Also no constant in DGP
      Z <- cbind(1, frac_diff(rep(1, n), -d))
      # Z <- as.matrix(frac_diff(rep(1, n), -d))
      if (!corr) {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, nu, b)$v)
      } else {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, Q, b, corr)$v)
      }
    } else if (det.type == "free") {
      d.det <- theta[1]
      if (d.det > d | d - d.det > .5) {
        return(.Machine$integer.max)
      }
      Z <- as.matrix(frac_diff(rep(1, n), -d.det))
      if (!corr) {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, nu, b)$v)
      } else {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, Q, b, corr)$v)
      }
    } else if (det.type == "const") {
      Z <- as.double(matrix(1, n, 1))
      if (!corr) {
        v2 <- fUC_comp(Z, d, nu, b)$v
      } else {
        v2 <- fUC_comp(Z, d, Q, b, corr)$v
      }
    } else if (det.type == "trend") {
      Z <- cbind(1, 1:n)
      # Z <- as.double(matrix(1:n, n, 1))
      if (!corr) {
        v2 <- fUC_comp(Z, d, nu, b)$v
      } else {
        v2 <- fUC_comp(Z, d, Q, b, corr)$v
      }
    } else if (det.type == "trend_corona") {
      Z <- cbind(1:n, 0)
      Z[294, 2] <- 1
      if (!corr) {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, nu, b)$v)
      } else {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, Q, b, corr)$v)
      }
    } else if (det.type == "trend_corona_break") {
      Z <- cbind((1:n), 0, c(rep(0, 215), 1:(n - 215)))
      Z[294, 2] <- 1
      if (!corr) {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, nu, b)$v)
      } else {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, Q, b, corr)$v)
      }
    }
    
    if (!GLS) {
      regmod <- lm(v ~ -1 + v2)
    } else {
      # GLS correction as in Harvey (1990, ch. 3.4.2)
      # gls.correction <- sqrt(diag(Pt)[-(1:(START))])
      gls.correction <- sqrt(diag(Pt)[-1])
      # regmod <- lm(v[-(1:(START - 1))] / gls.correction ~ -1 + I(v2[-(1:(START - 1))] / gls.correction))
      regmod <- lm(v / gls.correction ~ -1 + I(v2 / gls.correction))
    }
    mu <- summary(regmod)$coefficients[, 1]
    v <- residuals(regmod)
    # to harmonize output with other deterministics
    v2 <- as.double(v2)
  }
  
  
  #----------------------------#
  #### CSS or ML Objective  ####
  #----------------------------#
  
  
  if (!any(c(return.det, return.det.resid, return.filter.output))) {
    # In case of the CSS estimator
    if (nu.opt) {
      return(sum(v[START:n]^2))
      # In cse of QML
    } else {
      # Iterations for P of the Kalman filter
      # Construct the system matrices
      Tt <- rbind(
        -frac_diff(c(1, rep(0, n - 1)), d)[-1],
        cbind(diag(n - 2), 0)
      )
      if (!is.null(b)) {
        if (length(b) > m) {
          Tt_c <- rbind(-b[1:m], cbind(diag(m - 1), 0))
        } else {
          if (length(b) > 1) {
            Tt_c <- rbind(-b, cbind(diag(length(b) - 1), 0))
          } else {
            Tt_c <- matrix(-b, 1, 1)
          }
        }
        Tt <- bdiag(Tt, Tt_c)
      } else {
        Tt <- bdiag(Tt, 1)
      }
      Zt <- matrix(c(1, rep(0, n - 2)), ncol = n - 1)
      if (!is.null(b)) {
        Zt <- cbind(Zt, matrix(c(1, rep(0, ncol(Tt_c) - 1)), nrow = 1))
      } else {
        Zt <- cbind(Zt, 1)
      }
      Rt <- matrix(c(1, rep(0, NROW(Tt) - 1)), ncol = 1)
      if (!is.null(b)) {
        Rt <- cbind(Rt, c(rep(0, n - 1), 1, rep(0, ncol(Tt_c) - 1)))
      } else {
        Rt <- cbind(Rt, c(rep(0, n), 1))
      }
      # Initialize
      P_v <- matrix(NA, n, 1)
      if (diffuse) {
        if (length(b) < m) m <- length(b)
        A <- toComp(-b[1:m])$CompMat
        p <- m
        sigma_sq <- Q[2, 2]
        P1A <- matrix(solve(diag(p^2) - (A %x% A)) %*% c(sigma_sq, rep(0, p^2 - 1)), p, p)
        P1Inf <- bdiag(diag(n - 1) * .Machine$integer.max, P1A)
        Pt <- P1Inf
      } else {
        Pt <- matrix(0, NCOL(Tt), NCOL(Tt))
      }
      # Run the recursions
      for (t in 1:n) {
        P_tl <- tcrossprod(tcrossprod(Tt, Pt), Tt) + tcrossprod(tcrossprod(Rt, Q), Rt)
        F_t <- tcrossprod(tcrossprod(Zt, P_tl), Zt)
        Pt <- P_tl - tcrossprod(P_tl, Zt) %*% tcrossprod(Zt, P_tl) / c(F_t)
        P_v[t, ] <- F_t
        # Exploite the steady-state property
        if (approx & t > 1) {
          if (abs(P_v[t, 1] / P_v[t - 1, 1] - 1) < 0.001) {
            P_v[((t + 1):n), 1] <- P_v[t, 1]
            break_t <- t
            break
          }
        }
      }
      # Concentrated log lik according to p. 16 or Harvey ch. 3.4 (eq. 3.4.35) (Note ll already multiplied by -1)
      ll <- 1 / 2 * sum(log(P_v[START:n, ])) + 1 / 2 * sum(v[START:n]^2 / P_v[START:n, ]) + penalty
      return(ll)
    }
  }
  
  
  #----------------------------#
  #### Other output         ####
  #----------------------------#
  
  
  if (return.det) {
    return(mu)
  }
  
  if (return.det.resid) {
    outputList <- list(mu = mu, v_tc = resid.tc, v_det = v2, v = v)
    return(outputList)
  }
  
  if (return.filter.output) {
    # Remove the deterministic trend
    determ <- as.numeric(as.matrix(Z) %*% mu)
    y_tilde <- y - determ
    # Filter again
    if (!corr) {
      filterOutput <- fUC_comp(y_tilde, d, nu, b, corr = FALSE)
    } else {
      filterOutput <- fUC_comp(y_tilde, d, Q, b, corr = TRUE)
    }
    filterOutput$x <- filterOutput$x + determ
    return(filterOutput)
  }
}
