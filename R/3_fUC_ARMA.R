fUC_opt_ML_ARMA <- function(theta, y, nulim = c(0.05, 30), quiet = TRUE, approx = TRUE,
                            corr = FALSE, deterministics = FALSE, START = 1, diffuse = FALSE,
                            return.det = FALSE, nu.opt = TRUE, d.int = c(0, 2.5),
                            eta = NULL, Q.trans = "mlv", flip = F, neg = F, pq = c(2, 0), m = 30, penalty.corr = TRUE,
                            return.det.resid = FALSE, return.filter.output = FALSE, d = FALSE, rho = FALSE, GLS = FALSE) {
  # To allow for d to be held fixed in the optimization
  if (is.logical(d)) {
    d <- theta[1]
  } else {
    theta <- c(0, theta)
  }

  n <- length(y)


  #----------------------------#
  #### Set correlation      ####
  #----------------------------#


  # NO CORRELATION
  if (!corr) {
    # here even applicable? Level from Q defines Var of one-step-ahead error v??
    if (nu.opt) {
      nu <- exp(theta[2])
      Q <- diag(c(1, nu))
      if (nu[1] > nulim[2] | nu[1] < nulim[1]) {
        return(.Machine$integer.max)
      }
    } else {
      sigma_par <- exp(theta[2:3])
      Q <- diag(sigma_par)
      theta <- theta[-2]
      if (sigma_par[1] > nulim[2] | sigma_par[1] < nulim[1]) {
        return(.Machine$integer.max)
      }
      if (sigma_par[2] > nulim[2] | sigma_par[2] < nulim[1]) {
        return(.Machine$integer.max)
      }
      nu <- sigma_par[2] / sigma_par[1]
    }

    Ht <- matrix(0, 1, 1)
    if (length(theta) > 2) {
      if (pq[1] > 0) {
        ar <- theta[3:(2 + pq[1])]
      } else {
        ar <- NULL
      }

      if (pq[2] > 0) {
        ma <- theta[(3 + pq[1]):(2 + pq[1] + pq[2])]
      } else {
        ma <- NULL
      }
    }
    penalty <- 0
  } else {
    # CORRELATION

    Q <- pd_cov(theta[2:4], Q.trans, flip, neg)
    if (any(is.nan(Q)) | any(is.infinite(Q))) {
      return(.Machine$integer.max)
    }
    if (!is.logical(rho)) {
      Q[1, 2] <- Q[2, 1] <- rho * prod(sqrt(diag(Q)))
      if (any(eigen(Q)$values <= 0)) {
        return(.Machine$integer.max)
      }
    }

    # if(!is.null(eta)){
    #    if(any(eigen(cov2cor(Q))$values< 2*eta)) return(.Machine$integer.max)
    # }
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

    # Q <- trans2mat(theta[2:4])
    Ht <- matrix(0, 1, 1)
    nu <- c(Q[1, 1], Q[2, 2], Q[1, 2])
    if (length(theta) > 4) {
      if (pq[1] > 0) {
        ar <- theta[5:(4 + pq[1])]
      } else {
        ar <- NULL
      }
      if (pq[2] > 0) {
        ma <- theta[(5 + pq[1]):(4 + pq[1] + pq[2])]
      } else {
        ma <- NULL
      }
    } else {
      ar <- NULL
      ma <- NULL
    }
    if (nu[1] > nulim[2] | nu[1] < nulim[1] | nu[2] > nulim[2] | nu[2] < nulim[1]) {
      return(.Machine$integer.max)
    }
  }


  #----------------------------#
  #### ARMA to AR           ####
  #----------------------------#


  # to keep the integration order bounded
  if (d <= d.int[1] | d >= d.int[2]) {
    return(.Machine$integer.max)
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


  #---------------------------------------------------------------------------------------------------#
  #### Kalman filter-ish routine to obtain Var-Cov of Frankenstein one-step-ahead prediction error ####
  #---------------------------------------------------------------------------------------------------#


  if (!any(c(return.det, return.det.resid, return.filter.output))) {
    # iterations for P of the Kalman filter:
    Tt <- rbind(
      -frac_diff(c(1, rep(0, length(y) - 1)), d)[-1],
      cbind(diag(length(y) - 2), 0)
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

    Zt <- matrix(c(1, rep(0, length(y) - 2)), ncol = length(y) - 1)
    if (!is.null(b)) {
      Zt <- cbind(Zt, matrix(c(1, rep(0, ncol(Tt_c) - 1)), nrow = 1))
    } else {
      Zt <- cbind(Zt, 1)
    }

    Rt <- matrix(c(1, rep(0, NROW(Tt) - 1)), ncol = 1)
    if (!is.null(b)) {
      Rt <- cbind(Rt, c(rep(0, length(y) - 1), 1, rep(0, ncol(Tt_c) - 1)))
    } else {
      Rt <- cbind(Rt, c(rep(0, length(y)), 1))
    }

    P_v <- matrix(NA, length(y), 1)

    diffuse <- TRUE
    approx <- FALSE

    # initialization:
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
      #     Pt[1,1] <- Q[1,1]
      #     Pt[length(y), length(y)] <- Q[2,2]
    }

    # Error here if nu.opt == TRUE?

    for (t in 1:n) {
      P_tl <- tcrossprod(tcrossprod(Tt, Pt), Tt) + tcrossprod(tcrossprod(Rt, Q), Rt)
      F_t <- tcrossprod(tcrossprod(Zt, P_tl), Zt)
      Pt <- P_tl - tcrossprod(P_tl, Zt) %*% tcrossprod(Zt, P_tl) / c(F_t)
      P_v[t, ] <- F_t

      # I.e. exploiting the steady-state property
      if (approx & t > 1) {
        if (abs(P_v[t, 1] / P_v[t - 1, 1] - 1) < 0.001) {
          P_v[((t + 1):n), 1] <- P_v[t, 1]
          break_t <- t
          break
        }
      }
    }
  }

  #----------------------------#
  #### Deterministics       ####
  #----------------------------#


  # account for deterministic terms
  # default is "no deterministics". (FALSE)
  # prio2 is "deterministic trend of order d" (TRUE)
  # prio3 is "constant" (constant)
  # prio4 is "constant + linear trend" (trend)

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
      # Z <- cbind(1, frac_diff(rep(1, length(y)), -d))
      Z <- as.matrix(frac_diff(rep(1, length(y)), -d))
      if (!corr) {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, nu, b)$v)
      } else {
        v2 <- apply(Z, 2, function(x) fUC_comp(x, d, Q, b, corr)$v)
      }
    } else if (det.type == "free") {
      # Omit estimation of the constant. Also no constant in DGP
      # Z <- cbind(1, frac_diff(rep(1, length(y)), -d))
      d.det <- theta[length(theta)]
      if (d.det > d | d - d.det > .5) {
        return(.Machine$integer.max)
      }
      Z <- as.matrix(frac_diff(rep(1, length(y)), -d.det))
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
      Z <- as.double(matrix(1:n, n, 1))
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
      # gls.correction <- sqrt(diag(Pt)[-(1:(START))])
      gls.correction <- sqrt(diag(Pt)[-1])
      # regmod <- lm(v[-(1:(START - 1))] / gls.correction ~ -1 + I(v2[-(1:(START - 1))] / gls.correction))
      regmod <- lm(v / gls.correction ~ -1 + I(v2 / gls.correction))
    }
    mu <- summary(regmod)$coefficients[, 1]
    v <- residuals(regmod)
    # to harmonize output with other dets
    v2 <- as.double(v2)
  }

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


  #----------------------------#
  #### Calculate likelihood ####
  #----------------------------#


  # See formula on p. 16 or Harvey ch. 3.4 (eq. 3.4.35) for concentrated log lik (Note already multiplied by -1)
  ll <- 1 / 2 * sum(log(P_v[START:n, ])) + 1 / 2 * sum(v[START:n]^2 / P_v[START:n, ]) + penalty
  return(ll)
}
