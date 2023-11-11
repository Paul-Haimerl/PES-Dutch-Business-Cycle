#' @description Constructs a grid of randomly drawn initial parameters
#' @param nRandom number of parameter vectors to draw
#' @param d.estim boolean. If TRUE, a integration order parameter is drawn. Else not
#' @param det.type type of deterministic If "free", a parameter for the order of the deterministic term is drawn. Else not.
#' @param corr.ind Boolean. If TRUE, correlation parameter is inizialized
#' @return matrix with the initial parameters

thetaRandFree <- function(nRandom, d.estim = TRUE, det.type = "frac", corr.ind = TRUE, nu.opt = TRUE) {
  set.seed(2)

  Q <- sapply(1:nRandom, function(x) {
    Qmat <- diag(c(runif(1, 0, .7^2), runif(1, 0, .7^2)))

    if (corr.ind) {
      Qmat[1, 2] <- Qmat[2, 1] <- runif(1, -.4, 0) * sqrt(Qmat[1, 1] * Qmat[2, 2]) * as.numeric(corr.ind)
      if (nu.opt) {
        Q_output <- c(log(Qmat[2, 2] / Qmat[1, 1]), Qmat[1, 2] / Qmat[1, 1])
        names(Q_output) <- paste0("nu_", 1:2)
      } else {
        Q_output <- mlogvech(Qmat)
        names(Q_output) <- paste0("vech_", 1:3)
      }
    } else {
      if (nu.opt) {
        Q_output <- as.vector(log(Qmat[2, 2] / Qmat[1, 1]))
        names(Q_output) <- "nu"
      } else {
        Q_output <- log(diag(Qmat))
        names(Q_output) <- c("sd_eta", "sd_eps")
      }
    }
    return(Q_output)
  }) %>%
    t()
  if (nu.opt & !corr.ind) {
    Q <- matrix(Q, nc = 1)
    rownames(Q) <- NULL
    colnames(Q) <- "nu"
  }

  AR <- sapply(1:nRandom, function(x) {
    ar <- rep(NA, 2)
    ar[1] <- runif(1, -1.3, -.7)
    ar[2] <- runif(1, .3, .9)
    while (!toComp(-ar)$stable) {
      ar[1] <- runif(1, -1.3, -.7)
      ar[2] <- runif(1, .1, .7)
    }
    names(ar) <- paste0("ar_", 1:2)
    return(ar)
  }) %>%
    t()

  D <- sapply(1:nRandom, function(x) {
    D <- rep(NA, 2)
    D[1] <- runif(1, 1, 1.5)
    d_diff <- runif(1, 0, .5)
    D[2] <- D[1] - d_diff
    names(D) <- c("d", "d_det")
    return(D)
  }) %>%
    t()

  thetaRand <- cbind(
    d = D[, 1],
    as.matrix(Q),
    AR,
    d_det = D[, 2]
  )

  if (!d.estim) thetaRand <- thetaRand[, -1]
  if (det.type != "free") thetaRand <- thetaRand[, -NCOL(thetaRand)]
  return(thetaRand)
}


#' @description Executing the numerical (quasi-)ML optimization
#' @param par0 initial parameter vector
#' @param y numerical vector with observations
#' @param corr.ind logical indicator whether correlation should be estimated
#' @param START number of initial residuals to be neglected in the optimization
#' @param d integration order of stoch. trend if it is not estimated. Else FALSE
#' @param rho trend-cycle correlation is it is not estimated. Else FALSE
#' @return Vector holding the ll and parameters

optfn <- function(par0, y, corr.ind = TRUE, START = 10, det.type, d = FALSE, rho = FALSE, nu.opt = TRUE) {
  # Obtain length of the output
  outputVec.length <- length(par0) + 4
  if (!is.logical(d)) outputVec.length <- outputVec.length + 1
  # Names of the output
  if (corr.ind) {
    if (nu.opt) {
      Qmat.param.names <- paste0("nu_", 1:2)
    } else {
      Qmat.param.names <- paste0("Vech_", 1:3)
    }
  } else {
    if (nu.opt) {
      Qmat.param.names <- "nu"
    } else {
      Qmat.param.names <- c("var_eta", "var_eps")
    }
  }
  # maybe extend with for diverse pq values
  parNames <- c("ll", "d_hat", Qmat.param.names, "ar_1", "ar_2", "sd_eta", "sd_eps", "corr")
  if (det.type == "free") parNames <- c(parNames[1:(length(Qmat.param.names) + 4)], "d.det_hat", parNames[-(1:(length(Qmat.param.names) + 4))])

  # Perform the actual optimization
  tryCatch(
    {
      est <- optim(
        par = par0, control = list(maxit = 1000),
        # START: how many initial residuals are burned
        fn = fUC_opt_ML_ARMA, d.int = c(0, 2), START = START, pq = c(2, 0),
        corr = corr.ind, nu.opt = nu.opt, penalty.corr = TRUE,
        y = y, method = "BFGS", deterministics = det.type, d = d, rho = rho
      )

      # Select and transform the var and corr parameters according to the respective specification
      Qmat.params.ind <- 2:(length(Qmat.param.names) + 1)
      if (!is.logical(d)) Qmat.params.ind <- Qmat.params.ind - 1
      Qmat.params <- est$par[Qmat.params.ind]
      Qmat <- param2Qmat(Qmat.params = Qmat.params, corr.ind = corr.ind, nu.opt = nu.opt)
      if (!is.logical(rho)) {
        Qmat[2, 1] <- Qmat[1, 2] <- rho * prod(sqrt(diag(Qmat)))
      }
      corr <- cov2cor(Qmat)[2, 1]

      # Correct -ll to ll and put the rest of the output together
      outputVec <- c(-est$value, est$par, sqrt(Qmat[1, 1]), sqrt(Qmat[2, 2]), corr)
      # Insert the fixed d or rho values into the optimized parameter vector for the subsequent filter run
      if (!is.logical(d)) outputVec <- c(outputVec[1], d, outputVec[-1])
      if (corr.ind & !is.logical(rho)) {
        Qmat.param <- Qmat.transform(Q = Qmat, corr.ind = corr.ind, nu.opt = nu.opt)
        outputVec <- c(outputVec[1:2], Qmat.param, outputVec[-(1:(2 + length(Qmat.params)))])
      }
      names(outputVec) <- parNames
      return(outputVec)
    },
    error = function(e) {
      outputVec_NA <- rep(NA, outputVec.length)
      names(outputVec_NA) <- parNames
      return(outputVec_NA)
    }
  )
}


Qmat.transform <- function(Q, corr.ind, nu.opt) {
  if (corr.ind) {
    if (nu.opt) {
      Q.param <- c(log(Q[2, 2] / Q[1, 1]), Q[2, 1] / Q[1, 1])
    } else {
      Q.param <- mlogvech(Q)
    }
  } else {
    if (nu.opt) {
      Q.param <- log(Q[2, 2] / Q[1, 1])
    } else {
      Q.param <- log(diag(Q))
    }
  }
  return(Q.param)
}

param2Qmat <- function(Qmat.params, corr.ind, nu.opt) {
  if (corr.ind) {
    if (nu.opt) {
      Qmat <- diag(c(1, exp(Qmat.params[1])))
      Qmat[1, 2] <- Qmat[2, 1] <- Qmat.params[2]
    } else {
      Qmat <- mlogvech2mat(Qmat.params)
    }
  } else {
    if (nu.opt) {
      Qmat <- diag(c(1, exp(Qmat.params)))
    } else {
      Qmat <- diag(c(exp(Qmat.params)))
    }
  }
  return(Qmat)
}



RunUC <- function(y, theta, det.type, corr.ind = NULL, nu.opt) {
  if (is.null(corr.ind)) corr.ind <- ifelse("corr" %in% names(theta), TRUE, FALSE)

  filterOutput <- fUC_opt_ML_ARMA(
    theta = theta, y = y, nulim = c(0, Inf), pq = c(2, 0), corr = corr.ind, d.int = c(0, 2),
    nu.opt = nu.opt, penalty.corr = F, deterministics = det.type, return.filter.output = TRUE,
    Q.trans = "mlv"
  )

  outputMat <- cbind(trend = filterOutput$x, cycle = filterOutput$c)
  return(outputMat)
}



UCGridSearch <- function(y, det.type, corr.ind, d = FALSE, nRandom, outputPath, fileName = NULL, nu.opt) {
  d.estim <- ifelse(is.logical(d), TRUE, FALSE)

  thetaMat <- thetaRandFree(nRandom = nRandom, d.estim = d.estim, det.type = det.type, corr.ind = corr.ind, nu.opt = nu.opt)
  UCC_theta <- pbapply(thetaMat, 1, function(theta, y, det.type, corr.ind, d, nu.opt) {
    optfn(theta, y, det.type = det.type, corr.ind = corr.ind, START = 10, d = d, nu.opt = nu.opt)
  }, y = y, det.type = det.type, corr.ind = corr.ind, d = d, nu.opt = nu.opt) %>%
    t()
  UCC_ResultClean <- UCC_theta[!is.na(UCC_theta[, 1]), ]
  # Sort by -ll in decreasing order
  UCC_Result <- UCC_ResultClean[order(UCC_ResultClean[, 1], decreasing = TRUE), ]
  # Produce a table with the grid search results
  UCC_theta_hat <- UCC_Result[, !c(str_detect(colnames(UCC_Result), paste(c("sd", "ll", "corr"),
    collapse = "|"
  )))]
  muVec <- apply(UCC_theta_hat, 1, function(theta, y, det.type, corr.ind, nu.opt) {
    fUC_opt_ML_ARMA(
      theta = theta, y = y, nulim = c(0, Inf), pq = c(2, 0), corr = corr.ind, d.int = c(0, 2.5),
      nu.opt = nu.opt, penalty.corr = F, deterministics = det.type, return.det = TRUE, Q.trans = "mlv"
    )
  }, y = y, det.type = det.type, corr.ind = corr.ind, nu.opt = nu.opt)

  UCC_Result_tib <- as_tibble(UCC_Result)
  if (is.matrix(muVec)) {
    rownames(muVec) <- paste0("mu_", 1:NROW(muVec))
    UCC_Result_tib <- cbind(UCC_Result_tib, t(muVec)) %>%
      as_tibble()
  } else {
    UCC_Result_tib$mu <- muVec
  }
  write.xlsx(UCC_Result_tib, file = paste0(outputPath, fileName, ".xlsx"))

  theta <- UCC_theta_hat[1, ]
  return(theta)
}


UC_Wrapper <- function(y, det.type, corr.ind, d = FALSE, nRandom, outputPath, fileName = NULL, theta = NULL) {
  CSS <- TRUE
  if (is.null(theta)) {
    theta <- UCGridSearch(
      y = y, det.type = det.type, corr.ind = corr.ind, d = d, nRandom = nRandom,
      outputPath = outputPath, fileName = fileName, nu.opt = CSS
    )
  }
  browser()
  output <- RunUC(y = y, theta = theta, det.type = det.type, corr.ind = corr.ind, nu.opt = CSS)
  return(output)
}
