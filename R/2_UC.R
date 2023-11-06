#' @description Constructs a grid of randomly drawn initial parameters
#' @param nRandom number of parameter vectors to draw
#' @param d.estim boolean. If TRUE, a integration order parameter is drawn. Else not
#' @param det.type type of deterministic If "free", a parameter for the order of the deterministic term is drawn. Else not.
#' @param corr.ind Boolean. If TRUE, correlation parameter is inizialized
#' @return matrix with the initial parameters

thetaRandFree <- function(nRandom, d.estim = TRUE, det.type = "frac", corr.ind = TRUE) {
  set.seed(2)


  Q <- sapply(1:nRandom, function(x) {
    Qmat <- diag(c(runif(1, 0, .7^2), runif(1, 0, .7^2)))
    if (corr.ind) {
      Qmat[1, 2] <- Qmat[2, 1] <- runif(1, -.4, 0) * sqrt(Qmat[1, 1] * Qmat[2, 2]) * as.numeric(corr.ind)
      constrainParam <- mlogvech(Qmat)
      names(constrainParam) <- paste0("vech_", 1:3)
    } else {
      constrainParam <- log(diag(Qmat))
      names(constrainParam) <- c("sd_eta", "sd_eps")
    }
    return(constrainParam)
  }) %>%
    t()

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
    Q,
    AR,
    d_det = D[, 2]
  )

  if (!d.estim) thetaRand <- thetaRand[, -1]
  if (det.type != "free") thetaRand <- thetaRand[, -NCOL(thetaRand)]
  return(thetaRand)
}



#' @param y numerical vector with observations
#' @param corr.ind logical indicator whether correlation should be estimated
#' @param START number of initial residuals to be neglected in the optimization
#' @param d integration order of stoch. trend if it is not estimated. Else FALSE
#' @param rho trend-cycle correlation is it is not estimated. Else FALSE
#' @return Vector holding the ll and parameters

optfn <- function(par0, y, corr.ind = TRUE, START = 10, det.type, d = FALSE, rho = FALSE) {
  # Obtain length of the output
  outputVec.length <- length(par0) + 4
  if (!is.logical(d)) outputVec.length <- outputVec.length + 1
  # Names of the output
  parNames_1 <- c("LL", "d_hat", "Vech1", "Vech2")
  parNames_2 <- c("ar1", "ar2")
  parNames_3 <- c("SdEta", "SdEps")
  if (corr.ind & det.type != "free") {
    parNames <- c(parNames_1, "Vech3", parNames_2, parNames_3, "corr")
  } else if (corr.ind & det.type == "free") {
    parNames <- c(parNames_1, "Vech3", parNames_2, "d.det_hat", parNames_3, "corr")
  } else if (!corr.ind & det.type == "free") {
    parNames <- c(parNames_1, parNames_2, "d.det_hat", parNames_3)
  } else {
    parNames <- c(parNames_1, parNames_2, parNames_3)
  }

  tryCatch(
    {
      est <- optim(
        par = par0, control = list(maxit = 1000),
        # START: how many initial residuals are burned
        fn = fUC_opt_ML_ARMA, d.int = c(0, 2), START = START, pq = c(2, 0),
        nulim = c(0, Inf), corr = corr.ind, nu.opt = FALSE, penalty.corr = TRUE,
        y = y, method = "BFGS", deterministics = det.type, d = d, rho = rho
      )
      if (corr.ind) {
        if (is.logical(d)) {
          Qmat.params <- est$par[2:4]
        } else {
          Qmat.params <- est$par[1:3]
        }
        Qmat <- mlogvech2mat(Qmat.params)
        if (!is.logical(rho)) {
          Qmat[2, 1] <- Qmat[1, 2] <- rho * prod(sqrt(diag(Qmat)))
        }
        corr <- cov2cor(Qmat)[2, 1]

        # Correct -ll to ll
        outputVec <- c(-est$value, est$par, sqrt(Qmat[1, 1]), sqrt(Qmat[2, 2]), corr)
      } else {
        if (is.logical(d)) {
          Qmat.params <- est$par[2:3]
        } else {
          Qmat.params <- est$par[1:2]
        }
        # Correct -ll to ll
        outputVec <- c(-est$value, est$par, sqrt(exp(Qmat.params)))
      }

      # Insert the fixed values into the optimized parameter vector for the subsequent filter run
      if (!is.logical(d)) outputVec <- c(outputVec[1], d, outputVec[-1])
      if (!is.logical(rho)) outputVec <- c(outputVec[1:2], mlogvech(Qmat), outputVec[-c(1:5)])

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



RunUC <- function(y, theta, det.type, corr.ind = NULL) {
  if (is.null(corr.ind)) corr.ind <- ifelse("corr" %in% names(theta), TRUE, FALSE)

  filterOutput <- fUC_opt_ML_ARMA(
    theta = theta, y = y, nulim = c(0, Inf), pq = c(2, 0), corr = corr.ind, d.int = c(0, 2),
    nu.opt = FALSE, penalty.corr = F, deterministics = det.type, return.filter.output = TRUE,
    Q.trans = "mlv"
  )

  outputMat <- cbind(trend = filterOutput$x, cycle = filterOutput$c)
  return(outputMat)
}



UCGridSearch <- function(y, det.type, corr.ind, d = FALSE, nRandom, outputPath, fileName = NULL) {
  d.estim <- ifelse(is.logical(d), FALSE, TRUE)

  thetaMat <- thetaRandFree(nRandom = nRandom, d.estim = d.estim, det.type = det.type, corr.ind = corr.ind)
  UCC_theta <- pbapply(thetaMat, 1, function(theta, y, det.type, corr.ind, d) {
    optfn(theta, y, det.type = det.type, corr.ind = corr.ind, START = 10, d = d)
  }, y = y, det.type = det.type, corr.ind = corr.ind, d) %>%
    t()

  UCC_ResultClean <- UCC_theta[!is.na(UCC_theta[, 1]), ]
  # Sort by -ll in decreasing order
  UCC_Result <- UCC_ResultClean[order(UCC_ResultClean[, 1], decreasing = TRUE), ]

  # Produce a table with the grid search results
  UCC_theta_hat <- UCC_Result[, !c(str_detect(colnames(UCC_Result), paste(c("Sd", "LL", "corr"),
    collapse = "|"
  )))]
  mu_vec <- apply(UCC_theta_hat, 1, function(theta, y, det.type, corr.ind) {
    fUC_opt_ML_ARMA(
      theta = theta, y = y, nulim = c(0, Inf), pq = c(2, 0), corr = corr.ind, d.int = c(0, 2.5),
      nu.opt = FALSE, penalty.corr = F, deterministics = det.type, return.det = TRUE, Q.trans = "mlv"
    )
  }, y = y, det.type = det.type, corr.ind = corr.ind)

  as_tibble(UCC_Result) %>%
    mutate(mu = mu_vec) %>%
    write.xlsx(., file = paste0(outputPath, fileName, ".xlsx"))

  theta <- UCC_theta_hat[1, ]
  return(theta)
}


UC_Wrapper <- function(y, det.type, corr.ind, d = FALSE, nRandom, outputPath, fileName = NULL){
  theta <- UCGridSearch(y = y, det.type = det.type, corr.ind = corr.ind, d = d, nRandom = nRandom, 
               outputPath = outputPath, fileName = fileName)
  output <- RunUC(y = y, theta = theta, det.type = det.type, corr.ind = corr.ind)
  return(output)
}
