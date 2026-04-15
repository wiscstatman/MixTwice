mixtwice_nloptr <- function(thetaHat, s2, Btheta, Bsigma2, df) {
  
  theta0 <- thetaHat
  s20 <- s2
  p0 <- length(theta0)
  
  theta <- theta0
  s2 <- s20
  p <- length(theta)

  cc <- max(abs(theta)) * 1.1
  if (!is.finite(cc) || cc <= 0) cc <- 1
  
  gridtheta <- (cc / Btheta) * seq(-Btheta, Btheta, by = 1)
  
  sd_min <- sqrt(min(s2))
  sd_max <- sqrt(max(s2))
  if (sd_max == sd_min) {
    gridsigma <- rep(sd_min, Bsigma2 + 1)
  } else {
    gridsigma <- seq(sd_min, sd_max, by = (sd_max - sd_min) / Bsigma2)
  }
  
  ltheta <- length(gridtheta)
  lsigma <- length(gridsigma)
  
  grid1 <- rep(gridtheta, each = lsigma)
  grid2 <- rep(gridsigma, ltheta)
  

  likelihood.theta <- matrix(NA_real_, nrow = p, ncol = length(grid1))
  likelihood.s2 <- matrix(NA_real_, nrow = p, ncol = length(grid2))
  
  for (i in seq_len(p)) {
    likelihood.theta[i, ] <- stats::dnorm(theta[i], mean = grid1, sd = grid2)
    likelihood.s2[i, ] <- stats::dchisq(df * s2[i] / grid2^2, df = df) * df / grid2^2
  }
  
  likelihood.conditional <- likelihood.theta * likelihood.s2
  

  L <- function(x) {
    xtheta <- x[1:ltheta]
    xsigma <- x[(ltheta + 1):(ltheta + lsigma)]
    
    yy <- array(xsigma %o% xtheta)
    dens <- pmax(as.vector(yy %*% t(likelihood.conditional)), 1e-300)
    -sum(log(dens))
  }
  
  G <- function(x) {
    g <- numeric(ltheta)
    h <- numeric(lsigma)
    
    xtheta <- x[1:ltheta]
    xsigma <- x[(ltheta + 1):(ltheta + lsigma)]
    
    yy <- array(xsigma %o% xtheta)
    d <- pmax(as.vector(yy %*% t(likelihood.conditional)), 1e-300)
    
    for (i in seq_len(ltheta)) {
      cols <- (lsigma * (i - 1) + 1):(lsigma * i)
      g[i] <- -sum((xsigma %*% t(likelihood.conditional[, cols, drop = FALSE])) / d)
    }
    
    for (j in seq_len(lsigma)) {
      cols <- seq(j, j + (ltheta - 1) * lsigma, by = lsigma)
      h[j] <- -sum((xtheta %*% t(likelihood.conditional[, cols, drop = FALSE])) / d)
    }
    
    c(g, h)
  }
  

  idx0 <- Btheta + 1
  
  hin <- function(x) {
    g <- x[1:ltheta]
    c(
      g[1:Btheta] - g[2:(Btheta + 1)],                 # increasing to mode
      -g[(Btheta + 1):(ltheta - 1)] + g[(Btheta + 2):ltheta]  # decreasing after mode
    )
  }
  
  hinjac <- function(x) {
    J <- matrix(0, nrow = ltheta - 1, ncol = ltheta + lsigma)
    
    for (i in seq_len(Btheta)) {
      J[i, i] <- 1
      J[i, i + 1] <- -1
    }
    
    for (i in (Btheta + 1):(ltheta - 1)) {
      J[i, i] <- -1
      J[i, i + 1] <- 1
    }
    
    J
  }
  
  heq <- function(x) {
    c(
      sum(x[1:ltheta]) - 1,
      sum(x[(ltheta + 1):(ltheta + lsigma)]) - 1
    )
  }
  
  heqjac <- function(x) {
    J <- matrix(0, nrow = 2, ncol = ltheta + lsigma)
    J[1, 1:ltheta] <- 1
    J[2, (ltheta + 1):(ltheta + lsigma)] <- 1
    J
  }
  
  loglik <- log(pmax(likelihood.conditional, 1e-300))
  t1 <- rep(seq_len(ltheta), each = lsigma)
  t2 <- rep(seq_len(lsigma), ltheta)
  
  est.theta <- rep(1 / ltheta, ltheta)
  est.sigma <- rep(1 / lsigma, lsigma)
  
  for (iter in 1:20) {
    est <- as.numeric(array(est.sigma %o% est.theta))
    vv <- t(t(loglik) + log(pmax(est, 1e-300)))
    vv <- exp(vv - apply(vv, 1, max))
    vv <- vv / rowSums(vv)
    
    est <- colMeans(vv)
    est.theta <- as.numeric(tapply(est, t1, sum))
    est.sigma <- as.numeric(tapply(est, t2, sum))
    
    # PAVA unimodal projection at the center
    if (requireNamespace("Iso", quietly = TRUE)) {
      tmp <- Iso::ufit(y = est.theta, x = seq_along(est.theta), lmode = idx0)
      est.theta <- tmp$y
    } else {
      # fallback: triangular feasible init if Iso is unavailable
      left <- est.theta[1:idx0]
      right <- est.theta[idx0:ltheta]
      left <- cummax(left)
      right <- rev(cummax(rev(right)))
      est.theta <- c(left[1:(idx0 - 1)], min(left[idx0], right[1]), right[2:length(right)])
    }
    
    est.theta[est.theta < 0] <- 0
    est.sigma[est.sigma < 0] <- 0
    est.theta <- est.theta / sum(est.theta)
    est.sigma <- est.sigma / sum(est.sigma)
  }
  
  x0 <- c(est.theta, est.sigma)

  fit <- nloptr::slsqp(
    x0 = x0,
    fn = L,
    gr = G,
    lower = rep(0, ltheta + lsigma),
    upper = rep(1, ltheta + lsigma),
    hin = hin,
    hinjac = hinjac,
    heq = heq,
    heqjac = heqjac,
    control = list(
      xtol_rel = 1e-8,
      ftol_rel = 1e-10,
      maxeval = 3000
    ),
    deprecatedBehavior = FALSE,
    nl.info = FALSE
  )
  
  est.theta <- fit$par[1:ltheta]
  est.theta[est.theta < 0] <- 0
  est.theta <- est.theta / sum(est.theta)
  
  est.sigma <- fit$par[(ltheta + 1):(ltheta + lsigma)]
  est.sigma[est.sigma < 0] <- 0
  est.sigma <- est.sigma / sum(est.sigma)
  
  est.matrix <- outer(est.theta, est.sigma)
  est.array <- as.vector(t(est.matrix))
  

  likelihood.theta <- matrix(NA_real_, nrow = p0, ncol = length(grid1))
  likelihood.s2 <- matrix(NA_real_, nrow = p0, ncol = length(grid2))
  
  for (i in seq_len(p0)) {
    likelihood.theta[i, ] <- stats::dnorm(theta0[i], mean = grid1, sd = grid2)
    likelihood.s2[i, ] <- stats::dchisq(df * s20[i] / grid2^2, df = df) * df / grid2^2
  }
  
  likelihood.conditional <- likelihood.theta * likelihood.s2
  LFDR <- matrix(NA_real_, nrow = p0, ncol = ltheta)
  
  for (i in seq_len(p0)) {
    ddd <- likelihood.conditional[i, ] * est.array
    UUU <- numeric(ltheta)
    
    for (j in seq_len(ltheta)) {
      begin <- (j - 1) * lsigma + 1
      end <- j * lsigma
      UUU[j] <- sum(ddd[begin:end])
    }
    
    UUU <- UUU / sum(UUU)
    LFDR[i, ] <- UUU
  }
  
  lfdr <- LFDR[, idx0]
  
  lfsr <- numeric(p0)
  for (i in seq_len(p0)) {
    if (theta0[i] > 0) {
      lfsr[i] <- sum(LFDR[i, 1:idx0])
    } else if (theta0[i] < 0) {
      lfsr[i] <- sum(LFDR[i, idx0:ltheta])
    } else {
      lfsr[i] <- min(sum(LFDR[i, 1:idx0]), sum(LFDR[i, idx0:ltheta]))
    }
  }
  
  list(
    grid.theta = gridtheta,
    grid.sigma2 = gridsigma^2,
    mix.theta = est.theta,
    mix.sigma2 = est.sigma,
    lfdr = lfdr,
    lfsr = lfsr,
    LFDR = LFDR,
    fit = fit
  )
}
