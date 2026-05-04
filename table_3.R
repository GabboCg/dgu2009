#!/usr/bin/env Rscript
# ======================================================== #
#
#                  Replication of Table 3 
#
#                 Gabriel E. Cabrera-Guzmán
#                The University of Manchester
#
#                        Spring, 2026
#
#                https://gcabrerag.rbind.io
#
# ------------------------------ #
# email: gabriel.cabreraguzman@postgrad.manchester.ac.uk
# ======================================================== #

# Replication of Table 3 from DeMiguel, Garlappi, Uppal (RFS, 2009)
# "Optimal Versus Naive Diversification: How Inefficient is the 1/N Portfolio Strategy?"
#
# Table 3: Sharpe ratios for empirical data
# Settings: M = 120 months estimation window, gamma = 1, rolling window

# Load packages
library(quadprog)

# ==========================================
#               Configuration
# ------------------------------------------

DATA_DIR <- "data"
M_WINDOW <- 120 # estimation window length
GAMMA <- 1      # risk aversion coefficient

# Dataset names (correspond to .rds files in data/)
DATASET_NAMES <- c("SPSectors", "Industry", "International", "MKT_SMB_HML", "FF_1factor", "FF_4factor")

# ==========================================
#               Helper Functions
# ------------------------------------------

load_data <- function(ds_name) {
    
  filepath <- file.path(DATA_DIR, paste0(ds_name, ".rds"))
  dat <- readRDS(filepath)
  raw <- dat$data
  nFactors <- dat$nFactors
  rRiskfree <- raw[, 2]
  rRisky <- as.matrix(raw[, 3:ncol(raw)])
  nCols <- ncol(rRisky)
  
  if (nFactors > 0) {
      
    rFactor <- as.matrix(rRisky[, (nCols - nFactors + 1):nCols])
    
  } else {
    
    rFactor <- NULL
    
  }
  
  n <- 1 + nCols
  nRisky <- n - 1
  
  list(rRiskfree = rRiskfree, rRisky = rRisky, rFactor = rFactor, n = n, nRisky = nRisky, nFactors = nFactors)
  
}

# Out-of-sample return
out_sample_return <- function(alpha, rRisky, M, j) {
    
  drop(crossprod(alpha, rRisky[M + j,]))
    
}

# Sharpe ratio from return series
sharpe_ratio <- function(xrp) {
    
  mn <- mean(xrp)
  sd_val <- sd(xrp)
  sr <- mn / sd_val
  
  list(mean = mn, sd = sd_val, sr = sr)
  
}

# Jobson-Korkie test statistic
jk_stat <- function(r1, r2) {
    
  mu1 <- mean(r1)
  mu2 <- mean(r2)
  S   <- cov(cbind(r1, r2))
  s1  <- sqrt(S[1, 1])
  s2  <- sqrt(S[2, 2])
  s12 <- S[1, 2]
  TT  <- length(r1)
  theta <- (1 / TT) * (2 * s1 ^ 2 * s2 ^ 2 - 2 * s1 * s2 * s12 + 0.5 * mu1 ^ 2 * s2 ^ 2 + 0.5 * mu2 ^ 2 * s1 ^ 2 - (mu1 * mu2 / (s1 * s2)) * s12 ^ 2)
  
  (s2 * mu1 - s1 * mu2) / max(.Machine$double.eps, sqrt(theta))
  
}

# p-value for JK stat (one-sided, matching Memmel 2003 / MATLAB code)
jk_pvalue <- function(r1, r2) {
    
  z <- jk_stat(r1, r2)
  1 - pnorm(abs(z))
  
}

# Minimum variance constrained (long-only, weights sum to 1)
min_var_constrained <- function(Sigma) {
    
  nn <- nrow(Sigma)
  Dmat <- Sigma
  dvec <- rep(0, nn)
  
  # Constraints: sum(w) = 1, w >= 0
  Amat <- cbind(rep(1, nn), diag(nn))
  bvec <- c(1, rep(0, nn))
  meq <- 1
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
  sol$solution
  
}

# Jagannathan-Ma: constrained min variance with lb = 1/(2n)
jagannathan_ma <- function(Sigma) {
    
  nn <- nrow(Sigma)
  lb <- 1 / (2 * nn)
  Dmat <- Sigma
  dvec <- rep(0, nn)
  
  # sum(w) = 1, w >= lb, w <= 1
  Amat <- cbind(rep(1, nn), diag(nn), -diag(nn))
  bvec <- c(1, rep(lb, nn), rep(-1, nn))
  meq <- 1
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
  sol$solution
  
}

# Mean-variance constrained (long-only, weights sum <= 1, w_i <= 1)
mv_constrained <- function(mu_risky, Sigma, gamma) {
    
  nn <- length(mu_risky)
  Dmat <- gamma * Sigma
  dvec <- mu_risky
  # Constraints: sum(w) <= 1, w >= 0, w <= 1
  # solve.QP requires A'x >= b
  Amat <- cbind(-rep(1, nn), diag(nn), -diag(nn))
  bvec <- c(-1, rep(0, nn), rep(-1, nn))
  meq <- 0
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
  sol$solution
  
}

# Kan-Zhou portfolio (eq. 67 in Kan-Zhou 2007)
kan_zhou <- function(N, TT, muhat, invSigma, gamma) {
    
  ones_N <- rep(1, N)
  mug <- drop(crossprod(muhat, invSigma %*% ones_N)) / drop(crossprod(ones_N, invSigma %*% ones_N))
  diff_mu <- muhat - mug
  phiHat2 <- drop(crossprod(diff_mu, invSigma %*% diff_mu))
  x <- phiHat2 / (1 + phiHat2)
  B <- pbeta(x, (N - 1) / 2, (TT - N + 1) / 2)
  phiHat2a <- ((TT - N - 1) * phiHat2 - (N - 1)) / TT + (2 * phiHat2 ^ ((N - 1) / 2) * (1 + phiHat2) ^ (-(TT - 2) / 2)) / (TT * B)
  c3 <- ((TT - N - 1) * (TT - N - 4)) / (TT * (TT - 2))
  denom <- phiHat2a + N / TT
  alpha <- (c3 / gamma) * ((phiHat2a / denom) * invSigma %*% muhat + ((N / TT) / denom) * mug * invSigma %*% ones_N)
  
  drop(alpha)
  
}

# Kan-Zhou 1/N (ew-min combination)
kan_zhou_1overN <- function(N, TT, Sigma) {
    
  invSigma <- solve(Sigma)
  ones_N <- rep(1, N)
  esige <- drop(crossprod(ones_N, Sigma %*% ones_N))
  einvsige <- drop(crossprod(ones_N, invSigma %*% ones_N))
  k <- (TT^2 * (TT - 2)) / ((TT - N - 1) * (TT - N - 2) * (TT - N - 4))
  d_val <- ((TT - N - 2) * esige * einvsige - N ^ 2 * TT) / (N ^ 2 * (TT - N - 2) * k * einvsige - 2 * TT * N ^ 2 * einvsige + (TT - N - 2) * einvsige ^ 2 * esige)
  
  c_val <- 1 - d_val * einvsige
  alpha <- c_val * (1 / N) * ones_N + d_val * invSigma %*% ones_N
  
  drop(alpha)
  
}

# Pastor-MacKinlay portfolio (with Kan-Zhou approximation)
pastor_mackinlay <- function(nRisky, nPoints, muHat, SigmaHat, gamma) {
    
  Uhat <- SigmaHat + muHat %*% t(muHat)
  eig <- eigen(Uhat, symmetric = TRUE)
  l1 <- eig$values[1]
  q1 <- eig$vectors[, 1]
  muTilde <- drop(crossprod(q1, muHat)) * q1
  
  (1 / gamma) * muTilde / (l1 - drop(crossprod(muTilde, muTilde)))
  
}

# Pastor (Data-and-Model) approach
pastor_dm <- function(rRisky_noFac, rFactor, sigma_alpha_annual) {
    
  mu1hat <- colMeans(rRisky_noFac)
  mu2hat <- colMeans(rFactor)
  omega22hat <- cov(rFactor)
  TT <- nrow(rRisky_noFac)
  m <- ncol(rRisky_noFac)
  k <- ncol(rFactor)

  Shat2 <- drop(t(mu2hat) %*% solve(omega22hat) %*% mu2hat)
  delta_bar <- (TT * (TT - 2) + k) / (TT * (TT - k - 2)) - (k + 3) * Shat2 / (TT * (TT - k - 2) * (1 + Shat2))
  delta_hat <- (TT - 2) * (TT + 1) / (TT * (TT - k - 2))
  b_val <- (TT + 1) / (TT - k - 2)
  h_val <- TT / (TT - m - k - 1)

  # --- Regression WITH intercept ---
  X_int <- cbind(1, rFactor)
  beta_h <- solve(crossprod(X_int), crossprod(X_int, rRisky_noFac))
  R <- rRisky_noFac - X_int %*% beta_h
  beta_hat <- t(beta_h[2:(k + 1), , drop = FALSE])
  Sigma_hat <- cov(R)

  # --- Regression WITHOUT intercept ---
  beta_bar_raw <- solve(crossprod(rFactor), crossprod(rFactor, rRisky_noFac))
  R1 <- rRisky_noFac - rFactor %*% beta_bar_raw
  beta_bar <- t(beta_bar_raw)
  Sigma_bar <- cov(R1)

  s2 <- max(diag(Sigma_hat))

  s_alpha <- sigma_alpha_annual / 12
  omega_vec <- 1 / (1 + s_alpha ^ 2 * TT / ((1 + Shat2) * s2))

  mu_list    <- list()
  Sigma_list <- list()
  
  for (i in seq_along(sigma_alpha_annual)) {
      
    om <- omega_vec[i]
    beta_mix <- om * beta_bar + (1 - om) * beta_hat
    v11 <- b_val * beta_mix %*% omega22hat %*% t(beta_mix) + h_val * (om * delta_bar + (1 - om) * delta_hat) * (om * Sigma_bar + (1 - om) * Sigma_hat)
    v12 <- b_val * beta_mix %*% omega22hat
    mu_combined <- om * beta_bar %*% mu2hat + (1 - om) * mu1hat
    mu_fac <- mu2hat
    mu_all <- c(drop(mu_combined), drop(mu_fac))
    
    V_full <- rbind(cbind(v11, v12), cbind(t(v12), b_val * omega22hat))
    
    mu_list[[i]] <- mu_all
    Sigma_list[[i]] <- V_full
    
  }
  
  list(mu = mu_list, Sigma = Sigma_list)
  
}

# In-sample Sharpe ratio (using out-of-sample period data, risky-only)
# Matches MATLAB: uses data from upperM+1:end, computes MV weights on that,
# then normalizes to risky-only
in_sample_sr <- function(rRisky, M, nRisky, n, gamma) {
    
  # Use data from the out-of-sample period (after estimation window)
  risky <- rRisky[(M + 1):nrow(rRisky),]
  nPoints <- nrow(risky)
  mu <- colMeans(risky)
  
  # --- Unbiased estimator of Sigma ---
  Sigma <- cov(risky) * (nPoints - 1) / (nPoints - nRisky - 2)
  invSigma <- solve(Sigma)
  
  # --- MV weights ---
  w_mv <- drop((1 / gamma) * invSigma %*% mu)
  
  # Risky-only: divide by sum (not abs(sum))
  w_ro <- w_mv / sum(w_mv)
  mn_is <- drop(crossprod(w_ro, mu))
  sd_is <- sqrt(drop(t(w_ro) %*% Sigma %*% w_ro))
  abs(mn_is / sd_is) # abs per MATLAB code
  
}

# ==========================================
#            Main Analysis Function
# ------------------------------------------

run_analysis <- function(ds_name, M = M_WINDOW, gamma = GAMMA) {
    
  dat <- load_data(ds_name)
  rRiskfree <- dat$rRiskfree
  rRisky    <- dat$rRisky
  rFactor   <- dat$rFactor
  n         <- dat$n
  nRisky    <- dat$nRisky
  nFactors  <- dat$nFactors
  TT        <- nrow(rRisky)
  nSubsets  <- TT - M
  nPoints   <- M # estimation window size

  # Strategy names
  strat_names <- c("ew", "mv", "bs", "min", "vw", "mp", "mv_c", "bs_c", "min_c", "g_min_c", "mv_min", "ew_min")
  if (nFactors > 0) strat_names <- c(strat_names, "dm")

  # Storage for out-of-sample excess returns (risky-only)
  xrp <- setNames(lapply(strat_names, function(s) numeric(nSubsets)), strat_names)

  sigma_alpha_vals <- c(999, 0.02, 0.01, 0.005, 0)

  cat(sprintf("  Processing %d rolling windows...\n", nSubsets))

  for (j in seq_len(nSubsets)) {
      
    # Extract estimation window
    idx <- j:(M + j - 1)
    rRisky_sub    <- rRisky[idx,]
    rRiskfree_sub <- rRiskfree[idx]

    # ---- Moments ----
    mu <- c(mean(rRiskfree_sub), colMeans(rRisky_sub))
    mu_risky <- mu[2:n]

    # Unbiased estimator of Sigma (for MV, BS)
    Sigma <- cov(rRisky_sub) * (nPoints - 1) / (nPoints - nRisky - 2)
    invSigma <- solve(Sigma)

    # MLE estimator of Sigma (for MinVar, KZ, PM)
    SigmaMLE <- cov(rRisky_sub) * (nPoints - 1) / nPoints
    invSigmaMLE <- solve(SigmaMLE)
    AMLE <- drop(crossprod(rep(1, nRisky), invSigmaMLE %*% rep(1, nRisky)))

    # ---- Bayes-Stein moments (Jorion 1986) ----
    Y <- mu_risky
    Ahat <- drop(crossprod(rep(1, nRisky), invSigma %*% rep(1, nRisky)))
    Y0 <- drop(crossprod(rep(1, nRisky), invSigma %*% Y)) / Ahat
    w_bs <- drop((nRisky + 2) / ((nRisky + 2) + t(Y - Y0) %*% (nPoints * invSigma) %*% (Y - Y0)))
    lambda_bs <- drop((nRisky + 2) / (t(Y - Y0) %*% invSigma %*% (Y - Y0)))
    muBS <- c(mean(rRiskfree_sub), (1 - w_bs) * Y + w_bs * Y0)
    SigmaBS <- Sigma * (1 + 1 / (nPoints + lambda_bs)) + (lambda_bs / (nPoints * (nPoints + 1 + lambda_bs))) * matrix(1, nRisky, nRisky) / Ahat
    invSigmaBS <- solve(SigmaBS)

    # Helper: convert weights to risky-only (normalize by |sum(w)|)
    to_ro <- function(w) w / abs(sum(w))

    # ---- Portfolio Weights (original, then convert to risky-only) ----
    # 1/N: risky-only = 1/nRisky each
    alpha_ew_ro <- rep(1 / nRisky, nRisky)

    # Sample mean-variance
    alpha_mv <- drop((1 / gamma) * invSigma %*% mu_risky)
    alpha_mv_ro <- to_ro(alpha_mv)

    # Bayes-Stein
    alpha_bs <- drop((1 / gamma) * invSigmaBS %*% muBS[2:n])
    alpha_bs_ro <- to_ro(alpha_bs)

    # Minimum variance (unconstrained, uses MLE Sigma)
    alpha_min <- drop((1 / AMLE) * invSigmaMLE %*% rep(1, nRisky))
    # min weights already sum to 1, so risky-only = same
    alpha_min_ro <- alpha_min

    # Value-weighted: market portfolio (last non-factor column)
    # In risky-only, vw = single-asset return of the market index
    alpha_vw_ro <- rep(0, nRisky)
    
    if (nFactors > 0) {
        
      # Market factor is last column in rRisky
      alpha_vw_ro[nRisky] <- 1
      
    } else {
        
      alpha_vw_ro[nRisky] <- 1
      
    }

    # MacKinlay-Pastor (mp)
    alpha_mp <- pastor_mackinlay(nRisky, nPoints, mu_risky, SigmaMLE, gamma)
    alpha_mp_ro <- to_ro(alpha_mp)

    # MV-Constrained (long-only)
    alpha_mv_c <- mv_constrained(mu_risky, Sigma, gamma)
    alpha_mv_c_ro <- to_ro(alpha_mv_c)

    # BS-Constrained (long-only)
    alpha_bs_c <- mv_constrained(muBS[2:n], SigmaBS, gamma)
    alpha_bs_c_ro <- to_ro(alpha_bs_c)

    # Min-Variance Constrained (long-only, weights already sum to 1)
    alpha_min_c <- min_var_constrained(SigmaMLE)
    alpha_min_c_ro <- alpha_min_c

    # Jagannathan-Ma (g-min-c, weights already sum to 1)
    alpha_gminc <- jagannathan_ma(Sigma)
    alpha_gminc_ro <- alpha_gminc

    # Kan-Zhou mv-min (three-fund)
    alpha_kz <- kan_zhou(nRisky, nPoints, mu_risky, invSigmaMLE, gamma)
    alpha_kz_ro <- to_ro(alpha_kz)

    # Kan-Zhou 1/N (ew-min)
    alpha_ewmin <- kan_zhou_1overN(nRisky, M, Sigma)
    alpha_ewmin_ro <- to_ro(alpha_ewmin)

    # Data-and-Model (dm) with sigma_alpha = 1% per year (index 3)
    alpha_dm_ro <- NULL
    if (nFactors > 0) {
        
      rFactor_sub <- rFactor[idx, , drop = FALSE]
      nRisky_noFac <- nRisky - nFactors
      rRisky_noFac_sub <- rRisky_sub[, 1:nRisky_noFac]
      dm_result <- pastor_dm(rRisky_noFac_sub, rFactor_sub, sigma_alpha_vals)
      mu_dm <- dm_result$mu[[3]]
      Sigma_dm <- dm_result$Sigma[[3]]
      alpha_dm <- drop((1 / gamma) * solve(Sigma_dm) %*% mu_dm)
      
      # DM weights cover nRisky_noFac + nFactors assets
      alpha_dm_full <- c(alpha_dm, rep(0, nRisky - length(alpha_dm)))
      alpha_dm_ro <- to_ro(alpha_dm_full)
      
    }

    # ---- Out-of-sample returns (risky-only) ----
    ret_next <- rRisky[M + j,]
    xrp$ew[j]      <- drop(crossprod(alpha_ew_ro, ret_next))
    xrp$mv[j]      <- drop(crossprod(alpha_mv_ro, ret_next))
    xrp$bs[j]      <- drop(crossprod(alpha_bs_ro, ret_next))
    xrp$min[j]     <- drop(crossprod(alpha_min_ro, ret_next))
    xrp$vw[j]      <- drop(crossprod(alpha_vw_ro, ret_next))
    xrp$mp[j]      <- drop(crossprod(alpha_mp_ro, ret_next))
    xrp$mv_c[j]    <- drop(crossprod(alpha_mv_c_ro, ret_next))
    xrp$bs_c[j]    <- drop(crossprod(alpha_bs_c_ro, ret_next))
    xrp$min_c[j]   <- drop(crossprod(alpha_min_c_ro, ret_next))
    xrp$g_min_c[j] <- drop(crossprod(alpha_gminc_ro, ret_next))
    xrp$mv_min[j]  <- drop(crossprod(alpha_kz_ro, ret_next))
    xrp$ew_min[j]  <- drop(crossprod(alpha_ewmin_ro, ret_next))
    
    if (!is.null(alpha_dm_ro)) {
        
      xrp$dm[j] <- drop(crossprod(alpha_dm_ro, ret_next))
      
    }
    
  }

  # ---- Compute Sharpe ratios and p-values ----
  results <- data.frame(strategy = character(), sr = numeric(),
                        pvalue = numeric(), stringsAsFactors = FALSE)

  for (s in strat_names) {
      
    sr_info <- sharpe_ratio(xrp[[s]])
    
    if (s == "ew") {
        
      pv <- NA  # benchmark, no p-value
      
    } else {
        
      pv <- jk_pvalue(xrp$ew, xrp[[s]])
      
    }
    
    results <- rbind(results, data.frame(strategy = s, sr = sr_info$sr, pvalue = pv, stringsAsFactors = FALSE))
    
  }

  # In-sample MV Sharpe ratio
  sr_is <- in_sample_sr(rRisky, M, nRisky, n, gamma)

  list(results = results, sr_insample = sr_is, n = n, nRisky = nRisky)
  
}

# ==========================================
#               Run All Datasets
# ------------------------------------------

cat("=================================================================\n")
cat("Replicating Table 3: Sharpe ratios for empirical data\n")
cat("DeMiguel, Garlappi, Uppal (RFS, 2009)\n")
cat("Estimation window M =", M_WINDOW, ", gamma =", GAMMA, "\n")
cat("=================================================================\n\n")

all_results <- list()

for (ds_name in DATASET_NAMES) {
    
  cat(sprintf("Dataset: %s\n", ds_name))
  res <- tryCatch(
    run_analysis(ds_name),
    error = function(e) {
        
      cat(sprintf("  ERROR: %s\n", e$message))
      NULL
      
    }
  )
  if (!is.null(res)) {
      
    all_results[[ds_name]] <- res
    cat(sprintf("  N = %d (risky), done.\n\n", res$nRisky))
    
  }
  
}

# ==========================================
#         Format and Display Table 3
# ------------------------------------------

# Map strategy codes to paper labels
label_map <- c(
  ew      = "1/N",
  mv_is   = "mv (in sample)",
  mv      = "mv",
  bs      = "bs",
  dm      = "dm (sa = 1.0%)",
  min     = "min",
  vw      = "vw",
  mp      = "mp",
  mv_c    = "mv-c",
  bs_c    = "bs-c",
  min_c   = "min-c",
  g_min_c = "g-min-c",
  mv_min  = "mv-min",
  ew_min  = "ew-min"
)

# Table 3 row order (matching paper)
row_order <- c("ew", "mv_is", "mv", "bs", "dm", "min", "vw", "mp",
               "mv_c", "bs_c", "min_c", "g_min_c", "mv_min", "ew_min")

# Column headers
ds_labels <- c(
  SPSectors   = "S&P sectors",
  Industry    = "Industry",
  International = "Inter'l",
  MKT_SMB_HML = "Mkt/SMB/HML",
  FF_1factor  = "FF 1-factor",
  FF_4factor  = "FF 4-factor"
)

cat("\n")
cat("Table 3\n")
cat("Sharpe ratios for empirical data\n")
cat(sprintf("%-20s", "Strategy"))
for (ds in names(all_results)) {
    
  cat(sprintf("%14s", ds_labels[ds]))
    
}

cat("\n")
cat(sprintf("%-20s", ""))

for (ds in names(all_results)) {
    
  cat(sprintf("%14s", paste0("N = ", all_results[[ds]]$nRisky)))
    
}

cat("\n")
cat(paste(rep("-", 20 + 14 * length(all_results)), collapse = ""), "\n")

for (strat in row_order) {
    
  lab <- label_map[strat]
  cat(sprintf("%-20s", lab))

  for (ds in names(all_results)) {
      
    r <- all_results[[ds]]
    
    if (strat == "mv_is") {
        
      cat(sprintf("%14.4f", r$sr_insample))
        
    } else {
        
      idx <- which(r$results$strategy == strat)
      
      if (length(idx) > 0) {
          
        cat(sprintf("%14.4f", r$results$sr[idx]))
          
      } else {
          
        cat(sprintf("%14s", "---"))
          
      }
      
    }
    
  }
  
  cat("\n")

  # p-value row (skip for 1/N and mv in-sample)
  if (!(strat %in% c("ew", "mv_is"))) {
      
    cat(sprintf("%-20s", ""))
    
    for (ds in names(all_results)) {
        
      r <- all_results[[ds]]
      idx <- which(r$results$strategy == strat)
      
      if (length(idx) > 0 && !is.na(r$results$pvalue[idx])) {
          
        cat(sprintf("%14s", sprintf("(%.2f)", r$results$pvalue[idx])))
          
      } else {
          
        cat(sprintf("%14s", ""))
          
      }
      
    }
      
    cat("\n")
    
  }
  
}

cat(paste(rep("-", 20 + 14 * length(all_results)), collapse = ""), "\n")
cat("\nNote: p-values (in parentheses) are for the JK test of the difference\n")
cat("in Sharpe ratio from the 1/N benchmark.\n")
