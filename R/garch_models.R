# =============================================================================
# Covariance forecasts: DCC-GARCH (normal / Student-t) and full BEKK(1,1)
# Every function returns list(Sigma = N x N matrix, fallback = TRUE/FALSE).
# fallback = TRUE means the model failed and the sample covariance was used.
# Fallbacks are counted and reported, so a "GARCH" result can never silently
# be the static model in disguise.
# =============================================================================

# VAR(1) mean filter: returns residuals and the coefficients needed to compute
# residuals for new observations later on.
fit_var1 <- function(X) {
  Y  <- X[-1, , drop = FALSE]
  Z  <- cbind(1, X[-nrow(X), , drop = FALSE])
  B  <- solve(crossprod(Z), crossprod(Z, Y))       # (1+N) x N
  list(coef = B, resid = Y - Z %*% B)
}
var1_resid <- function(B, x_now, x_prev) as.numeric(x_now - c(1, x_prev) %*% B)

# ---- DCC ---------------------------------------------------------------------

make_dcc_spec <- function(n_assets, dist = c("mvnorm", "mvt")) {
  dist <- match.arg(dist)
  u <- rugarch::ugarchspec(
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    distribution.model = if (dist == "mvt") "std" else "norm")
  rmgarch::dccspec(uspec = rugarch::multispec(replicate(n_assets, u)),
                   dccOrder = c(1, 1), distribution = dist)
}

dcc_forecast <- function(X, spec) {
  fallback <- list(Sigma = stats::cov(X), fallback = TRUE)
  tryCatch({
    res <- fit_var1(X)$resid
    fit <- rmgarch::dccfit(spec, data = res, solver = c("solnp", "nlminb"))
    if (fit@mfit$convergence != 0) return(fallback)
    H <- rmgarch::rcov(rmgarch::dccforecast(fit, n.ahead = 1))[[1]][, , 1]
    if (any(!is.finite(H))) return(fallback)
    list(Sigma = H, fallback = FALSE)
  }, error = function(e) fallback)
}

# ---- BEKK --------------------------------------------------------------------
# mgarchBEKK estimates the model but has no forecast method, so the recursion is
# implemented here. Because the package's matrix convention (C vs C', A vs A')
# is easy to get wrong, the code tries the four possible conventions and keeps
# the one that reproduces the package's own fitted H_T. If none does, the
# window is flagged as a fallback instead of silently using wrong forecasts.

bekk_step <- function(H, e, C, A, B, conv) {
  CC <- if (conv %in% c(1, 3)) C %*% t(C) else t(C) %*% C
  ee <- e %*% t(e)
  if (conv %in% c(1, 2)) CC + t(A) %*% ee %*% A + t(B) %*% H %*% B
  else                   CC + A %*% ee %*% t(A) + B %*% H %*% t(B)
}

bekk_fit <- function(resid) {
  fit <- mgarchBEKK::BEKK(resid, order = c(1, 1), method = "BFGS", verbose = FALSE)
  P <- fit$est.params
  C <- P[[1]]; A <- P[[2]]; B <- P[[3]]
  Hs <- fit$H.estimated
  get_H <- function(k) if (is.list(Hs)) Hs[[k]] else Hs[, , k]
  T_ <- if (is.list(Hs)) length(Hs) else dim(Hs)[3]
  H_prev <- get_H(T_ - 1); H_last <- get_H(T_)
  e_prev <- resid[nrow(resid) - 1, ]
  errs <- sapply(1:4, function(cv)
    norm(bekk_step(H_prev, e_prev, C, A, B, cv) - H_last, "F") / norm(H_last, "F"))
  conv <- which.min(errs)
  if (errs[conv] > 1e-3) stop("BEKK convention check failed (rel. error ",
                              signif(errs[conv], 3), ")")
  list(C = C, A = A, B = B, conv = conv, H_last = H_last,
       e_last = resid[nrow(resid), ])
}

# Estimates BEKK once on the window ending at t0 - 1, then produces one-step
# covariance forecasts for each rebalance date in `rebal_block`, updating H with
# the realised VAR residuals between re-estimations (parameters stay fixed).
bekk_block <- function(X_all, t0, rebal_block, window_len) {
  idx <- (t0 - window_len):(t0 - 1)
  X <- X_all[idx, , drop = FALSE]
  sample_cov <- stats::cov(X)
  out <- lapply(rebal_block, function(t) list(t = t, Sigma = sample_cov, fallback = TRUE))
  est <- tryCatch({
    v <- fit_var1(X)
    list(var = v, bekk = bekk_fit(v$resid))
  }, error = function(e) NULL)
  if (is.null(est)) return(out)

  H <- est$bekk$H_last; e <- est$bekk$e_last; s <- t0   # s = next day to forecast
  for (k in seq_along(rebal_block)) {
    t <- rebal_block[k]
    # roll the recursion forward through realised days s .. t-1
    while (s < t) {
      H <- bekk_step(H, e, est$bekk$C, est$bekk$A, est$bekk$B, est$bekk$conv)
      e <- var1_resid(est$var$coef, X_all[s, ], X_all[s - 1, ])
      s <- s + 1
    }
    H_fc <- bekk_step(H, e, est$bekk$C, est$bekk$A, est$bekk$B, est$bekk$conv)
    if (all(is.finite(H_fc))) out[[k]] <- list(t = t, Sigma = H_fc, fallback = FALSE)
  }
  out
}
