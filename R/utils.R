# =============================================================================
# Shared helpers: data loading, portfolio weights, backtest engine, metrics
# =============================================================================

TICKERS   <- c("PZU", "BNP", "CDR", "OPL", "PKN")      # stooq symbols
LABELS    <- c("PZU", "BNP", "CDR", "ORANGE", "ORLEN")  # names used in output

# ---- Data -------------------------------------------------------------------

# Reads stooq daily CSVs from data/. If a file is missing, tries to download it
# from stooq (free, no API key). Returns an xts of closing prices.
load_prices <- function(data_dir = "data", tickers = TICKERS, labels = LABELS,
                        download_if_missing = TRUE) {
  dir.create(data_dir, showWarnings = FALSE)
  series <- vector("list", length(tickers))
  for (k in seq_along(tickers)) {
    path <- file.path(data_dir, paste0(labels[k], ".csv"))
    if (!file.exists(path) && download_if_missing) {
      url <- sprintf("https://stooq.pl/q/d/l/?s=%s&i=d", tolower(tickers[k]))
      message("Downloading ", labels[k], " from stooq ...")
      utils::download.file(url, path, quiet = TRUE, mode = "wb")
    }
    if (!file.exists(path)) stop("Missing data file: ", path)
    df <- utils::read.csv(path)
    # stooq uses Polish headers (Data, Zamkniecie) on .pl and English on .com
    date_col  <- intersect(c("Data", "Date"), names(df))[1]
    close_col <- intersect(c("Zamkniecie", "Close"), names(df))[1]
    series[[k]] <- xts::xts(df[[close_col]], order.by = as.Date(df[[date_col]]))
    colnames(series[[k]]) <- labels[k]
  }
  stats::na.omit(do.call(merge, series))
}

# Daily simple and log returns from prices
to_returns <- function(prices) {
  simple <- stats::na.omit(prices / stats::lag(prices) - 1)
  log_r  <- stats::na.omit(log(prices / stats::lag(prices)))
  list(simple = simple, log = log_r)
}

# ---- Portfolio weights --------------------------------------------------------

# Long-only minimum-variance weights (sum = 1, w >= 0) via quadratic programming.
# A small ridge keeps the problem well posed if Sigma is near-singular.
mvp_weights_long_only <- function(Sigma, ridge = 1e-10) {
  n <- ncol(Sigma)
  Dmat <- 2 * (Sigma + diag(ridge, n))
  sol <- quadprog::solve.QP(Dmat = Dmat, dvec = rep(0, n),
                            Amat = cbind(rep(1, n), diag(n)),
                            bvec = c(1, rep(0, n)), meq = 1)
  w <- pmax(sol$solution, 0)
  w / sum(w)
}

# ---- Backtest engine ----------------------------------------------------------

# Generic rolling backtest.
#   simple_ret : T x N matrix of simple returns
#   dates      : length-T Date vector
#   rebal_idx  : indices t at which new weights are set (applied from t onwards)
#   weight_fun : function(t) -> list(w = weights, fallback = TRUE/FALSE),
#                using only data strictly before t
#   hold       : holding period in days
# Weights drift with prices between rebalances (buy-and-hold inside the period),
# which is what a real portfolio does. Returns daily portfolio simple returns.
run_backtest <- function(simple_ret, dates, rebal_idx, weight_fun, hold) {
  n_days <- nrow(simple_ret)
  out_ret <- numeric(0); out_date <- as.Date(character(0))
  w_hist <- list(); fallback <- logical(0)
  for (t in rebal_idx) {
    res <- weight_fun(t)
    w <- res$w; fallback <- c(fallback, isTRUE(res$fallback))
    t_end <- min(t + hold - 1, n_days)
    R <- simple_ret[t:t_end, , drop = FALSE]
    # value of each sleeve, starting from target weights
    growth <- apply(1 + R, 2, cumprod)
    if (is.null(dim(growth))) growth <- matrix(growth, nrow = 1)
    value <- as.numeric(growth %*% w)          # portfolio value path (start = 1)
    daily <- value / c(1, head(value, -1)) - 1
    out_ret <- c(out_ret, daily)
    out_date <- c(out_date, dates[t:t_end])
    w_row <- as.data.frame(as.list(setNames(w, colnames(simple_ret))))
    w_hist[[length(w_hist) + 1]] <- data.frame(Date = dates[t], w_row, check.names = FALSE)
  }
  list(returns = xts::xts(out_ret, order.by = out_date),
       weights = do.call(rbind, w_hist),
       fallback_rate = mean(fallback))
}

# ---- Performance metrics ------------------------------------------------------

max_drawdown <- function(r) {
  wealth <- cumprod(1 + as.numeric(r))
  min(wealth / cummax(wealth) - 1)
}

# All metrics computed on daily SIMPLE returns. Risk-free rate = 0, so the
# "Sharpe" column is a return/volatility ratio; this is stated in the README.
perf_stats <- function(r) {
  r <- as.numeric(stats::na.omit(r))
  n <- length(r)
  ann_ret  <- prod(1 + r)^(252 / n) - 1                 # CAGR
  ann_vol  <- stats::sd(r) * sqrt(252)
  down_vol <- sqrt(mean(pmin(r, 0)^2)) * sqrt(252)
  mdd      <- max_drawdown(r)
  var95    <- -stats::quantile(r, 0.05, names = FALSE)  # historical, positive = loss
  es95     <- -mean(r[r <= stats::quantile(r, 0.05)])
  c(CAGR = ann_ret, Volatility = ann_vol,
    Sharpe = ann_ret / ann_vol, Sortino = ann_ret / down_vol,
    MaxDrawdown = mdd, Calmar = ann_ret / abs(mdd),
    VaR95 = var95, ES95 = es95)
}

stats_table <- function(ret_xts) {
  tab <- t(apply(ret_xts, 2, perf_stats))
  data.frame(Strategy = rownames(tab), tab, row.names = NULL, check.names = FALSE)
}

yearly_table <- function(ret_xts) {
  years <- unique(format(zoo::index(ret_xts), "%Y"))
  do.call(rbind, lapply(years, function(y) {
    tab <- stats_table(ret_xts[y])
    cbind(Year = y, tab)
  }))
}

# Percent formatting for markdown tables
fmt_pct <- function(x, d = 2) sprintf(paste0("%.", d, "f%%"), 100 * x)
