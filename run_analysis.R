# =============================================================================
# Minimum-variance portfolio on 5 WSE stocks: static vs DCC-GARCH vs BEKK
# Run from the repository root:  Rscript run_analysis.R
# Outputs: results/*.csv, results/summary.md, figures/*.png
# =============================================================================

pkgs <- c("xts", "zoo", "quadprog", "rugarch", "rmgarch", "mgarchBEKK",
          "ggplot2", "tidyr", "dplyr", "scales", "parallel")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

source("R/utils.R")
source("R/garch_models.R")

# ---- Settings -----------------------------------------------------------------
WINDOW     <- 500                    # estimation window (trading days)
HOLD       <- 5                      # rebalance every 5 days
BEKK_EVERY <- 50                     # BEKK re-estimated every 50 days
START      <- as.Date("2019-01-01")  # first out-of-sample day
END        <- as.Date("2025-12-31")
N_CORES    <- max(1, parallel::detectCores() - 1)
RUN_GARCH  <- TRUE                   # FALSE = only static MVP vs EW (seconds)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- Data ---------------------------------------------------------------------
prices <- load_prices("data")
rets   <- to_returns(prices)
R_s    <- as.matrix(rets$simple)              # simple returns for P&L
R_l    <- as.matrix(rets$log) * 100           # log returns in % for estimation
dates  <- zoo::index(rets$simple)
N      <- ncol(R_s)

start_idx <- which(dates >= START)[1]
end_idx   <- max(which(dates <= END))
stopifnot(start_idx > WINDOW)
rebal_idx <- seq(start_idx, end_idx, by = HOLD)
cat(sprintf("Data: %s to %s | out-of-sample from %s | %d rebalances\n",
            min(dates), max(dates), dates[start_idx], length(rebal_idx)))

# ---- Covariance forecasts per rebalance date (cached) --------------------------
window_of <- function(t) R_l[(t - WINDOW):(t - 1), , drop = FALSE]

sigma_static <- lapply(rebal_idx, function(t)
  list(t = t, Sigma = cov(window_of(t)), fallback = FALSE))

cached <- function(file, expr) {
  if (file.exists(file)) return(readRDS(file))
  res <- expr; saveRDS(res, file); res
}

if (RUN_GARCH) {
  cl <- parallel::makeCluster(N_CORES)
  parallel::clusterEvalQ(cl, { library(rugarch); library(rmgarch); library(mgarchBEKK) })
  parallel::clusterExport(cl, c("R_l", "WINDOW", "window_of", "fit_var1",
                                "var1_resid", "dcc_forecast", "make_dcc_spec",
                                "bekk_step", "bekk_fit", "bekk_block", "N"))

  for (d in c("mvnorm", "mvt")) {
    f <- sprintf("results/cache_dcc_%s.rds", d)
    if (!file.exists(f)) cat("Estimating DCC", d, "(this takes a while) ...\n")
    parallel::clusterExport(cl, "d", envir = environment())
    assign(paste0("sigma_dcc_", d), cached(f, parallel::parLapply(cl, rebal_idx, function(t) {
      fc <- dcc_forecast(window_of(t), make_dcc_spec(N, d))
      c(list(t = t), fc)
    })))
  }

  est_points <- seq(start_idx, end_idx, by = BEKK_EVERY)
  blocks <- lapply(est_points, function(t0) rebal_idx[rebal_idx >= t0 & rebal_idx < t0 + BEKK_EVERY])
  parallel::clusterExport(cl, c("est_points", "blocks"))
  if (!file.exists("results/cache_bekk.rds")) cat("Estimating BEKK (slowest step) ...\n")
  sigma_bekk <- cached("results/cache_bekk.rds", unlist(parallel::parLapply(
    cl, seq_along(est_points), function(k) bekk_block(R_l, est_points[k], blocks[[k]], WINDOW)
  ), recursive = FALSE))
  parallel::stopCluster(cl)
}

# ---- Backtests ----------------------------------------------------------------
lookup_weights <- function(sig_list) {
  by_t <- setNames(sig_list, sapply(sig_list, `[[`, "t"))
  function(t) {
    s <- by_t[[as.character(t)]]
    list(w = mvp_weights_long_only(s$Sigma), fallback = s$fallback)
  }
}

models <- list(MVP_Static = sigma_static)
if (RUN_GARCH) models <- c(models, list(DCC_Normal = sigma_dcc_mvnorm,
                                        DCC_Student = sigma_dcc_mvt,
                                        BEKK = sigma_bekk))

bt <- lapply(models, function(s) run_backtest(R_s, dates, rebal_idx, lookup_weights(s), HOLD))
bt$EW <- run_backtest(R_s, dates, rebal_idx, function(t) list(w = rep(1 / N, N), fallback = FALSE), HOLD)

strat_ret <- do.call(merge, lapply(names(bt), function(nm) setNames(bt[[nm]]$returns, nm)))
single    <- rets$simple[zoo::index(strat_ret)]
colnames(single) <- paste0("Single_", colnames(single))
all_ret   <- merge(strat_ret, single)

fallbacks <- data.frame(Model = names(bt), FallbackRate = sapply(bt, `[[`, "fallback_rate"))
print(fallbacks)

# ---- Tables -------------------------------------------------------------------
overall <- stats_table(all_ret)
yearly  <- yearly_table(all_ret)
weights_static <- bt$MVP_Static$weights
w_summary <- data.frame(
  Asset = LABELS,
  Mean = colMeans(weights_static[, -1]),
  Min  = apply(weights_static[, -1], 2, min),
  Max  = apply(weights_static[, -1], 2, max), row.names = NULL)

write.csv(overall,   "results/stats_overall.csv", row.names = FALSE)
write.csv(yearly,    "results/stats_yearly.csv",  row.names = FALSE)
write.csv(w_summary, "results/mvp_weights_summary.csv", row.names = FALSE)
write.csv(fallbacks, "results/model_fallback_rates.csv", row.names = FALSE)
write.csv(data.frame(Date = zoo::index(all_ret), coredata(all_ret)),
          "results/daily_returns.csv", row.names = FALSE)
print(overall, digits = 3)

# Markdown summary to paste into README
md_tab <- function(df) {
  num <- sapply(df, is.numeric)
  df[num] <- lapply(names(df)[num], function(nm)
    if (nm %in% c("Sharpe", "Sortino", "Calmar")) sprintf("%.2f", df[[nm]]) else fmt_pct(df[[nm]]))
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
writeLines(c("## Full period", md_tab(overall), "",
             "## Average MVP weights (static model)", md_tab(w_summary), "",
             "## Model fallback rates", md_tab(fallbacks)), "results/summary.md")

# ---- Figures ------------------------------------------------------------------
theme_set(theme_minimal(base_size = 12) +
          theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
nice <- c(MVP_Static = "MVP (sample covariance)", DCC_Normal = "MVP (DCC, normal)",
          DCC_Student = "MVP (DCC, Student-t)", BEKK = "MVP (BEKK)", EW = "Equal weight",
          Single_PZU = "PZU", Single_BNP = "BNP Paribas", Single_CDR = "CD Projekt",
          Single_ORANGE = "Orange Polska", Single_ORLEN = "Orlen")
key <- c("MVP_Static", "DCC_Normal", "DCC_Student", "BEKK", "EW")
pal <- setNames(c("#1b3a6b", "#2a9d8f", "#8ab17d", "#9b5de5", "#e76f51"), nice[key])
relabel <- function(x) factor(nice[x], levels = nice[key])

wealth <- data.frame(Date = zoo::index(all_ret), apply(1 + coredata(all_ret[, key]), 2, cumprod)) |>
  pivot_longer(-Date, names_to = "Strategy", values_to = "Value") |> mutate(Strategy = relabel(Strategy))
p1 <- ggplot(wealth, aes(Date, Value, colour = Strategy)) + geom_line(linewidth = 0.8) +
  scale_colour_manual(values = pal) + scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(title = "Growth of 1 PLN, 2019-2025", x = NULL, y = "Portfolio value (PLN)", colour = NULL) +
  guides(colour = guide_legend(nrow = 2))
ggsave("figures/equity_curves.png", p1, width = 9, height = 5, dpi = 150)

dd <- data.frame(Date = zoo::index(all_ret), apply(1 + coredata(all_ret[, key]), 2,
                 function(x) { w <- cumprod(x); w / cummax(w) - 1 })) |>
  pivot_longer(-Date, names_to = "Strategy", values_to = "DD") |> mutate(Strategy = relabel(Strategy))
p2 <- ggplot(dd, aes(Date, DD, colour = Strategy)) + geom_line(linewidth = 0.6) +
  scale_colour_manual(values = pal) + scale_y_continuous(labels = percent) +
  labs(title = "Drawdown: fall from the previous peak", x = NULL, y = NULL, colour = NULL) +
  guides(colour = guide_legend(nrow = 2))
ggsave("figures/drawdowns.png", p2, width = 9, height = 4.5, dpi = 150)

w_long <- weights_static |> pivot_longer(-Date, names_to = "Asset", values_to = "Weight") |>
  mutate(Asset = factor(Asset, levels = LABELS,
                        labels = c("PZU", "BNP Paribas", "CD Projekt", "Orange Polska", "Orlen")))
p3 <- ggplot(w_long, aes(Date, Weight, fill = Asset)) + geom_area(alpha = 0.9) +
  scale_fill_brewer(palette = "Set2") + scale_y_continuous(labels = percent) +
  labs(title = "What the minimum-variance portfolio held over time",
       subtitle = "Sample-covariance model, long-only, rebuilt every 5 trading days",
       x = NULL, y = NULL, fill = NULL)
ggsave("figures/mvp_weights.png", p3, width = 9, height = 4.5, dpi = 150)

rr <- overall |> mutate(Type = ifelse(Strategy %in% key, "Portfolio", "Single stock"),
                        Label = nice[Strategy])
p4 <- ggplot(rr, aes(Volatility, CAGR, colour = Type, label = Label)) +
  geom_point(size = 3) + geom_text(hjust = -0.12, size = 3.2, show.legend = FALSE) +
  scale_x_continuous(labels = percent, expand = expansion(mult = c(0.05, 0.25))) +
  scale_y_continuous(labels = percent) +
  scale_colour_manual(values = c(Portfolio = "#1b3a6b", `Single stock` = "grey55")) +
  labs(title = "Risk vs return, 2019-2025", x = "Annual volatility (risk)",
       y = "Average annual return (CAGR)", colour = NULL)
ggsave("figures/risk_return.png", p4, width = 9, height = 5.5, dpi = 150)

yr <- yearly |> filter(Strategy %in% key) |> mutate(Strategy = relabel(Strategy))
p5 <- ggplot(yr, aes(Year, CAGR, fill = Strategy)) + geom_col(position = "dodge") +
  scale_fill_manual(values = pal) + scale_y_continuous(labels = percent) +
  labs(title = "Return in each calendar year", x = NULL, y = NULL, fill = NULL) +
  guides(fill = guide_legend(nrow = 2))
ggsave("figures/yearly_returns.png", p5, width = 9, height = 4.5, dpi = 150)

cat("\nDone. See results/summary.md and figures/.\n")
