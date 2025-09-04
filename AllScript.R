# ── Libraries ─────────────────────────────────────────────────────────
library(tidyverse)
library(zoo)

# ── Tunables (trimming aggressiveness) ────────────────────────────────
SPLINE_SPAR    <- 0.4     # higher = smoother (often trims earlier)
RUN_LEN        <- 10       # consecutive calm points required for plateau
REL_PROP       <- 0.008    # fraction of max post-peak down-slope allowed
ABS_THR        <- 0.0003  # absolute slope threshold floor
WINDOW_LEN     <- 4       # rolling range window length for flatness
LEVEL_DELTA    <- 0.0025  # allowed vertical wiggle (range) within window
FLAT_RANGE_OK  <- 0.05    # skip series if total range below this

INPUT_CSV      <- "Oxygen_Data_Long.csv"
TRIMMED_CSV    <- "Oxygen_Data_Trimmed.csv"
FILTERED_CSV   <- "Oxygen_Data_Filtered.csv"
SKIPPED_CSV    <- "Skipped_Series_Log.csv"
DIAG_PDF       <- "oxygen_trimming_diagnostics.pdf"

# ── Helper functions ──────────────────────────────────────────────────
find_plateau <- function(o2_vec, peak_idx,
                         run_len = 2, rel_prop = 0.01, abs_thr = 0.0003,
                         window_len = 4, level_delta = 0.0025) {
  n <- length(o2_vec); if (peak_idx >= n) return(n)
  slopes <- diff(o2_vec)
  max_down <- abs(min(slopes[peak_idx:(n - 1)], na.rm = TRUE))
  slope_thr <- max(abs_thr, rel_prop * max_down)
  
  rng <- rep(NA_real_, n)
  if (n >= window_len) {
    rng[1:(n - window_len + 1)] <-
      rollapply(o2_vec, window_len,
                FUN = function(x) max(x) - min(x), align = "left")
  }
  
  for (i in seq(peak_idx, n - run_len)) {
    if (all(abs(slopes[i:(i + run_len - 1)]) <= slope_thr, na.rm = TRUE) &&
        all(rng[i:(i + run_len - 1)] <= level_delta, na.rm = TRUE)) {
      return(i)
    }
  }
  n
}

find_second_inflection <- function(o2_fit, idx_peak, idx_plate) {
  slopes <- diff(o2_fit); curves <- diff(slopes)
  post_pk <- seq(idx_peak + 1, length(curves))
  if (length(post_pk) == 0) return(NA_integer_)
  steepest <- which.min(slopes[post_pk]) + idx_peak
  win <- seq(steepest + 1, idx_plate - 1)
  win <- win[!is.na(curves[win])]
  if (length(win) < 3) return(NA_integer_)
  win[which.max(curves[win])]
}

# ── Load data ─────────────────────────────────────────────────────────
raw <- read_csv(INPUT_CSV, show_col_types = FALSE) %>%
  mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate),
    series_id = paste(Taxon, "Rep=", Replicate, sep = " ")
  )

series_ids  <- unique(raw$series_id)
trimmed_lst <- vector("list", length(series_ids))
names(trimmed_lst) <- series_ids
skipped_log <- tibble(series_id = character(), reason = character())

# ── Main loop ─────────────────────────────────────────────────────────
for (sid in series_ids) {
  df <- raw %>% filter(series_id == sid) %>% arrange(Time)
  
  if (nrow(df) < 5) {
    skipped_log <- add_row(skipped_log, series_id = sid, reason = "Too few rows")
    next
  }
  
  # Smoothing used only for detection; not plotted
  df <- df %>%
    mutate(O2_fit = suppressWarnings(
      predict(smooth.spline(Time, Oxygen, spar = SPLINE_SPAR), x = Time)$y
    ))
  
  if (all(is.na(df$O2_fit))) {
    skipped_log <- add_row(skipped_log, series_id = sid, reason = "Spline failed")
    next
  }
  
  if ((max(df$O2_fit, na.rm = TRUE) - min(df$O2_fit, na.rm = TRUE)) < FLAT_RANGE_OK) {
    skipped_log <- add_row(skipped_log, series_id = sid, reason = "Too flat")
    next
  }
  
  # Peak on raw Oxygen (true peak)
  idx_peak <- which.max(df$Oxygen)
  
  # Plateau & second inflection on smoothed curve
  idx_plate  <- min(
    find_plateau(df$O2_fit, idx_peak,
                 run_len = RUN_LEN, rel_prop = REL_PROP, abs_thr = ABS_THR,
                 window_len = WINDOW_LEN, level_delta = LEVEL_DELTA),
    nrow(df)
  )
  idx_second <- find_second_inflection(df$O2_fit, idx_peak, idx_plate)
  
  if (is.na(idx_second) || idx_second <= idx_peak || idx_second > nrow(df)) {
    idx_end <- nrow(df); used_full <- TRUE
  } else {
    idx_end <- idx_second; used_full <- FALSE
  }
  
  trimmed_lst[[sid]] <- df[idx_peak:idx_end, ] %>%
    mutate(
      peak_time           = Time[idx_peak],
      plateau_time        = if (idx_plate <= nrow(df)) Time[idx_plate] else NA_real_,
      second_inflect_time = if (!is.na(idx_second)) Time[idx_second] else NA_real_,
      used_full_series    = used_full,
      series_id           = sid
    )
}

# ── Save CSV outputs ──────────────────────────────────────────────────
trimmed <- bind_rows(trimmed_lst)
write_csv(trimmed,  TRIMMED_CSV)

filtered <- trimmed %>%
  select(Taxon, Replicate, Time, Oxygen)
write_csv(filtered, FILTERED_CSV)

write_csv(skipped_log, SKIPPED_CSV)

# ── Diagnostics PDF (no spline or marker overlays) ────────────────────
pdf(DIAG_PDF, 7, 5)

for (sid in names(trimmed_lst)) {
  df_trim <- trimmed_lst[[sid]]
  if (is.null(df_trim) || nrow(df_trim) == 0) next
  
  raw_df <- raw %>%
    filter(series_id == sid) %>%
    arrange(Time)
  
  # Robust bounds for the orange rectangle
  xmin_val <- suppressWarnings(min(df_trim$Time, na.rm = TRUE))
  xmax_val <- suppressWarnings(max(df_trim$Time, na.rm = TRUE))
  if (!is.finite(xmin_val) || !is.finite(xmax_val)) next
  
  rect_df <- tibble(
    xmin = xmin_val,
    xmax = xmax_val,
    ymin = -Inf, ymax = Inf
  )
  
  p <- ggplot(raw_df, aes(Time, Oxygen)) +
    geom_rect(
      data = rect_df, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "orange", alpha = 0.08
    ) +
    geom_line(colour = "grey60", linewidth = 0.8) +
    labs(
      title = sid,
      subtitle = "Trimmed zone (orange)",
      x = "Time (min)", y = "O₂ (mg L⁻¹)"
    ) +
    theme_classic(base_size = 11)
  
  print(p)
}

dev.off()


################################################################################
#  O2‐based growth & respiration analysis
#  FULLY UPDATED: N0 from Cell_Counts.csv, no cell‐carbon constants
################################################################################

## ───────────────────────── 1.  Libraries ───────────────────────────────── ##
suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, ggplot2, purrr, readr, etc.
  library(minpack.lm)  # Levenberg–Marquardt nls
  library(patchwork)   # for wrap_plots()
  library(zoo)         # rolling means
})

## ──────────────────── 2.  QC / Fit thresholds ───────────────────────────── ##
REL_SE_THRESHOLD <- 0.15
PVAL_THRESHOLD   <- 0.001
R2_THRESHOLD     <- 0.90
MAX_RESID_RANGE  <- 0.20
AIC_IMPROVEMENT  <- 10
MAPE_MAX         <- 0.15

## ───────────────────────── 3.  Load & merge N0 ───────────────────────────── ##
# Read FC counts and convert events/µL → cells/L
fc_data <- read_csv("Cell_Counts.csv") %>%
  mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate),
    N0        = FC_Initial * 1e6  # cells/L
  ) %>%
  select(Taxon, Replicate, N0)

## ───────────────────────── 4.  Load oxygen data ──────────────────────────── ##
oxygen_data <- read_csv("Oxygen_Data_Filtered.csv") %>%
  mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  ) %>%
  inner_join(fc_data, by = c("Taxon", "Replicate"))

## ───────────────────────── 5.  Model function ────────────────────────────── ##
resp_with_growth <- function(N0, r, resp_rate, t, O2_0) {
  O2_0 + (resp_rate / r) * N0 * (1 - exp(r * t))
}

## ───────────────────────── 6.  Plot theme ───────────────────────────────── ##
isme_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14),
      axis.line  = element_line(linewidth = .8),
      axis.ticks = element_line(linewidth = .6),
      legend.position = "none"
    )
}

## ───────────────────────── 7.  Prepare outputs ──────────────────────────── ##
results    <- tibble(
  Taxon               = character(),
  Replicate           = character(),
  N0                  = numeric(),
  r_per_minute        = numeric(),
  r_per_hour          = numeric(),
  resp_rate           = numeric(),
  O2_0                = numeric(),
  AICc                = numeric(),
  lnO2_change_per_min = numeric(),
  pseudo_R2           = numeric(),
  fit_ok              = logical()
)
plots_list <- list()

## ───────────────────────── 8.  Fitting loop ─────────────────────────────── ##
combos <- group_keys(group_by(oxygen_data, Taxon, Replicate))

for (i in seq_len(nrow(combos))) {
  cfg <- combos[i, ]
  Tax <- cfg$Taxon
  Rep <- cfg$Replicate
  
  df <- oxygen_data %>%
    filter(Taxon == Tax, Replicate == Rep) %>%
    arrange(Time)
  if (nrow(df) < 5) next
  
  # detect onset of decline via 3‐point rolling mean of dO2/dt
  df <- df %>%
    mutate(
      dO2   = c(NA, diff(Oxygen)),
      dt    = c(NA, diff(Time)),
      dO2dt = dO2 / dt,
      sm    = rollmean(dO2dt, 3, fill = NA, align = "right")
    )
  idx <- which(df$sm < -1e-7)[1]
  idx <- if (is.na(idx) || idx > nrow(df) - 2) 1 else idx + 5
  df  <- df[idx:nrow(df), ] %>%
    mutate(Time0 = Time - min(Time, na.rm = TRUE))
  
  # normalize O2
  O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)
  df <- df %>% mutate(Oxygen_norm = Oxygen / O0)
  
  # starting values
  r_start <- {
    seg    <- head(df, max(3, floor(.3 * nrow(df))))
    slopes <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
    pmin(pmax(max(slopes, na.rm = TRUE), 1e-4), 5e-2)
  }
  resp_rate_start <- {
    slope <- abs(min(diff(df$Oxygen_norm) / diff(df$Time0), na.rm = TRUE))
    pmin(pmax(slope / unique(df$N0), 1e-12), 1e-6)
  }
  
  fit <- tryCatch(
    nlsLM(
      Oxygen_norm ~ resp_with_growth(N0, r, resp_rate, Time0, O2_0),
      data    = df,
      start   = list(r = r_start, resp_rate = resp_rate_start, O2_0 = 1),
      lower   = c(r = 1e-4, resp_rate = 1e-12, O2_0 = .8),
      upper   = c(r = .1,   resp_rate = 1e-6,  O2_0 = 1.2),
      control = nls.lm.control(maxiter = 300)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  
  # prune outliers and refit
  pred    <- predict(fit, df)
  keep_ix <- abs(df$Oxygen_norm - pred) < 2 * sd(df$Oxygen_norm - pred, na.rm = TRUE)
  df_kept <- df[keep_ix, ]
  fit     <- tryCatch(update(fit, data = df_kept), error = function(e) fit)
  
  pars      <- coef(summary(fit))
  r_est     <- pars["r", "Estimate"]
  C_est     <- pars["resp_rate", "Estimate"]
  pseudo_R2 <- 1 - sum(residuals(fit)^2) /
    sum((df_kept$Oxygen_norm - mean(df_kept$Oxygen_norm))^2)
  fit_ok    <- all(
    abs(pars[,"Std. Error"] / pars[,"Estimate"]) < REL_SE_THRESHOLD,
    pars[,"Pr(>|t|)"] < PVAL_THRESHOLD,
    pseudo_R2 >= R2_THRESHOLD,
    diff(range(residuals(fit), na.rm = TRUE)) < MAX_RESID_RANGE,
    between(C_est, 1e-12, 1e-6),
    between(r_est, 1e-4, 0.1),
    AIC(lm(Oxygen_norm ~ 1, data = df_kept)) - AIC(fit) >= AIC_IMPROVEMENT,
    mean(abs(residuals(fit) / df_kept$Oxygen_norm)) < MAPE_MAX
  )
  if (!fit_ok) next
  
  # collect results
  results <- results %>%
    add_row(
      Taxon               = Tax,
      Replicate           = Rep,
      N0                  = unique(df$N0),
      r_per_minute        = r_est,
      r_per_hour          = r_est * 60,
      resp_rate           = C_est,
      O2_0                = coef(fit)["O2_0"],
      AICc                = AIC(fit),
      lnO2_change_per_min = (log(last(df_kept$Oxygen_norm)) -
                               log(first(df_kept$Oxygen_norm))) /
        (last(df_kept$Time0) - first(df_kept$Time0)),
      pseudo_R2           = pseudo_R2,
      fit_ok              = TRUE
    )
  
  # diagnostic plot
  df_kept <- df_kept %>% mutate(Pred = predict(fit, df_kept))
  plot_key <- paste(Tax, Rep, sep = "_")
  plots_list[[plot_key]] <-
    ggplot(df_kept, aes(Time0, Oxygen_norm)) +
    geom_point(size = 2.5) +
    geom_line(aes(y = Pred), linewidth = .8, colour = "red") +
    labs(
      title = plot_key,
      x     = expression(Time~"(min)"),
      y     = expression(Normalised~O[2])
    ) +
    isme_theme()
}

## ──────────────────── 9.  Save outputs ─────────────────────────── ##
write_csv(results, "oxygen_model_results.csv")

if (length(plots_list) > 0) {
  pdf("oxygen_dynamics_all_models.pdf", width = 14, height = 10)
  print(wrap_plots(plots_list))
  dev.off()
  
  pdf("oxygen_dynamics_fullsize_per_page.pdf", width = 8, height = 6)
  walk(plots_list, print)
  dev.off()
}


################################################################################
#  4-panel facet plot (labels under curve, no overflow, fixed sci notation)
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

# ========== TUNABLES ==========================================================
TEXT_SIZE    <- 3            # annotation font size
X_ANCHOR     <- 0.00001      # relative x-position (0..1 from left)
X_MARGIN     <- 0.05         # min horizontal margin fraction
Y_GAP_FRAC   <- 0.80         # vertical gap below curve (fraction of y-range)
Y_MARGIN     <- 0.05         # min vertical margin fraction
ONE_LINE     <- FALSE        # one-line labels (TRUE) vs stacked atop() (FALSE)

# ========== helper: scientific notation for plotmath ==========================
sci_pm <- function(x, digits = 2) {
  out <- rep(NA_character_, length(x))
  ok  <- is.finite(x)
  if (any(ok)) {
    s <- formatC(x[ok], format = "e", digits = digits - 1)   # e.g. "2.4e-11"
    m <- sub("e[+-].*$", "", s)                              # "2.4"
    e <- sub("^.*e([+-]?)([0-9]+)$", "\\1\\2", s)            # "-11"
    e <- sub("^\\+", "", e)
    out[ok] <- sprintf("%s%%*%%10^{%s}", m, e)               # "2.4%*%10^{-11}"
  }
  out
}

# ========== choose 4 panels ===================================================
selected_combos <- tribble(
  ~Taxon,          ~Replicate,
  "Bacillus", "R4",
  "Burkholderia",   "R2",
  "Arthrobacter",  "R2",
  "Yersinia",      "R1"
)
letters_vec <- letters[seq_len(nrow(selected_combos))]

# ========== rebuild raw + fit data from plots_list ============================
facet_data <- map2_dfr(
  selected_combos$Taxon, selected_combos$Replicate,
  ~{
    key <- paste(.x, .y, sep = "_")
    p   <- plots_list[[key]]
    if (is.null(p)) return(NULL)
    
    pts <- layer_data(p, 1) %>% select(x, y) %>% rename(Time = x, Oxygen = y)
    fit <- layer_data(p, 2) %>% select(x, y) %>% rename(Time = x, Predicted_O2 = y)
    
    full_join(pts, fit, by = "Time") %>%
      mutate(Taxon = .x, Replicate = .y, series_id = key)
  }
)
stopifnot(all(c("Time","Oxygen","Predicted_O2","Taxon","Replicate") %in% names(facet_data)))

# ========== build annotation info (safe math labels) ==========================
annot_info <- results %>%
  filter(fit_ok) %>%
  inner_join(selected_combos, by = c("Taxon","Replicate")) %>%
  mutate(
    FacetLabel = paste0(
      "(", letters_vec[
        match(paste(Taxon, Replicate),
              paste(selected_combos$Taxon, selected_combos$Replicate))
      ], ")~italic('", Taxon, "')"
    ),
    label_text = if (ONE_LINE) {
      paste0(
        "italic(r)==", sprintf("%.2f", r_per_hour), "~h^{-1}*','~~",
        "italic(R)==", sci_pm(resp_rate, digits = 2),
        "~mg~O[2]~cell^{-1}~min^{-1}"
      )
    } else {
      paste0(
        "atop(",
        "italic(r)==", sprintf("%.2f", r_per_hour), "~h^{-1},",
        "italic(R)==", sci_pm(resp_rate, digits = 2),
        "~mg~O[2]~cell^{-1}~min^{-1})"
      )
    }
  ) %>%
  select(Taxon, Replicate, FacetLabel, label_text)

facet_data <- facet_data %>%
  inner_join(annot_info, by = c("Taxon","Replicate"))

# ========== compute label positions (under curve, clamped) ====================
label_positions <- facet_data %>%
  group_by(FacetLabel) %>%
  group_modify(~{
    df <- arrange(.x, Time) %>% distinct(Time, .keep_all = TRUE)
    
    xmin <- min(df$Time,    na.rm = TRUE); xmax <- max(df$Time,    na.rm = TRUE)
    ymin <- min(df$Oxygen,  na.rm = TRUE); ymax <- max(df$Oxygen,  na.rm = TRUE)
    xr   <- xmax - xmin;                 yr   <- ymax - ymin
    
    # horizontal: anchor very left but keep at least X_MARGIN inside
    x_pos <- xmin + X_ANCHOR * xr
    x_pos <- max(xmin + X_MARGIN * xr, min(xmax - X_MARGIN * xr, x_pos))
    
    # curve value at x_pos
    ok <- is.finite(df$Predicted_O2)
    curveY <- if (sum(ok) >= 2) approx(df$Time[ok], df$Predicted_O2[ok],
                                       xout = x_pos, rule = 2)$y
    else ymin + 0.90 * yr
    
    # vertical: large gap below curve, then clamp to margins
    gap   <- Y_GAP_FRAC * yr
    y_raw <- curveY - gap
    y_pos <- max(ymin + Y_MARGIN * yr, min(ymax - Y_MARGIN * yr, y_raw))
    
    tibble(x_pos = x_pos, y_pos = y_pos)
  }) %>%
  ungroup()

annotations <- facet_data %>%
  distinct(FacetLabel, label_text) %>%
  left_join(label_positions, by = "FacetLabel")

# ========== build plot ========================================================
facet_plot <- ggplot(facet_data, aes(x = Time)) +
  geom_point(aes(y = Oxygen), size = 2.2, alpha = 0.85) +
  geom_line(aes(y = Predicted_O2), linewidth = 1.2) +
  facet_wrap(~ FacetLabel, nrow = 1, labeller = label_parsed) +
  geom_text(
    data = annotations,
    aes(x = x_pos, y = y_pos, label = label_text),
    inherit.aes = FALSE,
    parse = TRUE, hjust = 0, vjust = 1.05, size = TEXT_SIZE
  ) +
  coord_cartesian(clip = "on") +  # no overflow outside panel
  labs(
    x = "Time (minutes)",
    y = expression("Normalised Oxygen (" * O[2] / O[2*","*0] * ")")
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 16, face = "italic"),
    axis.title       = element_text(size = 17),
    axis.text        = element_text(size = 14),
    panel.spacing    = unit(1.2, "lines"),
    plot.margin      = margin(12, 20, 12, 20),
    axis.line        = element_line(linewidth = 0.7),
    axis.ticks       = element_line(linewidth = 0.6)
  )


ggsave("oxygen_dynamics_facet4_no_overflow.pdf", facet_plot, width = 14, height = 4.5)

################################################################################
# Supplementary (NORMALISED): color by Replicate (R1–R5), facet by TaxonFull
# Supports plot keys:
#   "Flavobacterium_A_R1"  or  "Acinetobacter_R2"
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(forcats)
})

all_fits_df_norm <- imap_dfr(
  plots_list,
  ~{
    key <- .y
    p   <- .x
    
    # Expect: layer 1 = points, layer 2 = fitted line
    pts <- tryCatch(layer_data(p, 1) %>% select(x, y) %>% rename(Time = x, Oxygen = y),
                    error = function(e) NULL)
    fit <- tryCatch(layer_data(p, 2) %>% select(x, y) %>% rename(Time = x, Predicted = y),
                    error = function(e) NULL)
    if (is.null(pts) || is.null(fit)) return(NULL)
    
    df <- full_join(pts, fit, by = "Time") %>% arrange(Time)
    
    # Parse "<TaxonFull>_<Replicate>" by splitting at the LAST underscore
    m1 <- str_match(key, "^(.*)_(R[0-9]+)$")
    TaxonFull <- m1[,2]
    Replicate <- m1[,3]
    
    # Defensive normalisation (match per-page figs)
    O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)
    if (is.finite(O0) && O0 > 0) {
      df <- df %>%
        mutate(
          Oxygen_n = Oxygen   / O0,
          Pred_n   = Predicted / O0
        )
    } else {
      df <- df %>% mutate(Oxygen_n = Oxygen, Pred_n = Predicted)
    }
    
    tibble(
      TaxonFull = TaxonFull,
      Replicate = Replicate,
      Time      = df$Time,
      Oxygen_n  = df$Oxygen_n,
      Pred_n    = df$Pred_n
    )
  }
) %>%
  filter(!is.na(TaxonFull), !is.na(Replicate)) %>%
  mutate(
    # Order legend as R1..R5 (keeps others if present)
    Replicate = factor(Replicate, levels = paste0("R", 1:5))
  )

# Plot: grey points, colored fit lines by Replicate
supp_plot_norm_rep <- ggplot(all_fits_df_norm, aes(x = Time, y = Oxygen_n)) +
  geom_point(color = "grey60", size = 1, alpha = 0.55) +
  geom_line(
    data = all_fits_df_norm %>% arrange(TaxonFull, Replicate, Time),
    aes(y = Pred_n, color = Replicate, group = Replicate),
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  facet_wrap(~ TaxonFull, scales = "free_y") +
  labs(
    title = "Supplementary: Normalised O₂ Dynamics (colour = Replicate)",
    x = "Time (minutes)",
    y = expression("Normalised Oxygen ("*O[2]*"/"*O[2][0]*")"),
    color = "Replicate"
  ) +
  scale_color_brewer(palette = "Dark2", drop = FALSE) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )

ggsave(
  filename = "supp_oxygen_all_replicates_NORMALISED_by_replicate.pdf",
  plot     = supp_plot_norm_rep,
  width    = 14,
  height   = 10
)


################################################################################
# BOXPLOT with ns lines on right panel
################################################################################

# ────────────────────── Load Required Packages ──────────────────────
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)
library(ggsignif)   # << for significance bars

# ─────── Read Growth Durations and Initial/Final Counts ───────
od_fc_data <- read_csv("OD_r_FC_r.csv", show_col_types = FALSE) %>%
  mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate),
    Duration  = Time
  )

# ─────── Calculate Growth Rates from OD600 and FC ───────
OD_FC <- od_fc_data %>%
  mutate(
    r_OD600 = (log(OD_Final) - log(OD_Initial)) / Duration,
    r_FC    = (log(FC_Final) - log(FC_Initial)) / Duration
  ) %>%
  dplyr::select(Taxon, Replicate, Duration, r_OD600, r_FC)

# ─────── Combine with O₂-based Growth Rates ───────
growth <- results %>%
  dplyr::select(Taxon, Replicate, r_per_minute) %>%
  rename(r_O2 = r_per_minute) %>%
  left_join(OD_FC, by = c("Taxon", "Replicate")) %>%
  dplyr::select(Taxon, Replicate, r_O2, r_OD600, r_FC)

# ─────── Calculate Doubling Times ───────
growth <- growth %>%
  mutate(
    doubling_time_O2    = log(2) / r_O2,
    doubling_time_OD600 = log(2) / r_OD600,
    doubling_time_FC    = log(2) / r_FC
  )

write_csv(growth, "growth_rates_combined.csv")

# ─────── Plot: Growth Rate Comparison ───────
cb_colors <- c(
  "Oxygen"         = "#E69F00",
  "OD600"          = "#56B4E9",
  "Flow Cytometry" = "#009E73"
)

# 1. Reshape for plotting
growth_long <- growth %>%
  pivot_longer(
    cols = c(r_O2, r_OD600, r_FC),
    names_to = "Method",
    values_to = "Growth_Rate"
  ) %>%
  mutate(
    Method = factor(
      Method,
      levels = c("r_O2", "r_OD600", "r_FC"),
      labels = c("Oxygen", "OD600", "Flow Cytometry")
    )
  )

# 2. Boxplot by Taxon
growth_comparison <- ggplot(growth_long, aes(x = Taxon, y = Growth_Rate, fill = Method)) +
  geom_vline(
    xintercept = seq(1.5, length(unique(growth_long$Taxon)) - 0.5),
    color      = "gray80", linetype = "dashed"
  ) +
  geom_boxplot(outlier.size = 1.2, outlier.shape = 16) +
  scale_fill_manual(values = cb_colors) +
  theme_classic(base_size = 16) +
  labs(
    x    = "Taxon",
    y    = expression(paste("Growth rate (min"^{-1},")")),
    fill = "Method"
  ) +
  theme(
    axis.title      = element_text(size = 18),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 14, face = "italic"),
    axis.text.y     = element_text(size = 14),
    legend.title    = element_text(size = 16),
    legend.text     = element_text(size = 14),
    legend.position = "top",
    strip.text      = element_text(size = 16, face = "italic")
  ) +
  facet_grid(. ~ "By Taxon")

# 3. Global Method Comparison + ns lines
global_comparison <- ggplot(growth_long, aes(x = Method, y = Growth_Rate, fill = Method)) +
  geom_boxplot(outlier.size = 1.2, outlier.shape = 16) +
  geom_signif(
    comparisons = list(
      c("Oxygen", "OD600"),
      c("Oxygen", "Flow Cytometry"),
      c("OD600", "Flow Cytometry")
    ),
    annotations = c("ns", "ns", "ns"),   # force "ns"
    step_increase = 0.1,
    tip_length = 0.01
  ) +
  scale_fill_manual(values = cb_colors) +
  theme_classic(base_size = 16) +
  labs(
    x = "Method",
    y = expression(paste("Growth rate (min"^{-1},")"))
  ) +
  theme(
    axis.title  = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "plain"),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  facet_grid(. ~ "All Taxa")

# 4. Combine + Save Plot
combined_comparison <- growth_comparison + global_comparison + plot_layout(widths = c(3, 1))

ggsave(
  filename = "growth_rate_comparison.pdf",
  plot     = combined_comparison,
  width    = 14,
  height   = 6
)


# ─────────────────────────── ANOVA: MIXED EFFECTS ────────────────────────────

library(lme4)
library(multcomp)

# Fit mixed effects models: method as fixed effect, taxon as random intercept
mm_method   <- lmer(Growth_Rate ~ 1 + Method + (1 | Taxon), REML = FALSE, data = growth_long)
mm_method_1 <- lmer(Growth_Rate ~ 1 + (1 | Taxon), REML = FALSE, data = growth_long)

anova(mm_method, mm_method_1)

# Compare model AICs
aic_comparison <- AIC(mm_method, mm_method_1)
print(aic_comparison)

# Model summary and CIs
print(summary(mm_method))
ci_method <- confint(mm_method, method = "Wald")
coefficients <- fixef(mm_method)

# Post hoc pairwise comparisons with Bonferroni correction
mc      <- glht(mm_method, linfct = mcp(Method = "Tukey"), test = adjusted("bonferroni"))
summary_mc <- summary(mc)
mc_ci      <- confint(mc)

# Build dataframe for plot: estimates and CIs
method_effects <- data.frame(
  Method = c("Oxygen", "OD600", "Flow Cytometry"),
  Estimate = c(coefficients[1], 
               coefficients[1] + coefficients["MethodOD600"],
               coefficients[1] + coefficients["MethodFlow Cytometry"]),
  CI_lower = c(ci_method["(Intercept)", 1],
               ci_method["(Intercept)", 1] + ci_method["MethodOD600", 1],
               ci_method["(Intercept)", 1] + ci_method["MethodFlow Cytometry", 1]),
  CI_upper = c(ci_method["(Intercept)", 2],
               ci_method["(Intercept)", 2] + ci_method["MethodOD600", 2],
               ci_method["(Intercept)", 2] + ci_method["MethodFlow Cytometry", 2])
)

# Significance stars from pairwise p-values
pairs <- data.frame(
  y.position = max(method_effects$CI_upper, na.rm = TRUE) + c(0.0005, 0.001, 0.0015),
  xmin = c(1, 1, 2),
  xmax = c(2, 3, 3),
  p.value = summary_mc$test$pvalues
) %>%
  mutate(stars = case_when(
    p.value > 0.05                     ~ "ns",
    p.value <= 0.05 & p.value > 0.01  ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001                  ~ "***"
  ))

# Plot: Estimated effects with confidence intervals
method_comparison <- ggplot(method_effects, aes(x = Method, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_segment(data = pairs, aes(x = xmin, xend = xmax, y = y.position, yend = y.position)) +
  geom_text(data = pairs, aes(x = (xmin + xmax) / 2, y = y.position, label = stars),
            vjust = -0.5, size = 5) +
  labs(
    x = "Method", 
    y = expression("Growth rate ("*min^{-1}*")"),
    title = "Comparison of Method Effects (Mixed Model)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the CI plot
ggsave("method_effects_CI.pdf", method_comparison, width = 8, height = 6)

# Also export key results
write_csv(method_effects, "method_effects_estimates.csv")
write_csv(pairs, "method_effects_significance.csv")

# Print summaries to console
print(method_effects)
print(summary_mc)


# ───────────────────────────── LME: O₂ vs OD600 / FC ─────────────────────────────

library(lme4)
library(ggplot2)
library(patchwork)
library(dplyr)

# Fit models: O2-based growth ~ OD600 / FC + (1 | Taxon)
mm_OD <- lmer(r_O2 ~ r_OD600 + (1 | Taxon), data = growth)
mm_FC <- lmer(r_O2 ~ r_FC + (1 | Taxon), data = growth)

# Helper function to normalize by random intercepts
create_norm_data <- function(model, x_var) {
  ranef_taxon <- ranef(model)$Taxon
  fixed_coef <- fixef(model)
  
  rand_effects <- data.frame(
    Taxon = rownames(ranef_taxon),
    rand_int = ranef_taxon[, 1]
  )
  
  growth_norm <- growth %>%
    left_join(rand_effects, by = "Taxon") %>%
    mutate(r_O2_norm = r_O2 - rand_int)
  
  r2_mixed <- cor(fitted(model), growth$r_O2, use = "complete.obs")^2
  
  eq_text <- sprintf("y = %.3f + %.3fx\nR² = %.3f", 
                     fixed_coef[1], 
                     fixed_coef[2],
                     r2_mixed)
  
  list(
    data = growth_norm,
    fixed_coef = fixed_coef,
    eq_text = eq_text,
    x_var = x_var
  )
}

# Normalize and extract equations
od_norm <- create_norm_data(mm_OD, "r_OD600")
fc_norm <- create_norm_data(mm_FC, "r_FC")

# Shared axis limits
rate_range_od <- range(c(od_norm$data$r_O2_norm, od_norm$data$r_OD600), na.rm = TRUE)
rate_range_fc <- range(c(fc_norm$data$r_O2_norm, fc_norm$data$r_FC), na.rm = TRUE)
rate_range <- range(c(rate_range_od, rate_range_fc))

# Okabe–Ito extended palette (15 colors)
okabe_ito_extended <- c(
  "#E69F00",   # orange
  "#56B4E9",   # sky blue
  "#009E73",   # bluish green
  "#F0E442",   # yellow
  "#0072B2",   # blue
  "#D55E00",   # vermillion
  "#CC79A7",   # reddish purple
  "#000000",   # black
  "#E69F00B0", # transparent orange
  "#56B4E9B0", # transparent sky blue
  "#009E73B0", # transparent bluish green
  "#0072B2B0", # transparent blue
  "#D55E00B0", # transparent vermillion
  "#999999",   # soft gray
  "#AA4499"    # purple-pink
)

# OD600 plot
p_od <- ggplot(od_norm$data, aes(x = r_OD600, y = r_O2_norm)) +
  geom_point(aes(color = Taxon), size = 3) +
  # 1:1 line (black dashed)
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  # model regression line (solid black)
  geom_abline(intercept = od_norm$fixed_coef[1], slope = od_norm$fixed_coef[2],
              color = "black", linetype = "solid") +
  scale_color_manual(values = okabe_ito_extended) +
  theme_classic() +
  labs(
    x = expression(paste("Growth Rate - OD"[600], " (min"^-1, ")")),
    y = expression(paste("Growth Rate - O"[2], " (min"^-1, ")"))
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) +
  coord_cartesian(xlim = rate_range, ylim = rate_range) +
  annotate("text",
           x = min(rate_range) + 0.1 * diff(rate_range),
           y = max(rate_range) - 0.1 * diff(rate_range),
           label = od_norm$eq_text, hjust = 0, vjust = 1, size = 5) +
  ggtitle("(A)")

# FC plot
p_fc <- ggplot(fc_norm$data, aes(x = r_FC, y = r_O2_norm)) +
  geom_point(aes(color = Taxon), size = 3) +
  # 1:1 line (black dashed)
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  # model regression line (solid black)
  geom_abline(intercept = fc_norm$fixed_coef[1], slope = fc_norm$fixed_coef[2],
              color = "black", linetype = "solid") +
  scale_color_manual(values = okabe_ito_extended) +
  theme_classic() +
  labs(
    x = expression(paste("Growth Rate - FC (min"^-1, ")")),
    y = expression(paste("Growth Rate - O"[2], " (min"^-1, ")"))
  ) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = rate_range, ylim = rate_range) +
  annotate("text",
           x = min(rate_range) + 0.1 * diff(rate_range),
           y = max(rate_range) - 0.1 * diff(rate_range),
           label = fc_norm$eq_text, hjust = 0, vjust = 1, size = 5) +
  ggtitle("(B)")

# Combine and save
combined_norm <- p_od / p_fc

ggsave(
  filename = "growth_rate_regression_normalized_combined.pdf",
  plot     = combined_norm,
  width    = 10,
  height   = 16
)

# Also export model summaries
capture.output(summary(mm_OD), file = "mixed_model_OD600_summary.txt")
capture.output(summary(mm_FC), file = "mixed_model_FC_summary.txt")


################################################################################
# Bland–Altman Plots (All Replicates) with LME Regression Inside LoA
# - Uses all replicates (no per-taxon averaging)
# - LoA computed on all points
# - Regression inside LoA: lmer(diff ~ avg + (1 | Taxon))  (stats printed to console only)
# - No regression line, no R²/slope/p/formula in the plots
################################################################################

install.packages(c("lme4","lmerTest"))

# -------- Packages --------
library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)
library(lme4)
library(lmerTest)  # provides p-values for lmer
library(glue)

# -------- Load data --------
# Expecting columns: Taxon, r_O2, r_OD600, r_FC at replicate level
df_raw <- read_csv("growth_rates_combined.csv", show_col_types = FALSE) %>%
  rename(
    Oxygen_r = r_O2,
    OD_r     = r_OD600,
    FC_r     = r_FC
  )

if (!"Taxon" %in% names(df_raw)) stop("Input must contain a 'Taxon' column.")
df_raw <- df_raw %>% mutate(Taxon = as.factor(Taxon))

# -------- Helper: fallback-safe extractors --------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && is.finite(a)) a else b

get_lmer_slope_p <- function(fit) {
  # lmerTest adds p-values to coef summary table
  tab <- coef(summary(fit))
  if (is.null(tab)) return(NA_real_)
  rn <- rownames(tab)
  if (is.null(rn)) return(NA_real_)
  idx <- which(rn == "avg")
  if (!length(idx)) return(NA_real_)
  val <- tab[idx, "Pr(>|t|)"]
  if (length(val) == 0) return(NA_real_) else as.numeric(val)
}

# -------- Bland–Altman with LME --------
bland_altman_plot_lme <- function(df, method1, method2, label1, label2, panel_letter) {
  # Replicate-level BA data
  df_ba <- df %>%
    mutate(
      method1 = {{ method1 }},
      method2 = {{ method2 }},
      avg  = (method1 + method2) / 2,
      diff = method1 - method2
    ) %>%
    filter(is.finite(avg), is.finite(diff))
  
  # Bias & LoA across all points
  bias    <- mean(df_ba$diff, na.rm = TRUE)
  sd_diff <- sd(df_ba$diff,   na.rm = TRUE)
  upper   <- bias + 1.96 * sd_diff
  lower   <- bias - 1.96 * sd_diff
  within_limits <- mean(df_ba$diff >= lower & df_ba$diff <= upper, na.rm = TRUE) * 100
  
  # Restrict regression to within LoA (for stats ONLY; not plotted)
  df_reg <- df_ba %>% filter(diff >= lower, diff <= upper)
  
  # Choose model (fallback to lm if needed)
  use_lmer <- dplyr::n_distinct(df_reg$Taxon) >= 2 && nrow(df_reg) >= 5
  
  if (use_lmer) {
    mm_fit <- lmer(diff ~ avg + (1 | Taxon), data = df_reg, REML = TRUE)
    fe <- fixef(mm_fit)
    intercept <- unname(fe[1]) %||% NA_real_
    slope     <- unname(fe[2]) %||% NA_real_
    r2_mixed  <- cor(fitted(mm_fit), df_reg$diff, use = "complete.obs")^2
    p_value   <- get_lmer_slope_p(mm_fit)
    model_used <- "lmer(diff ~ avg + (1 | Taxon))"
  } else {
    mm_fit <- lm(diff ~ avg, data = df_reg)
    coefs <- coef(mm_fit)
    intercept <- unname(coefs[1]) %||% NA_real_
    slope     <- unname(coefs[2]) %||% NA_real_
    r2_mixed  <- cor(fitted(mm_fit), df_reg$diff, use = "complete.obs")^2
    slope_row <- summary(mm_fit)$coefficients
    p_value   <- if (!is.null(slope_row) && nrow(slope_row) >= 2) slope_row["avg", "Pr(>|t|)"] else NA_real_
    model_used <- "lm(diff ~ avg)"
  }
  
  # ---- Console output of regression + LoA stats ----
  cat(glue(
    "\n[{panel_letter}] Bland–Altman {label1} vs {label2}\n",
    "  Model: {model_used}\n",
    "  N (all points): {nrow(df_ba)},  N (within LoA for regression): {nrow(df_reg)}\n",
    "  Bias: {round(bias, 6)}\n",
    "  SD(diff): {round(sd_diff, 6)}\n",
    "  Upper LoA: {round(upper, 6)}\n",
    "  Lower LoA: {round(lower, 6)}\n",
    "  % within LoA: {round(within_limits, 1)}%\n",
    "  Intercept: {round(intercept, 6)}\n",
    "  Slope: {round(slope, 6)}\n",
    "  R^2 (fitted vs. observed diffs, within LoA): {round(r2_mixed, 6)}\n",
    "  p-value for slope: {ifelse(is.na(p_value), 'NA', format.pval(p_value, digits = 4))}\n"
  ))
  
  # -------- Plot (no regression elements or stats) --------
  ggplot(df_ba, aes(x = avg, y = diff)) +
    # LoA band
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper,
             fill = "#D55E00", alpha = 0.15) +
    # All points (all replicates)
    geom_point(size = 2.5, alpha = 0.9, color = "black") +
    # Reference lines: bias & LoA
    geom_hline(yintercept = bias,  color = "black",    linewidth = 0.7) +
    geom_hline(yintercept = upper, color = "firebrick", linewidth = 0.8, linetype = "dashed") +
    geom_hline(yintercept = lower, color = "firebrick", linewidth = 0.8, linetype = "dashed") +
    # Labels (LoA + bias + % within limits only)
    annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = upper,
             label = sprintf("Upper LoA = %.3f", upper),
             hjust = 0, vjust = -0.8, size = 4, color = "black", fontface = "bold") +
    annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = lower,
             label = sprintf("Lower LoA = %.3f", lower),
             hjust = 0, vjust = 1.8, size = 4, color = "black", fontface = "bold") +
    annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = bias,
             label = sprintf("Bias = %.3f", bias),
             hjust = 0, vjust = -2.6, size = 4, color = "black", fontface = "bold") +
    annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = lower - 0.005,
             label = sprintf("%% within limits: %.1f%%", within_limits),
             hjust = 0, size = 4, color = "black", fontface = "italic") +
    labs(
      title = glue("({panel_letter}) {label1} vs {label2} — Bland–Altman (all replicates)"),
      x = "Mean growth rate (per replicate)",
      y = "Difference in growth rate (per replicate)"
    ) +
    theme_classic(base_size = 14)
}

# -------- Generate & Save --------
plot1 <- bland_altman_plot_lme(df_raw, Oxygen_r, OD_r, "Oxygen", "OD600", "A")
plot2 <- bland_altman_plot_lme(df_raw, Oxygen_r, FC_r, "Oxygen", "Flow Cytometry", "B")

combined_plot <- grid.arrange(plot1, plot2, ncol = 2)
ggsave("BlandAltman_AllReplicates_LME.pdf", combined_plot, width = 14, height = 6)

