################################################################################
#  O2-based growth & respiration analysis – Temperature gradient
#  Now with carbon-unit conversions (fg C h⁻¹ & fg C min⁻¹)
#  Figures styled for ISME, including Carbon Use Efficiency
################################################################################

## ───────────────────────── 1.  Libraries ──────────────────────────────── ##
suppressPackageStartupMessages({
  library(tidyverse)
  library(minpack.lm)
  library(patchwork)
  library(zoo)
})

## ───────────────────────── 2.  QC / Fit thresholds ───────────────────── ##
REL_SE_THRESHOLD <- 0.15
PVAL_THRESHOLD   <- 0.001
R2_THRESHOLD     <- 0.90
MAX_RESID_RANGE  <- 0.20
AIC_IMPROVEMENT  <- 10
MAPE_MAX         <- 0.15

## ───────────────────────── 3.  Fixed inoculum ────────────────────────── ##
N0_fixed <- 6.4e6    # cells mL⁻¹

## ───────────────────────── 4.  Cell-carbon constants ─────────────────── ##
# Rod geometry
cell_width  <- 0.65           # µm
cell_length <- 2.25           # µm
cell_radius <- cell_width / 2
cell_vol_um3 <- pi * cell_radius^2 * (cell_length - cell_width) +    # cylinder
  (4/3) * pi * cell_radius^3                           # caps
C_density_fg_per_um3 <- 100
init_biomass_fgC <- N0_fixed * cell_vol_um3 * C_density_fg_per_um3   # fg C mL⁻¹

# Helper: mg O₂ L⁻¹ → mol O₂ mL⁻¹  (MW O₂ = 32 g)
mgL_to_mol_per_mL <- function(mg_per_L) (mg_per_L * 1e-3) / 32 / 1000

## ───────────────────────── 5.  Load data ─────────────────────────────── ##
oxygen_data <- read_csv("Oxygen_Data_Filtered_CUE.csv") %>%
  mutate(
    Taxon       = as.character(Taxon),
    Temperature = as.numeric(Temperature),
    Replicate   = as.character(Replicate)
  )

## ───────────────────────── 6.  Model function ────────────────────────── ##
resp_with_growth <- function(N0, r, resp_rate, t, O2_0) {
  O2_0 + (resp_rate / r) * N0 * (1 - exp(r * t))
}

## ───────────────────────── 7.  Plot theme ────────────────────────────── ##
isme_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14),
      axis.line  = element_line(size = 0.8),
      axis.ticks = element_line(size = 0.6),
      legend.position = "none"
    )
}

## ───────────────────────── 8.  Prepare outputs ───────────────────────── ##
results <- tibble(
  Taxon = character(), Temperature = numeric(), Replicate = character(),
  N0 = numeric(),
  r_per_minute = numeric(), r_per_hour = numeric(),
  resp_rate = numeric(), O2_0 = numeric(), AICc = numeric(),
  lnO2_change_per_min = numeric(), pseudo_R2 = numeric(), fit_ok = logical(),
  growth_C_fg_per_hr  = numeric(), growth_C_fg_per_min = numeric(),
  resp_C_fg_per_hr    = numeric(), resp_C_fg_per_min = numeric(),
  resp_to_growth_C    = numeric(),
  carbon_use_efficiency = numeric()
)
plots_list <- list()

## ───────────────────────── 9.  Fitting loop ─────────────────────────── ##
grouped <- oxygen_data %>% group_by(Taxon, Temperature, Replicate)
combos  <- group_keys(grouped)

for (i in seq_len(nrow(combos))) {
  Tax  <- combos$Taxon[i]
  Temp <- combos$Temperature[i]
  Rep  <- combos$Replicate[i]
  
  df <- grouped %>%
    filter(Taxon == Tax, Temperature == Temp, Replicate == Rep) %>%
    arrange(Time)
  if (nrow(df) < 5) next
  
  # Compute derivative & rolling mean to find fitting window
  df <- df %>%
    mutate(
      dO2   = c(NA, diff(Oxygen)),
      dt    = c(NA, diff(Time)),
      dO2dt = dO2 / dt,
      sm    = rollmean(dO2dt, 3, fill = NA, align = "right")
    )
  idx <- which(df$sm < -1e-7)[1]
  if (!is.na(idx)) idx <- idx + 15
  if (is.na(idx) || idx > nrow(df) - 2) idx <- 1
  
  df <- df[idx:nrow(df), ] %>% mutate(Time0 = Time - min(Time, na.rm = TRUE))
  
  # Normalize O2
  O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)  # mg O2 L⁻¹
  df <- df %>% mutate(Oxygen_norm = Oxygen / O0)
  
  # Starting values
  r_start <- {
    seg <- head(df, max(3, floor(0.3 * nrow(df))))
    slopes <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
    pmin(pmax(max(slopes, na.rm = TRUE), 1e-4), 5e-2)
  }
  resp_rate_start <- {
    slope <- abs(min(diff(df$Oxygen_norm) / diff(df$Time0), na.rm = TRUE))
    pmin(pmax(slope / N0_fixed, 1e-12), 1e-6)
  }
  
  # Nonlinear fit
  fit <- tryCatch(
    nlsLM(
      Oxygen_norm ~ resp_with_growth(N0_fixed, r, resp_rate, Time0, O2_0),
      data    = df,
      start   = list(r = r_start, resp_rate = resp_rate_start, O2_0 = 1),
      lower   = c(r = 1e-4, resp_rate = 1e-12, O2_0 = 0.8),
      upper   = c(r = 0.1,  resp_rate = 1e-6,  O2_0 = 1.2),
      control = nls.lm.control(maxiter = 300)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  
  # Outlier pruning and re-fit
  pred <- predict(fit, df)
  df_kept <- df[abs(df$Oxygen_norm - pred) < 2 * sd(df$Oxygen_norm - pred, na.rm = TRUE), ]
  fit <- tryCatch(update(fit, data = df_kept), error = function(e) fit)
  
  # Extract parameters
  pars  <- coef(summary(fit))
  r_est <- pars["r", "Estimate"]           # min⁻¹
  C_est <- pars["resp_rate", "Estimate"]   # fraction O₂ cell⁻¹ min⁻¹
  
  pseudo_R2 <- 1 - sum(residuals(fit)^2) /
    sum((df_kept$Oxygen_norm - mean(df_kept$Oxygen_norm))^2)
  
  fit_ok <- all(
    abs(pars[, "Std. Error"] / pars[, "Estimate"]) < REL_SE_THRESHOLD,
    pars[, "Pr(>|t|)"] < PVAL_THRESHOLD,
    pseudo_R2 >= R2_THRESHOLD,
    diff(range(residuals(fit), na.rm = TRUE)) < MAX_RESID_RANGE,
    between(C_est, 1e-12, 1e-6),
    between(r_est, 1e-4, 0.1),
    AIC(lm(Oxygen_norm ~ 1, data = df_kept)) - AIC(fit) >= AIC_IMPROVEMENT,
    mean(abs(residuals(fit) / df_kept$Oxygen_norm)) < MAPE_MAX
  )
  if (!fit_ok) next
  
  # Log-change
  lnchg <- (log(last(df_kept$Oxygen_norm)) - log(first(df_kept$Oxygen_norm))) /
    (last(df_kept$Time0) - first(df_kept$Time0))
  
  # Carbon-unit conversions
  growth_C_fg_per_hr  <- r_est * 60 * init_biomass_fgC
  growth_C_fg_per_min <- growth_C_fg_per_hr / 60
  
  O0_mol_per_mL     <- mgL_to_mol_per_mL(O0)
  resp_molO2_per_hr <- C_est * N0_fixed * O0_mol_per_mL * 60
  resp_C_fg_per_hr  <- resp_molO2_per_hr * 12e15
  resp_C_fg_per_min <- resp_C_fg_per_hr / 60
  
  # Carbon Use Efficiency
  carbon_use_efficiency <- growth_C_fg_per_hr / 
    (growth_C_fg_per_hr + resp_C_fg_per_hr)
  
  resp_to_growth_C <- resp_C_fg_per_hr / growth_C_fg_per_hr
  
  # Store results
  results <- add_row(results,
                     Taxon                  = Tax,
                     Temperature            = Temp,
                     Replicate              = Rep,
                     N0                     = N0_fixed,
                     r_per_minute           = r_est,
                     r_per_hour             = r_est * 60,
                     resp_rate              = C_est,
                     O2_0                   = coef(fit)["O2_0"],
                     AICc                   = AIC(fit),
                     lnO2_change_per_min    = lnchg,
                     pseudo_R2              = pseudo_R2,
                     fit_ok                 = TRUE,
                     growth_C_fg_per_hr     = growth_C_fg_per_hr,
                     growth_C_fg_per_min    = growth_C_fg_per_min,
                     resp_C_fg_per_hr       = resp_C_fg_per_hr,
                     resp_C_fg_per_min      = resp_C_fg_per_min,
                     resp_to_growth_C       = resp_to_growth_C,
                     carbon_use_efficiency  = carbon_use_efficiency
  )
  
  # Plot fit
  df_kept$Pred <- predict(fit, df_kept)
  label <- paste(Tax, Temp, Rep, sep = "_")
  plots_list[[label]] <-
    ggplot(df_kept, aes(Time0, Oxygen_norm)) +
    geom_point(size = 2.5) +
    geom_line(aes(y = Pred), colour = "red", size = 0.8) +
    labs(
      title = label,
      x     = expression(Time ~ "(min)"),
      y     = expression(Normalised ~ O[2])
    ) +
    isme_theme()
}

## ───────────────────────── 10.  Save results & plots ─────────────────── ##
write_csv(results, "oxygen_model_results_good_only.csv")

if (length(plots_list) > 0) {
  pdf("oxygen_dynamics_all_models.pdf", width = 14, height = 10)
  print(wrap_plots(plots_list))
  dev.off()
  
  pdf("oxygen_dynamics_fullsize_per_page.pdf", width = 8, height = 6)
  walk(plots_list, print)
  dev.off()
}

## ───────────────────────── 11.  Temperature-response fits (3 panels) ───────── ##
#   • First try TWO-sided Sharpe–Schoolfield  (6 parameters: E, El, Eh, Tl, Th)
#   • If that fails: one-sided SS  (4 parameters)
#   • If that fails: Arrhenius / Boltzmann
#   • Safe exp() so geom_line() always has >1 point
#   • Output: SharpeSchoolfield_Temperature_Fits.pdf
##################################################################################

# 11.1 Constant ------------------------------------------------------------------
k_B <- 8.617e-5    # Boltzmann constant (eV K⁻¹)

# 11.2 Safe helpers --------------------------------------------------------------
safe_exp <- function(z) exp(pmin(700, z))        # prevents Inf
predict_safe <- function(fit, newdata) {
  out <- try(predict(fit, newdata = newdata), silent = TRUE)
  if (inherits(out, "try-error")) rep(NA_real_, nrow(newdata)) else out
}

# 11.3 Model functions -----------------------------------------------------------
# Two-sided Sharpe–Schoolfield
ln_SS_two <- function(T_C, lnR0, E, El, Tl, Eh, Th) {
  T   <- T_C + 273.15
  TlK <- Tl  + 273.15
  ThK <- Th  + 273.15
  lnR0 - E /(k_B *   T) -
    log1p( safe_exp(El / k_B * (1/T    - 1/TlK) ) ) -
    log1p( safe_exp(Eh / k_B * (1/ThK  - 1/T   ) ) )
}
# One-sided (high-T only)
ln_SS_one <- function(T_C, lnR0, E, Eh, Th) {
  T   <- T_C + 273.15
  ThK <- Th  + 273.15
  lnR0 - E /(k_B * T) -
    log1p( safe_exp(Eh / k_B * (1/ThK - 1/T)) )
}
# Boltzmann / Arrhenius
ln_boltz <- function(T_C, lnR0, E) {
  T <- T_C + 273.15
  lnR0 - E /(k_B * T)
}

# 11.4 Fit-and-plot helper -------------------------------------------------------
fit_T_plot <- function(df, yvar, ylab, show_x = TRUE) {
  library(minpack.lm)
  
  df <- df %>% filter(is.finite(.data[[yvar]]), .data[[yvar]] > 0)
  base <- ggplot(df, aes(Temperature, .data[[yvar]])) +
    geom_point(size = 2.5) +
    labs(x = if (show_x) expression(Temperature~"(°C)") else NULL,
         y = ylab) +
    isme_theme()
  if (nrow(df) < 4) return(base)
  
  T_min <- min(df$Temperature); T_max <- max(df$Temperature)
  T_opt <- df$Temperature[which.max(df[[yvar]])]
  
  # ---------- 2-sided SS --------------------------------------------------------
  start_two <- list(
    lnR0 = log(max(df[[yvar]])),
    E  = 0.6,
    El = 0.4,  Tl = T_min + 2,    # low-T branch
    Eh = 1.5,  Th = T_opt + 3     # high-T branch
  )
  fit_two <- try(
    nlsLM(
      as.formula(
        paste0("log(", yvar,
               ") ~ ln_SS_two(Temperature, lnR0, E, El, Tl, Eh, Th)")
      ),
      data = df, start = start_two,
      lower = c(lnR0=-Inf, E=0.1, El=0.1, Tl= 0,  Eh=0.1, Th= 0),
      upper = c(lnR0= Inf, E=2.5, El=2.5, Tl=60, Eh=5.0, Th=60),
      control = nls.lm.control(maxiter = 600)
    ),
    silent = TRUE
  )
  
  # ---------- 1-sided SS --------------------------------------------------------
  start_one <- list(
    lnR0 = log(max(df[[yvar]])),
    E  = 0.6,
    Eh = if (yvar == "growth_C_fg_per_hr") 3.0 else 1.5,
    Th = if (yvar == "growth_C_fg_per_hr") T_opt + 2 else T_opt + 5
  )
  fit_one <- try(
    nlsLM(
      as.formula(
        paste0("log(", yvar,
               ") ~ ln_SS_one(Temperature, lnR0, E, Eh, Th)")
      ),
      data = df, start = start_one,
      lower = c(lnR0=-Inf, E=0.1, Eh=0.5, Th=0),
      upper = c(lnR0= Inf, E=2.5, Eh=6.0, Th=60),
      control = nls.lm.control(maxiter = 500)
    ),
    silent = TRUE
  )
  
  # ---------- Boltzmann fallback -----------------------------------------------
  fit_bol <- try(
    nlsLM(
      as.formula(
        paste0("log(", yvar,
               ") ~ ln_boltz(Temperature, lnR0, E)")
      ),
      data  = df,
      start = list(lnR0 = log(max(df[[yvar]])), E = 0.6),
      lower = c(lnR0=-Inf, E=0.1),
      upper = c(lnR0= Inf,  E=2.5),
      control = nls.lm.control(maxiter = 400)
    ),
    silent = TRUE
  )
  
  # ---------- choose best model (lowest AIC) ------------------------------------
  get_AIC <- function(f) if (inherits(f, "try-error")) Inf else AIC(f)
  best_fit <- list(two = fit_two, one = fit_one, bol = fit_bol) |>
    purrr::imap(~ list(f = .x, a = get_AIC(.x))) |>
    purrr::compact() |>
    purrr::reduce(~ if (.x$a <= .y$a) .x else .y)
  if (best_fit$a == Inf) {
    message("⚠ ", yvar, ": no model converged.")
    return(base)
  }
  
  fit <- best_fit$f
  grid <- tibble(Temperature = seq(T_min, T_max, length.out = 300))
  grid$Pred <- safe_exp(predict_safe(fit, grid))
  grid <- grid %>% filter(is.finite(Pred))
  if (nrow(grid) > 1)
    base <- base + geom_line(data = grid, aes(Temperature, Pred), linewidth = 1)
  
  model_lab <- names(best_fit)[1]
  message("✅ ", yvar, ": using ", model_lab, "-sided model.")
  base
}

# 11.5 Prepare data --------------------------------------------------------------
results_filtered <- results %>%
  filter(fit_ok,
         growth_C_fg_per_hr > 0,
         resp_C_fg_per_hr   > 0)

# 11.6 Build panels --------------------------------------------------------------
p_growth <- fit_T_plot(results_filtered, "growth_C_fg_per_hr",
                       expression(Growth~(fg~C~h^{-1})), TRUE)
p_resp   <- fit_T_plot(results_filtered, "resp_C_fg_per_hr",
                       expression(Respiration~(fg~C~h^{-1})), TRUE)
p_cue    <- fit_T_plot(results_filtered, "carbon_use_efficiency",
                       "Carbon Use Efficiency", TRUE)

# 11.7 Combine & export ----------------------------------------------------------
combo_plot <- wrap_plots(p_growth, p_resp, p_cue, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.margin = margin(5, 10, 5, 10))

ggsave("SharpeSchoolfield_Temperature_Fits.pdf",
       combo_plot, width = 12, height = 4, dpi = 600)

##################################################################################
#                                   End Section 11                               #
##################################################################################

