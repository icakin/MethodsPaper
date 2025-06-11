# -------------------- LOAD LIBRARIES --------------------
library(readr)
library(dplyr)
library(ggplot2)
library(minpack.lm)
library(patchwork)
library(tidyr)
library(lme4)
library(multcomp)
library(ggtext)       # still used for theme element_markdown
library(gridExtra)

# -------------------- LOAD DATA --------------------
oxygen_data <- read_csv("Raw_Oxygen.csv")
fc_data     <- read_csv("Cell_Counts.csv")

# Convert events/µL → cells/L
fc_data <- fc_data %>%
  mutate(
    FC_Initial = FC_Initial * 1e6,
    FC_Final   = FC_Final   * 1e6
  )

# Merge FC_Initial into oxygen_data
merged_data <- oxygen_data %>%
  inner_join(fc_data, by = c("Taxon", "Replicate"))

# -------------------- MODEL FUNCTION --------------------
resp_with_growth <- function(N0, r, resp_rate, t, O2_0, t0) {
  O2 <- O2_0 + (resp_rate / r) * N0 * (1 - exp(r * (t - t0)))
  return(O2)
}

# -------------------- INITIALIZE CONTAINERS --------------------
results     <- data.frame()
plots_list  <- list()
all_fits_df <- data.frame()

# -------------------- FIT NONLINEAR MODEL PER Taxon × Replicate --------------------
grouped    <- merged_data %>% group_by(Taxon, Replicate)
group_keys <- grouped %>% group_keys()

for (i in seq_len(nrow(group_keys))) {
  this_taxon <- group_keys$Taxon[i]
  this_repl  <- group_keys$Replicate[i]
  group_data <- grouped %>%
    filter(Taxon == this_taxon, Replicate == this_repl) %>%
    ungroup()
  
  N0 <- unique(group_data$FC_Initial)
  label <- paste(this_taxon, this_repl, sep = "_")
  cat("Processing:", label, "\n")
  
  if (length(N0) != 1 || is.na(N0) || nrow(group_data) < 5) {
    message(" Skipping: ", label, " (invalid N0 or <5 rows)")
    next
  }
  
  O2_0_start      <- max(group_data$Oxygen, na.rm = TRUE)
  resp_rate_start <- abs(mean(diff(group_data$Oxygen), na.rm = TRUE)) / mean(N0)
  r_start         <- 0.001
  
  fit <- tryCatch({
    nlsLM(
      Oxygen ~ resp_with_growth(
        N0 = N0,
        r,
        resp_rate,
        t   = Time,
        O2_0,
        t0  = 0
      ),
      data = group_data,
      start = list(r = r_start, resp_rate = resp_rate_start, O2_0 = O2_0_start),
      lower = c(r = 0.0001, resp_rate = 1e-12, O2_0 = 1),
      upper = c(r = 0.05,    resp_rate = 1e-6,  O2_0 = 10)
    )
  }, error = function(e) {
    message(" Fit failed for ", label, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(fit)) {
    pars <- summary(fit)$parameters
    r_per_minute <- pars["r", "Estimate"]
    r_per_hour   <- r_per_minute * 60
    
    O2_initial <- group_data$Oxygen[1]
    O2_final   <- group_data$Oxygen[nrow(group_data)]
    delta_t    <- group_data$Time[nrow(group_data)] - group_data$Time[1]
    lnO2_change_per_min <- if (O2_initial > 0 && O2_final > 0 && delta_t > 0) {
      (log(O2_final) - log(O2_initial)) / delta_t
    } else {
      NA
    }
    
    group_data$Predicted_O2 <- predict(fit, newdata = group_data)
    y_obs <- group_data$Oxygen
    y_pred <- group_data$Predicted_O2
    SS_res <- sum((y_obs - y_pred)^2)
    SS_tot <- sum((y_obs - mean(y_obs))^2)
    pseudo_R2 <- if (SS_tot > 0) 1 - (SS_res / SS_tot) else NA
    
    results <- rbind(
      results,
      data.frame(
        Taxon               = this_taxon,
        Replicate           = this_repl,
        N0                  = N0,
        r_per_minute        = r_per_minute,
        r_per_hour          = r_per_hour,
        resp_rate           = pars["resp_rate", "Estimate"],
        O2_0                = pars["O2_0", "Estimate"],
        AICc                = AIC(fit),
        lnO2_change_per_min = lnO2_change_per_min,
        pseudo_R2           = pseudo_R2,
        stringsAsFactors    = FALSE
      )
    )
    
    all_fits_df <- bind_rows(all_fits_df, group_data)
    
    p <- ggplot(group_data, aes(x = Time)) +
      geom_point(aes(y = Oxygen), color = "blue", size = 2) +
      geom_line(aes(y = Predicted_O2), color = "red", linewidth = 1) +
      labs(title = label, x = "Time (minutes)", y = "Oxygen (mg/L)") +
      theme_classic()
    
    plots_list[[label]] <- p
  }
}

# -------------------- SAVE PLOTS --------------------
if (length(plots_list) > 0) {
  pdf("oxygen_dynamics_all_replicates.pdf", width = 14, height = 10)
  print(wrap_plots(plots_list))
  dev.off()
  
  pdf("oxygen_dynamics_fullsize_per_page.pdf", width = 8, height = 6)
  for (nm in names(plots_list)) print(plots_list[[nm]])
  dev.off()
}

# -------------------- SAVE RESULTS TABLE --------------------
results <- results %>% arrange(Taxon, Replicate)
write_csv(results, "oxygen_model_results.csv")
print(results)


# -------------------- 4-PANEL FACET PLOT --------------------
selected_combos_4 <- tibble::tribble(
  ~Taxon,               ~Replicate,
  "Acinetobacter",       "R2",
  "Aeromonas strain B",  "R1",
  "Arthrobacter",        "R2",
  "Yersinia",            "R1"
)
letters_vec_4 <- letters[1:nrow(selected_combos_4)]

facet_data_4 <- all_fits_df %>%
  inner_join(selected_combos_4, by = c("Taxon", "Replicate")) %>%
  left_join(results, by = c("Taxon", "Replicate")) %>%
  mutate(
    FacetLabel = paste0(
      "(", letters_vec_4[
        match(
          paste(Taxon, Replicate),
          paste(selected_combos_4$Taxon, selected_combos_4$Replicate)
        )
      ], ")~italic('", Taxon, "')"
    ),
    label_text = sprintf("r = %.3f\nR = %.2e", r_per_minute, resp_rate)
  )

annotations_4 <- facet_data_4 %>% distinct(FacetLabel, label_text)

facet_plot_4 <- ggplot(facet_data_4, aes(x = Time)) +
  geom_point(aes(y = Oxygen), color = "black", size = 1.8, alpha = 0.8) +
  geom_line(aes(y = Predicted_O2), color = "black", linewidth = 1) +
  facet_wrap(~ FacetLabel, nrow = 1, labeller = label_parsed) +
  geom_text(
    data        = annotations_4,
    aes(x = Inf, y = Inf, label = label_text),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.1,
    size = 4
  ) +
  labs(
    x = "Time (minutes)",
    y = expression(paste("Oxygen (mg ", L^{-1}, ")"))
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 12, face = "italic"),
    axis.title       = element_text(size = 13),
    axis.text        = element_text(size = 11),
    panel.spacing    = unit(1, "lines"),
    plot.margin      = margin(5, 5, 5, 5),
    axis.line        = element_line(size = 0.6),
    axis.ticks       = element_line(size = 0.5)
  )

ggsave(
  filename = "oxygen_dynamics_facet4_row_custom.pdf",
  plot     = facet_plot_4,
  width    = 14,
  height   = 4.5
)

# Prepare data for supplementary plot
all_fits_df <- merged_data %>%
  inner_join(results, by = c("Taxon", "Replicate")) %>%
  group_by(Taxon, Replicate) %>%
  mutate(
    Predicted_O2 = resp_with_growth(
      N0 = unique(FC_Initial),
      r = unique(r_per_minute),
      resp_rate = unique(resp_rate),
      t = Time,
      O2_0 = unique(O2_0),
      t0 = 0
    )
  ) %>%
  ungroup()

# Supplementary figure: All replicates by Taxon, colored lines
supp_plot <- ggplot(all_fits_df, aes(x = Time, y = Oxygen)) +
  geom_point(color = "grey50", size = 1, alpha = 0.5) +
  geom_line(aes(y = Predicted_O2, color = Replicate), linewidth = 0.8) +
  facet_wrap(~ Taxon, scales = "free_y") +
  labs(
    title = "Supplementary Figure: Oxygen Dynamics Across All Replicates",
    x = "Time (minutes)",
    y = expression(paste("Oxygen (mg ", L^{-1}, ")"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Save supplementary figure
ggsave(
  filename = paste0( "supp_oxygen_all_replicates_facet_by_taxon.pdf"),
  plot = supp_plot,
  width = 14,
  height = 10
)

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#### ------------ DERIVE RATES USING FC & OD600 DATA ------------ ####

# Arrange merged data by Taxon, Replicate, and Time
merged_data <- merged_data %>%
  arrange(Taxon, Replicate, Time)

# Calculate experiment duration and keep only last time point
OD_FC <- merged_data %>%
  group_by(Taxon, Replicate) %>%
  mutate(Duration = max(Time) - min(Time) - 10) %>%  # Subtract 10 minutes assumed for initial sampling
  filter(Time == max(Time)) %>%
  ungroup()

# Calculate growth rates from OD600 and Flow Cytometry
OD_FC <- OD_FC %>%
  mutate(
    r_OD600 = (log(OD_Final) - log(OD_Initial)) / Duration,
    r_FC = (log(FC_Final) - log(FC_Initial)) / Duration
  ) %>%
  dplyr::select(Taxon, Replicate, Duration, r_OD600, r_FC)

#### ------------ BRING GROWTH ESTIMATES TOGETHER ------------ ####

growth <- results %>%
  dplyr::select(Taxon, Replicate, r_per_minute) %>%
  rename(r_O2 = r_per_minute) %>%
  left_join(OD_FC, by = c("Taxon", "Replicate")) %>%
  dplyr::select(Taxon, Replicate, r_O2, r_OD600, r_FC)

# Derive doubling times from each growth rate estimate
growth <- growth %>%
  mutate(
    doubling_time_O2 = log(2) / r_O2,
    doubling_time_OD600 = log(2) / r_OD600,
    doubling_time_FC = log(2) / r_FC
  )

# ------------ PLOT: RESHAPE AND VISUALIZE ------------

# Colorblind-friendly palette (Okabe–Ito)
cb_colors <- c(
  "Oxygen"         = "#E69F00",  # orange
  "OD600"          = "#56B4E9",  # sky blue
  "Flow Cytometry" = "#009E73"   # bluish green
)

# Reshape data for plotting
growth_long <- growth %>%
  dplyr::select(Taxon, Replicate, r_O2, r_OD600, r_FC) %>%
  pivot_longer(
    cols = c(r_O2, r_OD600, r_FC),
    names_to = "Method",
    values_to = "Growth_Rate"
  ) %>%
  mutate(Method = factor(Method, 
                         levels = c("r_O2", "r_OD600", "r_FC"),
                         labels = c("Oxygen", "OD600", "Flow Cytometry")))

# Create boxplot by taxon
growth_comparison <- ggplot(growth_long, aes(x = Taxon, y = Growth_Rate, fill = Method)) +
  geom_vline(xintercept = seq(1.5, length(unique(growth_long$Taxon)) - 0.5), 
             color = "gray80", linetype = "dashed") +
  geom_boxplot() +
  scale_fill_manual(values = cb_colors) +
  theme_classic() +
  labs(x = "Taxon", 
       y = "Growth Rate (per minute)", 
       fill = "Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  facet_grid(. ~ "By Taxon")

# Create global comparison boxplot
global_comparison <- ggplot(growth_long, aes(x = Method, y = Growth_Rate, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = cb_colors) +
  theme_classic() +
  labs(x = "Method", 
       y = "Growth Rate (per minute)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  facet_grid(. ~ "All Taxa")

# Combine both plots using patchwork
combined_comparison <- growth_comparison + global_comparison +
  plot_layout(widths = c(3, 1))

# Save the combined plot
ggsave(
  filename = "growth_rate_comparison.pdf", 
  plot     = combined_comparison, 
  width    = 14, 
  height   = 6
)

#### ------------ ANOVA MIXED EFFECTS MODEL ------------ ####

mm_method <- lmer(Growth_Rate ~ 1 + Method + (1|Taxon), REML = FALSE, data = growth_long)
mm_method_1 <- lmer(Growth_Rate ~ 1 + (1|Taxon), REML = FALSE, data = growth_long)
AIC(mm_method, mm_method_1)

summary(mm_method)
confint(mm_method)

# Calculate confidence intervals
ci_method <- confint(mm_method, method="Wald")
coefficients <- fixef(mm_method)

# Perform multiple comparisons with Bonferroni correction
mc <- glht(mm_method, linfct = mcp(Method = "Tukey"), alternative = "two.sided",
           test = adjusted("bonferroni"))
summary_mc <- summary(mc)
mc_ci <- confint(mc)

# Create dataframe for plotting
method_effects <- data.frame(
  Method = c("Oxygen", "OD600", "Flow Cytometry"),
  Estimate = c(coefficients[1], 
               coefficients[2] + coefficients[1],
               coefficients[3] + coefficients[1]),
  CI_lower = c(ci_method["(Intercept)", 1],
               ci_method["(Intercept)", 1] + ci_method["MethodOD600", 1],
               ci_method["(Intercept)", 1] + ci_method["MethodFlow Cytometry", 1]),
  CI_upper = c(ci_method["(Intercept)", 2],
               ci_method["(Intercept)", 2] + ci_method["MethodOD600", 2],
               ci_method["(Intercept)", 2] + ci_method["MethodFlow Cytometry", 2])
)

# Get significance for plotting with stars
pairs <- data.frame(
  y.position = max(method_effects$CI_upper) + c(0.0005, 0.001, 0.0015),
  xmin = c(1, 1, 2),
  xmax = c(2, 3, 3),
  p.value = summary_mc$test$pvalues) %>%
  mutate(stars = case_when(
    p.value > 0.05 ~ "ns",
    p.value <= 0.05 & p.value > 0.01 ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001 ~ "***"
  ))

# Create confidence interval plot with significance bars
method_comparison <- ggplot(method_effects, aes(x = Method, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  theme_classic() +
  geom_segment(data = pairs, aes(x = xmin, xend = xmax, 
                                 y = y.position, yend = y.position)) +
  geom_text(data = pairs, aes(x = (xmin + xmax)/2, y = y.position,
                              label = stars),
            vjust = -0.5) +
  labs(x = "Method", 
       y = "Growth Rate (per minute)",
       title = "Comparison of Method Effects with 95% CIs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave(paste0("method_effects_CI.pdf"), 
       method_comparison, 
       width = 8, 
       height = 6)

# Print the estimates, CIs and pairwise comparisons
print(method_effects)
print(summary_mc)

library(lme4)
library(ggplot2)
library(dplyr)
library(patchwork)

#### ------------ LINEAR MIXED EFFECTS MODEL ------------ ####

mm_OD <- lmer(r_O2 ~ r_OD600 + (1 | Taxon), data = growth)
mm_FC <- lmer(r_O2 ~ r_FC + (1 | Taxon), data = growth)

### Normalized plots for both OD600 and FC

# Function to create normalized data
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
  
  r2_mixed <- cor(fitted(model), growth$r_O2)^2
  
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

# Create normalized data for both methods
od_norm <- create_norm_data(mm_OD, "r_OD600")
fc_norm <- create_norm_data(mm_FC, "r_FC")

# Find common axis range for both plots
rate_range_od <- range(c(od_norm$data$r_O2_norm, od_norm$data$r_OD600), na.rm = TRUE)
rate_range_fc <- range(c(fc_norm$data$r_O2_norm, fc_norm$data$r_FC), na.rm = TRUE)
rate_range <- range(c(rate_range_od, rate_range_fc))

# Extended Okabe-Ito color palette for 13 taxa
okabe_ito_extended <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7", "#000000",
  "#E69F00B0", "#56B4E9B0", "#009E73B0", "#0072B2B0", "#D55E00B0"
)

# OD600 plot
p_od <- ggplot(od_norm$data, aes(x = r_OD600, y = r_O2_norm)) +
  geom_point(aes(color = Taxon), size = 3) +
  geom_abline(intercept = od_norm$fixed_coef[1], slope = od_norm$fixed_coef[2], 
              color = "black", linetype = "dashed") +
  scale_color_manual(values = okabe_ito_extended) +
  theme_classic() +
  labs(x = expression(paste("Growth Rate - OD"[600], " (min"^-1, ")")),
       y = expression(paste("Growth Rate - O"[2], " (min"^-1, ")"))) +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  coord_cartesian(xlim = rate_range, ylim = rate_range) +
  annotate("text", x = min(rate_range) + 0.1 * (diff(rate_range)), 
           y = max(rate_range) - 0.1 * (diff(rate_range)),
           label = od_norm$eq_text, hjust = 0, vjust = 1) +
  ggtitle("(A)")

# FC plot
p_fc <- ggplot(fc_norm$data, aes(x = r_FC, y = r_O2_norm)) +
  geom_point(aes(color = Taxon), size = 3) +
  geom_abline(intercept = fc_norm$fixed_coef[1], slope = fc_norm$fixed_coef[2], 
              color = "black", linetype = "dashed") +
  scale_color_manual(values = okabe_ito_extended) +
  theme_classic() +
  labs(x = expression(paste("Growth Rate - FC (min"^-1, ")")),
       y = expression(paste("Growth Rate - O"[2], " (min"^-1, ")"))) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = rate_range, ylim = rate_range) +
  annotate("text", x = min(rate_range) + 0.1 * (diff(rate_range)), 
           y = max(rate_range) - 0.1 * (diff(rate_range)),
           label = fc_norm$eq_text, hjust = 0, vjust = 1) +
  ggtitle("(B)")

# Combine plots vertically
combined_norm <- p_od / p_fc

# Save the combined plot
ggsave("growth_rate_regression_normalized_combined.pdf",
       combined_norm,
       width = 10,
       height = 16)


# -------------------- LOAD LIBRARIES --------------------
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)

# -------------------- LOAD DATA --------------------
df <- read_csv("OD_r_FC_r.csv")

# -------------------- CALCULATE GROWTH RATES --------------------
df <- df %>%
  mutate(
    OD_r = (log(OD_Final) - log(OD_Initial)) / Time,
    FC_r = (log(FC_Final) - log(FC_Initial)) / Time
  )

# -------------------- NATURE/ISME-LIKE THEME --------------------
theme_nature <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 0.8, color = "black"),
      plot.title = element_text(size = base_size + 2, face = "bold"),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      legend.position = "none",
      panel.border = element_blank()
    )
}

# -------------------- BLAND–ALTMAN FUNCTION --------------------
bland_altman_plot <- function(df, method1, method2, label1, label2, panel_letter) {
  df <- df %>%
    mutate(
      avg = ({{ method1 }} + {{ method2 }}) / 2,
      diff = {{ method1 }} - {{ method2 }}
    )
  
  bias <- mean(df$diff, na.rm = TRUE)
  sd_diff <- sd(df$diff, na.rm = TRUE)
  upper <- bias + 1.96 * sd_diff
  lower <- bias - 1.96 * sd_diff
  within_limits <- mean(df$diff >= lower & df$diff <= upper, na.rm = TRUE) * 100
  
  cat(paste0("\n", label1, " vs ", label2, ":\n"))
  cat(sprintf("Bias: %.5f\n", bias))
  cat(sprintf("Upper LoA: %.5f\n", upper))
  cat(sprintf("Lower LoA: %.5f\n", lower))
  cat(sprintf("%% within limits: %.1f%%\n", within_limits))
  
  ggplot(df, aes(x = avg, y = diff)) +
    # Fully filled orange LoA area
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = lower, ymax = upper,
             fill = "#D55E00", alpha = 0.15) +
    geom_point(size = 2.5, alpha = 0.9, color = "black") +
    geom_hline(yintercept = bias, color = "black", size = 0.7) +
    geom_hline(yintercept = upper, linetype = "dashed", color = "#D55E00", size = 0.8) +
    geom_hline(yintercept = lower, linetype = "dashed", color = "#D55E00", size = 0.8) +
    annotate("text", x = min(df$avg, na.rm = TRUE), y = upper,
             label = sprintf("Upper LoA = %.3f", upper),
             hjust = 0, vjust = -0.8, size = 4, color = "#D55E00", fontface = "bold") +
    annotate("text", x = min(df$avg, na.rm = TRUE), y = lower,
             label = sprintf("Lower LoA = %.3f", lower),
             hjust = 0, vjust = 1.8, size = 4, color = "#D55E00", fontface = "bold") +
    annotate("text", x = min(df$avg, na.rm = TRUE), y = bias,
             label = sprintf("Bias = %.3f", bias),
             hjust = 0, vjust = -1.3, size = 4, fontface = "bold") +
    annotate("text", x = min(df$avg, na.rm = TRUE), y = lower - 0.005,
             label = sprintf("%% within limits: %.1f%%", within_limits),
             hjust = 0, size = 4, fontface = "italic") +
    labs(
      title = paste0("(", panel_letter, ") ", label1, " vs ", label2),
      x = "Mean growth rate",
      y = "Difference in growth rate"
    ) +
    theme_nature()
}

# -------------------- GENERATE BOTH PLOTS --------------------
plot1 <- bland_altman_plot(df, Oxygen_r, OD_r, "Oxygen_r", "OD_r", "A")
plot2 <- bland_altman_plot(df, Oxygen_r, FC_r, "Oxygen_r", "FC_r", "B")

# -------------------- COMBINE INTO ONE FIGURE --------------------
combined_plot <- grid.arrange(plot1, plot2, ncol = 2)

# -------------------- SAVE TO PDF --------------------
ggsave("BlandAltman_SingleFigure_ISMEstyle.pdf", combined_plot, width = 14, height = 6)
