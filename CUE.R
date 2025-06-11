################################################################################
#  O2-based growth & respiration analysis – Temperature gradient
#  Figures styled for ISME / Nature / Science, including Carbon Use Efficiency
#  Full script with updated CUE fit (centre-shifted Arrhenius)
################################################################################

## ───────────────────────── 1.  Libraries ──────────────────────────────── ##
suppressPackageStartupMessages({
  library(tidyverse)   # includes dplyr, ggplot2, purrr, readr …
  library(minpack.lm)  # Levenberg-Marquardt nls
  library(patchwork)   # multi-panel layout
  library(zoo)         # rolling means
})

## ───────────────────────── 2.  QC / Fit thresholds ───────────────────── ##
REL_SE_THRESHOLD <- 0.15
PVAL_THRESHOLD   <- 0.001
R2_THRESHOLD     <- 0.90
MAX_RESID_RANGE  <- 0.20
AIC_IMPROVEMENT  <- 10
MAPE_MAX         <- 0.15

## ───────────────────────── 3.  Fixed inoculum ────────────────────────── ##
N0_fixed <- 6.4e6    # cells mL⁻¹ (starting density)

## ───────────────────────── 4.  Cell-carbon constants ─────────────────── ##
# Rod geometry
a <- 0.65 ; b <- 2.25                # width & length (µm)
r <- a/2
cell_vol_um3 <- pi*r^2*(b - a) + 4/3*pi*r^3
C_density_fg_um3 <- 100              # fg C µm⁻³ (≈0.1 pg)
init_biomass_fgC <- N0_fixed * cell_vol_um3 * C_density_fg_um3

# Conversion helper
mgL_to_mol_per_mL <- function(x_mg_L) (x_mg_L * 1e-3) / 32 / 1000   # mg O2 L-1 -> mol O2 mL-1

## ───────────────────────── 5.  Load data ─────────────────────────────── ##
oxygen_data <- read_csv("Oxygen_Data_Filtered.csv") %>%
  mutate(across(c(Taxon, Replicate), as.character),
         Temperature = as.numeric(Temperature))

## ───────────────────────── 6.  Model function ────────────────────────── ##
resp_with_growth <- function(N0, r, resp_rate, t, O2_0) {
  O2_0 + (resp_rate/r)*N0*(1 - exp(r*t))
}

## ───────────────────────── 7.  Plot theme ────────────────────────────── ##
isme_theme <- function() {
  theme_classic(base_size = 14) +
    theme(axis.title = element_text(size = 16),
          axis.text  = element_text(size = 14),
          axis.line  = element_line(linewidth = .8),
          axis.ticks = element_line(linewidth = .6),
          legend.position = "none")
}

## ───────────────────────── 8.  Prepare outputs ───────────────────────── ##
results <- tibble(
  Taxon=character(), Temperature=numeric(), Replicate=character(),
  N0=numeric(), r_per_minute=numeric(), r_per_hour=numeric(),
  resp_rate=numeric(), O2_0=numeric(), AICc=numeric(),
  lnO2_change_per_min=numeric(), pseudo_R2=numeric(), fit_ok=logical(),
  growth_C_fg_per_hr=numeric(), growth_C_fg_per_min=numeric(),
  resp_C_fg_per_hr=numeric(), resp_C_fg_per_min=numeric(),
  resp_to_growth_C=numeric(), carbon_use_efficiency=numeric())
plots_list <- list()

## ───────────────────────── 9.  Fitting loop ─────────────────────────── ##
combos <- group_keys(group_by(oxygen_data, Taxon, Temperature, Replicate))
for (i in seq_len(nrow(combos))) {
  cfg  <- combos[i, ]
  Tax  <- cfg$Taxon; Temp <- cfg$Temperature; Rep <- cfg$Replicate
  
  df <- oxygen_data %>%
    filter(Taxon==Tax, Temperature==Temp, Replicate==Rep) %>%
    arrange(Time)
  if(nrow(df) < 5) next
  # compute smoothed derivative to detect onset of decline
  df <- df %>%
    mutate(dO2=c(NA, diff(Oxygen)), dt=c(NA, diff(Time)), dO2dt=dO2/dt,
           sm=rollmean(dO2dt,3,fill=NA,align="right"))
  idx <- which(df$sm < -1e-7)[1]
  idx <- if(is.na(idx) || idx>nrow(df)-2) 1 else idx+15
  df  <- df[idx:nrow(df),] %>% mutate(Time0=Time-min(Time,na.rm=TRUE))
  # normalise O2
  O0 <- mean(head(df$Oxygen,3), na.rm=TRUE)
  df <- df %>% mutate(Oxygen_norm = Oxygen/O0)
  # start values
  r_start <- {
    seg <- head(df, max(3, floor(.3*nrow(df))))
    slopes <- abs(diff(log(pmax(seg$Oxygen_norm,1e-6)))/diff(seg$Time0))
    pmin(pmax(max(slopes, na.rm=TRUE),1e-4),5e-2)
  }
  resp_rate_start <- {
    slope <- abs(min(diff(df$Oxygen_norm)/diff(df$Time0), na.rm=TRUE))
    pmin(pmax(slope/N0_fixed,1e-12),1e-6)
  }
  fit <- tryCatch(nlsLM(Oxygen_norm ~ resp_with_growth(N0_fixed,r,resp_rate,Time0,O2_0),
                        data=df,
                        start=list(r=r_start, resp_rate=resp_rate_start, O2_0=1),
                        lower=c(r=1e-4,resp_rate=1e-12,O2_0=.8),
                        upper=c(r=.1,  resp_rate=1e-6, O2_0=1.2),
                        control=nls.lm.control(maxiter=300)),
                  error=function(e) NULL)
  if(is.null(fit)) next
  # prune outliers and refit
  pred <- predict(fit, df)
  df_kept <- df[abs(df$Oxygen_norm-pred) < 2*sd(df$Oxygen_norm-pred,na.rm=TRUE),]
  fit <- tryCatch(update(fit, data=df_kept), error=function(e) fit)
  pars <- coef(summary(fit))
  r_est <- pars["r","Estimate"]; C_est <- pars["resp_rate","Estimate"]
  pseudo_R2 <- 1 - sum(residuals(fit)^2)/sum((df_kept$Oxygen_norm-mean(df_kept$Oxygen_norm))^2)
  fit_ok <- all(abs(pars[,"Std. Error"]/pars[,"Estimate"])<REL_SE_THRESHOLD,
                pars[,"Pr(>|t|)"] < PVAL_THRESHOLD,
                pseudo_R2 >= R2_THRESHOLD,
                diff(range(residuals(fit),na.rm=TRUE)) < MAX_RESID_RANGE,
                between(C_est,1e-12,1e-6), between(r_est,1e-4,0.1),
                AIC(lm(Oxygen_norm~1, data=df_kept)) - AIC(fit) >= AIC_IMPROVEMENT,
                mean(abs(residuals(fit)/df_kept$Oxygen_norm)) < MAPE_MAX)
  if(!fit_ok) next
  # conversions
  growth_C_fg_hr  <- r_est*60*init_biomass_fgC
  resp_molO2_hr   <- C_est*N0_fixed*mgL_to_mol_per_mL(O0)*60
  resp_C_fg_hr    <- resp_molO2_hr*12e15
  CUE <- growth_C_fg_hr/(growth_C_fg_hr+resp_C_fg_hr)
  # collect
  results <- add_row(results,
                     Taxon=Tax, Temperature=Temp, Replicate=Rep, N0=N0_fixed,
                     r_per_minute=r_est, r_per_hour=r_est*60, resp_rate=C_est,
                     O2_0=coef(fit)["O2_0"], AICc=AIC(fit),
                     lnO2_change_per_min=(log(last(df_kept$Oxygen_norm))-log(first(df_kept$Oxygen_norm)))/(last(df_kept$Time0)-first(df_kept$Time0)),
                     pseudo_R2=pseudo_R2, fit_ok=TRUE,
                     growth_C_fg_per_hr=growth_C_fg_hr, growth_C_fg_per_min=growth_C_fg_hr/60,
                     resp_C_fg_per_hr=resp_C_fg_hr, resp_C_fg_per_min=resp_C_fg_hr/60,
                     resp_to_growth_C=resp_C_fg_hr/growth_C_fg_hr, carbon_use_efficiency=CUE)
  # store diagnostic plot
  df_kept$Pred <- predict(fit, df_kept)
  plots_list[[paste(Tax,Temp,Rep,sep="_")]] <- ggplot(df_kept, aes(Time0,Oxygen_norm))+
    geom_point(size=2.5)+geom_line(aes(y=Pred),colour="red",linewidth=.8)+
    labs(title=paste(Tax,Temp,Rep,sep="_"),x=expression(Time~"(min)"),y=expression(Normalised~O[2]))+
    isme_theme()
}

## ───────────────────────── 10.  Save run-wise outputs ────────────────── ##
write_csv(results, "oxygen_model_results_good_only.csv")
if(length(plots_list)) {
  pdf("oxygen_dynamics_all_models.pdf",14,10); print(wrap_plots(plots_list)); dev.off()
  pdf("oxygen_dynamics_fullsize_per_page.pdf",8,6); walk(plots_list,print); dev.off()
}

##################################################################################
#                                   Section 11                                  #
#                     Temperature-response fits (3 panels)                      #
##################################################################################

# 11.1 Constant
k_B <- 8.617e-5  # eV K⁻¹

# 11.2 Helpers
safe_exp <- function(z) exp(pmin(700,z))
predict_safe <- function(fit,new) {
  out <- try(predict(fit,newdata=new),silent=TRUE)
  if(inherits(out,"try-error")) rep(NA_real_,nrow(new)) else out
}

# 11.3 Model functions
ln_SS_two <- function(T_C,lnR0,E,El,Tl,Eh,Th) {
  T<-T_C+273.15; TlK<-Tl+273.15; ThK<-Th+273.15
  lnR0 - E/(k_B*T) - log1p(safe_exp(El/k_B*(1/T-1/TlK))) - log1p(safe_exp(Eh/k_B*(1/ThK-1/T)))
}
ln_SS_one <- function(T_C,lnR0,E,Eh,Th) {
  T<-T_C+273.15; ThK<-Th+273.15
  lnR0 - E/(k_B*T) - log1p(safe_exp(Eh/k_B*(1/ThK-1/T)))
}
ln_boltz <- function(T_C,lnR0,E) {
  T<-T_C+273.15
  lnR0 - E/(k_B*T)
}
# Centre-shifted Arrhenius for CUE
arrh_mod_CUE2 <- function(T_C,A,E,B,T0=20) {
  T<-T_C+273.15; T0K<-T0+273.15
  A*exp(-E/k_B*(1/T-1/T0K))+B
}

# 11.4 Fit-and-plot helper
fit_T_plot <- function(df,yvar,ylab,show_x=TRUE) {
  df <- df %>% filter(is.finite(.data[[yvar]]), .data[[yvar]]>0)
  base <- ggplot(df,aes(Temperature,.data[[yvar]]))+geom_point(size=2.5)+
    labs(x=if(show_x) expression(Temperature~"(°C)") else NULL, y=ylab)+isme_theme()
  if(nrow(df)<4) return(base)
  T_min<-min(df$Temperature); T_max<-max(df$Temperature); T_opt<-df$Temperature[which.max(df[[yvar]])]
  # ---------------- CUE ----------------
  if(yvar=="carbon_use_efficiency") {
    start_vals <- c(A = min(df[[yvar]]) - max(df[[yvar]]), E = 2, B = min(df[[yvar]]))
    best_fit <- NULL; best_AIC <- Inf
    for(j in seq_len(4)) {
      sv <- if (j == 1) start_vals else start_vals * exp(rnorm(length(start_vals), 0, 0.2))
      ft <- try(nlsLM(
        carbon_use_efficiency ~ arrh_mod_CUE2(Temperature, A, E, B),
        data=df, start=sv,
        lower=c(A=-Inf, E=0.05, B=0), upper=c(A=0, E=3, B=1),
        control=nls.lm.control(maxiter=800)
      ), silent=TRUE)
      if(!inherits(ft, "try-error") && AIC(ft) < best_AIC) {
        best_fit <- ft; best_AIC <- AIC(ft)
      }
    }
    if (!is.null(best_fit)) {
      grid <- tibble(Temperature = seq(T_min, T_max, length.out = 300))
      grid$Pred <- predict_safe(best_fit, grid)
      grid <- grid %>% filter(is.finite(Pred))
      if (nrow(grid) > 1) {
        base <- base + geom_line(data = grid, aes(Temperature, Pred), linewidth = 1)
      }
    } else {
      message("⚠ carbon_use_efficiency: model failed to converge—plotting fallback curve")
      grid <- tibble(Temperature = seq(T_min, T_max, length.out = 300))
      start_A <- start_vals["A"]; start_E <- start_vals["E"]; start_B <- start_vals["B"]
      grid$Pred <- arrh_mod_CUE2(grid$Temperature, start_A, start_E, start_B)
      if (nrow(grid) > 1) {
        base <- base + geom_line(data = grid, aes(Temperature, Pred), linewidth = 1, linetype = "dashed")
      }
    }
    return(base)
  }
  # ------------- Growth & Respiration -------------
  # two-sided Sharpe–Schoolfield
  start_two <- list(lnR0=log(max(df[[yvar]])), E=0.6,
                    El=0.4, Tl=T_min+2, Eh=1.5, Th=T_opt+3)
  fit_two <- try(nlsLM(
    formula=paste0("log(",yvar,") ~ ln_SS_two(Temperature, lnR0, E, El, Tl, Eh, Th)"),
    data=df, start=start_two,
    lower=c(-Inf,0.1,0.1,0,0.1,0), upper=c(Inf,2.5,2.5,60,5,60),
    control=nls.lm.control(maxiter=600)
  ), silent=TRUE)
  # one-sided
  start_one <- list(lnR0=log(max(df[[yvar]])), E=0.6,
                    Eh=if(yvar=="growth_C_fg_per_hr") 3 else 1.5,
                    Th=if(yvar=="growth_C_fg_per_hr") T_opt+2 else T_opt+5)
  fit_one <- try(nlsLM(
    formula=paste0("log(",yvar,") ~ ln_SS_one(Temperature, lnR0, E, Eh, Th)"),
    data=df, start=start_one,
    lower=c(-Inf,0.1,0.5,0), upper=c(Inf,2.5,6,60),
    control=nls.lm.control(maxiter=500)
  ), silent=TRUE)
  # Arrhenius fallback
  fit_bol <- try(nlsLM(
    formula=paste0("log(",yvar,") ~ ln_boltz(Temperature, lnR0, E)"),
    data=df, start=list(lnR0=log(max(df[[yvar]])),E=0.6),
    lower=c(-Inf,0.1), upper=c(Inf,2.5),
    control=nls.lm.control(maxiter=400)
  ), silent=TRUE)
  getAIC <- function(f) if(inherits(f,"try-error")) Inf else AIC(f)
  fits <- list(two=fit_two, one=fit_one, bol=fit_bol)
  best_name <- names(which.min(map_dbl(fits,getAIC)))
  fit <- fits[[best_name]]
  if(getAIC(fit)==Inf) { message("⚠ ",yvar,": no model converged"); return(base) }
  grid <- tibble(Temperature=seq(T_min,T_max,length.out=300))
  grid$Pred <- safe_exp(predict_safe(fit,grid))
  grid <- grid %>% filter(is.finite(Pred))
  if(nrow(grid)>1) base <- base + geom_line(data=grid, aes(Temperature,Pred), linewidth=1)
  message("✅ ",yvar,": using ",best_name)
  base
}

# 11.5 Filter acceptable fits
auto <- results %>% filter(fit_ok, growth_C_fg_per_hr>0, resp_C_fg_per_hr>0)

# 11.6 Build panels
p_growth <- fit_T_plot(auto, "growth_C_fg_per_hr", expression(Growth~(fg~C~h^{-1})), TRUE)
p_resp   <- fit_T_plot(auto, "resp_C_fg_per_hr",   expression(Respiration~(fg~C~h^{-1})), TRUE)
p_cue    <- fit_T_plot(auto, "carbon_use_efficiency", "Carbon Use Efficiency", TRUE)

# 11.7 Extract parameters & T_opt
extract_params <- function(fit, name, var, Topt) {
  if(is.null(fit) || inherits(fit,"try-error")) return(NULL)
  cf <- coef(fit)
  tibble(
    Trait=var, Model=name, AIC=AIC(fit), lnR0=cf["lnR0"] %||% NA,
    E=cf["E"] %||% NA, El=cf["El"] %||% NA, Tl=cf["Tl"] %||% NA,
    Eh=cf["Eh"] %||% NA, Th=cf["Th"] %||% NA,
    A=cf["A"] %||% NA, B=cf["B"] %||% NA, T_opt=Topt
  )
}
params_list <- list()
for(var in c("growth_C_fg_per_hr","resp_C_fg_per_hr","carbon_use_efficiency")) {
  df <- auto %>% filter(is.finite(.data[[var]]), .data[[var]]>0)
  if(nrow(df)<4) next
  Tmin <- min(df$Temperature); Tmax <- max(df$Temperature)
  Topt_data <- df$Temperature[which.max(df[[var]])]
  grid1000 <- tibble(Temperature=seq(Tmin,Tmax,length.out=1000))
  if(var=="carbon_use_efficiency") {
    fit <- try(nlsLM(
      carbon_use_efficiency ~ arrh_mod_CUE2(Temperature,A,E,B),
      data=df, start=list(A=min(df[[var]])-max(df[[var]]), E=0.6, B=min(df[[var]])),
      lower=c(A=-Inf,E=0.01,B=0), upper=c(A=0,E=10,B=1),
      control=nls.lm.control(maxiter=500)
    ), silent=TRUE)
    pred <- predict_safe(fit, grid1000)
    Topt_fit <- grid1000$Temperature[which.max(pred)]
    params_list[[var]] <- extract_params(fit, "Arrhenius+CUE", var, Topt_fit)
  } else {
    fits <- list(
      two = try(nlsLM(
        formula=paste0("log(",var,") ~ ln_SS_two(Temperature, lnR0, E, El, Tl, Eh, Th)"),
        data=df, start=list(lnR0=log(max(df[[var]])),E=0.6,El=0.4,Tl=Tmin+2,Eh=1.5,Th=Topt_data+3),
        lower=c(-Inf,0.1,0.1,0,0.1,0), upper=c(Inf,2.5,2.5,60,5,60), control=nls.lm.control(maxiter=500)
      ), silent=TRUE),
      one = try(nlsLM(
        formula=paste0("log(",var,") ~ ln_SS_one(Temperature, lnR0, E, Eh, Th)"),
        data=df, start=list(lnR0=log(max(df[[var]])),E=0.6,Eh=if(var=="growth_C_fg_per_hr")3 else 1.5,Th=if(var=="growth_C_fg_per_hr")Topt_data+2 else Topt_data+5),
        lower=c(-Inf,0.1,0.5,0), upper=c(Inf,2.5,6,60), control=nls.lm.control(maxiter=500)
      ), silent=TRUE),
      bol = try(nlsLM(
        formula=paste0("log(",var,") ~ ln_boltz(Temperature, lnR0, E)"),
        data=df, start=list(lnR0=log(max(df[[var]])),E=0.6), lower=c(-Inf,0.1), upper=c(Inf,2.5), control=nls.lm.control(maxiter=400)
      ), silent=TRUE)
    )
    AICs <- map_dbl(fits, ~ if(inherits(.x,"try-error")) Inf else AIC(.x))
    best <- names(which.min(AICs))
    fit <- fits[[best]]
    pred <- exp(predict_safe(fit, grid1000))
    Topt_fit <- grid1000$Temperature[which.max(pred)]
    params_list[[var]] <- extract_params(fit, best, var, Topt_fit)
  }
}
all_params <- bind_rows(params_list)
write_csv(all_params, "thermal_fit_parameters_with_Topt.csv")

# 11.8 Combine & export plots
combo_plot <- wrap_plots(p_growth, p_resp, p_cue, ncol=3) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels="A", tag_prefix="(", tag_suffix=")") &
  theme(plot.margin=margin(5,10,5,10))

ggsave("SharpeSchoolfield_Temperature_Fits.pdf", combo_plot,
       width=12, height=4, dpi=600)

################################################################################
#                                   End of Script                             #
################################################################################
