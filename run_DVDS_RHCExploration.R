
rerun_sensitivity = F 
quant_boots  = 500 

source("DVDS.R")

library(data.table)
library(ggplot2)
library(stargazer)
library(lfe)

if ("rhc.csv" %in% list.files()) {
  rhc = fread("rhc.csv")
} else {
  #Downloaded from Github --- should match Hirano and Imbens, "Estimation of Causal Effects using Propensity Score Weighting: An Application to Data on Right Heart Catheterization"
  library(BWSnooping)
  
  data(rhc)
  
  rhc = as.data.table(rhc)
  fwrite(rhc, "rhc.csv")
}


#Only 52 covariates in the dataset
ncol(rhc)

#"Survival" is actually non-survival (matching Connors et al., Table 2)
rhc[, survival := !survival]

#3551 controls (RHC == FALSE), 2184 treated (rhc == TRUE)
print(rhc[, .N, by = rhc])
print(dcast(rhc[, .N, by = .(rhc, survival)], survival ~ rhc, value.var = "N"))


#======================================================================
# Prep Data
#======================================================================

#Demean covariates (based on https://stackoverflow.com/a/30017723)
control_vars = setdiff(names(rhc), c("rhc", "survival", "propwt"))
demean_vars = c()
for (col in control_vars) {
  these_vals = rhc[, get(col)]
  if (class(these_vals) == "character") {
    all_vals = unique(these_vals)
    for (i in 1:(length(all_vals)-1)) {
      rhc[, paste0(col, i) := (these_vals == all_vals[i]) - mean((these_vals == all_vals[i]))]
      demean_vars = c(demean_vars, paste0(col, i))
    }
  } else {
    rhc[, c(col) := these_vals]
    demean_vars = c(demean_vars, col)
  }
}

#Variables will be demeaned for linear regression, not for logit (since it's weighted demeaning)
run_formula_lm    = "survival ~ rhc"
run_formula_logit = "rhc ~ 1"
for (el in demean_vars) {
  run_formula_lm    = paste0(run_formula_lm,    " + ", el, " + ", el, ":rhc")
  run_formula_logit = paste0(run_formula_logit, " + ", el)
}



#======================================================================
# Nominal Propensity Checks
#======================================================================


#Demean (weighted) and standardize (unweighted) 
ehat_x = predict(glm(as.formula(run_formula_logit), data = rhc, family = "binomial"), type = "response")
rhc[, propwt := (rhc/ehat_x + (1-rhc)/(1-ehat_x))]
rhc[, paste(demean_vars) := lapply(.SD, function(x) (as.numeric(x)-sum(as.numeric(x)*propwt)/sum(propwt)) / sd(as.numeric(x))), .SDcols=demean_vars]


#Attempt at Figure 2 (close but not quite right)
ggplot(data.frame(Z = rhc$rhc, ehat = ehat_x), aes(x = ehat)) + 
  geom_histogram(binwidth=1/30, center = 1/60) + facet_wrap(~ Z, ncol = 1, scales = "free_y")

#Plausible propensities? E[Z | ehat_x] seems to be about ehat_x 
ggplot(data.frame(Z = as.numeric(rhc$rhc), ehat = ehat_x), aes(x = ehat_x, y = Z)) + 
  geom_smooth() + geom_abline(intercept = 0, slope = 1, lty = 2)

#======================================================================
# Propensity Calibration 
#======================================================================

glm_estims = coef(summary(glm(as.formula(run_formula_logit), data = rhc, family = "binomial")))
#Drop categorical variables
glm_estims = glm_estims[!row.names(glm_estims) == "(Intercept)", 1:4]
for (this_col in setdiff(names(rhc), names(demean_vars))) {
  if (class(rhc[, get(this_col)]) == "character") {
    glm_estims = glm_estims[!row.names(glm_estims) %like% this_col, 1:4]
  }
}

#Can do benchmarking here 
glm_estims = glm_estims[!row.names(glm_estims) %like% "cat1" & !row.names(glm_estims) %like% "cat2" & !row.names(glm_estims) %like% "insurance" & !row.names(glm_estims) == "(Intercept)", 1:4]
glm_estims = glm_estims[order(-abs(glm_estims[,1])), 1:4]

#======================================================================
# Run Hirano and Imbens WLS (comparable to results in their Table 3)
#======================================================================

propwt = rhc$propwt

#Latter for SEs, former for prediction
lm_results   =   lm(as.formula(run_formula_lm), data = rhc, weights = propwt)
felm_results = felm(as.formula(run_formula_lm), data = rhc, weights = propwt)
stargazer(felm_results, se = list(coef(summary(felm_results, robust = T))[, 2]), type = "text", keep = 1)

muhat_1 = predict(lm_results, copy(rhc)[, rhc := 1])
muhat_0 = predict(lm_results, copy(rhc)[, rhc := 0])




#======================================================================
# Run AIPW 
#======================================================================

rhc[, AIPW_ATEi := muhat_1 - muhat_0 + rhc * (survival - muhat_1)/ehat_x - (1-rhc) * (survival - muhat_0)/(1-ehat_x)]

#AIPW ATE estimate: 0.063, similar to weighted regression approach --- reject zero 
rhc[, mean(AIPW_ATEi)]

coef(summary(felm_results))[2,1]

#IPW ATE estimate: 0.038, close to threshold of AIPW CI (a potential issue with IPW approaches)
rhc[, mean(rhc * survival /ehat_x - (1-rhc) * survival /(1-ehat_x))]


#==================================
# Sensitivity Analysis 
#==================================

lin_covs  = paste(c(demean_vars), collapse = " + ")

#What specs to run? 
spec_list = list()

for (this_component in c("linear", "boost")) {
  for (this_formula in c("zsb", "qb", "dvds")) {
    
    
    add_function = paste0(this_formula, "(", paste(rep(" ", 14 - nchar(this_formula)), collapse = ""), "Lambdas, data, lin_covs, ")
    
    #Propensity model: Boost or linear
    add_function = paste0(add_function, "'Z', 'Y', ", ifelse(this_component == "boost", "boostregn, boostregn_option_bin, ", "  linregn,   linregn_option_bin, "))
    
    
    if (this_formula != "zsb") {
      #Quantiles: linear always 
      add_function = paste0(add_function, " binaryExtrapolationquant, ")
      add_function = paste0(add_function, ifelse(this_component == "boost", "binaryExtrapolation_option_boost, ", "binaryExtrapolation_option_linear, "))
      
      #Kappa-hat estimation 
      if (this_formula == "dvds") {
        add_function = paste0(add_function, " binaryExtrapolationKappa, ")
        add_function = paste0(add_function, ifelse(this_component == "boost", "binaryExtrapolation_option_boost, ", "binaryExtrapolation_option_linear, "))
        add_function = paste0(add_function, "semiadaptive=T, ")
      }
    }
    
    add_function = paste0(
      add_function, 
      "boot_infer = T, boot_settings = list(reuse_boot = T, refit_propensities = ",
      ifelse(this_component == "linear", "T", "F"), 
      ", B = ", quant_boots, ", return_seed = TRUE), "
    )
    
    
    spec_list[[length(spec_list) + 1]] = paste0(
      "results.", this_formula, ".", this_component, paste(rep(" ", 17 - nchar(this_component) + 4 - nchar(this_formula)), collapse = ""), " = ", 
      add_function, "trim.type='clip', normalize = 'F')"
    )
    
    spec_list[[length(spec_list) + 1]] = paste0(
      "results.", this_formula, ".", this_component, ".norm0", paste(rep(" ", 12 - nchar(this_component) + 4 - nchar(this_formula)), collapse = ""), " = ", 
      add_function, "trim.type='clip', normalize = '0')"
    )
    
    spec_list[[length(spec_list) + 1]] = paste0(
      "results.", this_formula, ".", this_component, ".norm1", paste(rep(" ", 12 - nchar(this_component) + 4 - nchar(this_formula)), collapse = ""), " = ", 
      add_function, "trim.type='clip', normalize = '1')"
    )
    
    
  }
}

#======================================================================
# Run Sensitivity Analysis  
#======================================================================

Lambdas = sort(unique(c(seq(1, 1.5, length.out = 11), seq(1, 2, length.out = 11))))

if (rerun_sensitivity) {
  
  data = copy(rhc)
  setnames(data, c("rhc", "survival"), c("Z", "Y"))
  
  for (el in spec_list) {
    print("===================")
    cat(paste0("Running ", el, "\n"))
    print("===================")
    
    set.seed(16)
    eval(parse(text=el))
  }
  
  
  all_results_combined = data.frame()
  for (these_results in ls()[grepl("results\\.", ls())]) {
    to_add = get(these_results) %>% mutate(method = gsub("results.", "", these_results, fixed = T))
    all_results_combined = rbind(all_results_combined, to_add)
  }
  
  write.csv(all_results_combined, "RHCExploration_sensitivity.csv")
  
} else {
  all_results_combined = fread("RHCExploration_sensitivity.csv")
}

fwrite(coef(summary(felm_results, robust = T)), "RHCE_OLSResults.csv")

#===============================
# Plot Sensitivity Analysis
#===============================

nominal_estimate = coef(summary(felm_results))[2,1] 
nominal_CIUB     = coef(summary(felm_results))[2,1] + coef(summary(felm_results, robust = T))[2,2] * qnorm(0.95)
nominal_CILB     = coef(summary(felm_results))[2,1] + coef(summary(felm_results, robust = T))[2,2] * qnorm(0.05)

all_results_combined = as.data.table(all_results_combined)
if ("V1" %in% names(all_results_combined)) {
  all_results_combined[, V1 := NULL]
}

#90% CI
all_results_combined[, ':='(
  CI_IIF   = ifelse(upper, estimate + qnorm(0.95) * sterr_iif, estimate + qnorm(0.05) * sterr_iif),
  CI_Boot  = ifelse(upper, estimate + qnorm(0.95) * sterr_boot, estimate + qnorm(0.05) * sterr_boot),
  CI_PBoot = quantiles90
  )]

sensitivity_results = dcast(all_results_combined, Lambda + estimand + method ~ upper, value.var = c("estimate", "CI_IIF", "CI_Boot", "CI_PBoot"))
setnames(sensitivity_results, c("estimate_FALSE", "estimate_TRUE", "CI_IIF_FALSE", "CI_IIF_TRUE", "CI_Boot_FALSE", "CI_Boot_TRUE", "CI_PBoot_FALSE", "CI_PBoot_TRUE"), c("Low", "Hi", "q5_IIF", "q95_IIF", "q5_Boot", "q95_Boot", "q5_PBoot", "q95_PBoot"))
sensitivity_results[, ':='(
  q5_Joint  = ifelse(method %like% "dvds",  q5_IIF,  q5_PBoot),
  q95_Joint = ifelse(method %like% "dvds", q95_IIF, q95_PBoot)
)]
sensitivity_results[, ':='(
  Estimator = strsplit(method, "\\.")[[1]][1],
  Nuisance  = strsplit(method, "\\.")[[1]][2],
  Normlztn  = strsplit(method, "\\.")[[1]][3]
), by = method]
sensitivity_results[, Normlztn := ifelse(is.na(Normlztn), "None", gsub("norm", "", Normlztn))]
sensitivity_results[, MethodPrint := paste0(
  toupper(Estimator), " (",
  ifelse(Nuisance == "boost", "Boost Non-Quantile", ifelse(Nuisance == "linear", "Linear Nuisances", ifelse(Nuisance == "boostokprop", "Boost Kappa", "Nuisance Unclear"))),
  ")"
)]
setkeyv(sensitivity_results, c("Nuisance", "Estimator"))
sensitivity_results[, MethodPrint := factor(MethodPrint, levels = unique(MethodPrint))]


plotSensitivity = function(lambda_vals = 1 + (0:8)/20, this_estimand = "ATE", this_estimators = c("dvds", "qb", "zsb"), this_nuisances = c("boost", "linear"), this_norm = "None",
                           bar_width = 0.02, color_column = "MethodPrint", y_add = NULL, this_infer) {
  
  to_plot = sensitivity_results[round(Lambda, 3) %in% round(lambda_vals, 3) & 
                                  estimand == this_estimand & Estimator %in% this_estimators & Nuisance %in% this_nuisances & Normlztn == this_norm]
  to_plot = to_plot[, .(
    Lambda, ColorValue = get(color_column), Low, Hi,
    q5 = get(paste0("q5_", this_infer)), q95 = get(paste0("q95_", this_infer))
  )]
  
  return(
    ggplot(to_plot, aes(x = Lambda, color = ColorValue)) +
      geom_errorbar(aes(ymin = Low, ymax = Hi, width = bar_width), position = position_dodge()) +
      geom_errorbar(aes(ymin = q5, ymax = q95, width = bar_width), position = position_dodge(), 
                    alpha = 0.3, linetype = 1) +
      theme_bw() + theme(panel.grid.major = element_blank()) +
      labs(x = "Lambda", y = paste0("Bounds", ifelse(is.null(y_add), "", y_add)), color = "Method", x = expression(paste(Lambda))) +
      geom_hline(yintercept = 0, alpha = 0.9, lty = 1) +
      #geom_hline(yintercept = nominal_estimate, alpha = 0.7, lty = 2) + 
      #geom_ribbon(aes(ymin = nominal_CILB, ymax = nominal_CIUB), alpha = 0.05, color = NA) +
      scale_x_continuous(breaks = seq(0.1 * floor(  min(to_plot$Lambda*10)), 0.1 * ceiling(max(to_plot$Lambda*10)), by = 0.1))
  )
}

pdf(file = "Results/Figs/RHC_InitialSensitivity.pdf", width = 6, height = 4)
print(plotSensitivity(this_infer = "Joint"))
dev.off()

pdf(file = "Results/Figs/RHC_InitialSensitivity_boostOnly.pdf", width = 4, height = 3)
print(plotSensitivity(this_nuisances = "boost", color_column = "Estimator", y_add = " (Boost Nuisance)", this_infer = "Joint"))
dev.off()

pdf(file = "Results/Figs/RHC_InitialSensitivity_linearOnly.pdf", width = 4, height = 3)
print(plotSensitivity(this_nuisances = "linear", color_column = "Estimator", y_add = " (Logit Nuisance)", this_infer = "Joint"))
dev.off()



#============================
# Facet Plot
#============================


to_plot = sensitivity_results[round(Lambda, 3) %in% round(1 + (0:8)/20, 3) &  estimand == "ATE" & Estimator %in% c("dvds", "qb", "zsb") & Nuisance %in% c("boost", "linear") & Normlztn == "None", .(
  Lambda, Nuisance, ColorValue = get("Estimator"), Low, Hi,
  q5 = q5_Joint, q95 = q95_Joint
)]
ggplot(to_plot, aes(x = Lambda, color = ColorValue)) + 
  facet_wrap(~ Nuisance) + theme_bw() + theme(panel.grid.major = element_blank()) +
  labs(x = "Lambda", y = "Bounds", color = "Method", x = expression(paste(Lambda))) +
  geom_hline(yintercept = 0, alpha = 0.9, lty = 1) +
  scale_x_continuous(breaks = seq(0.1 * floor(  min(to_plot$Lambda*10)), 0.1 * ceiling(max(to_plot$Lambda*10)), by = 0.1)) +
  geom_line(aes(y = Low), alpha = 1.0) + geom_line(aes(y = Hi),  alpha = 1.0) +
  geom_line(aes(y = q5),  alpha = 0.3) + geom_line(aes(y = q95), alpha = 0.3)

rm(to_plot)

#======================================================================
# Case Matching Estimates (Connors)
#======================================================================

library(MatchIt)

rhc_match = copy(rhc)[, .(survival, rhc, ehat_x = ehat_x, cat1, cat2)]

matched_vals = matchit(rhc ~ ehat_x, rhc_match, exact = c("cat1", "cat2"), replace = F)

rhc_match[, match_id := matched_vals$subclass]

pairs_diffs = rhc_match[!is.na(match_id), .(survival_diff = sum(survival * (2 * rhc - 1))), by = match_id]
print(pairs_diffs[, .N, by = survival_diff][, .(survival_diff, Pairs = N, Rate = round(100 * N / sum(N), 1))])


#=====================================================
# Summarize calibration_dt (if created in debugging) 
#=====================================================

if (exists("calibration_df")) {
  calibration_df$ZWrite = ifelse(calibration_df$Z, "Z=1", "Z=0")
  
  pdf(file = "Results/Figs/boost_propensity_hist.pdf", width = 6, height = 4)
  print(
    ggplot(calibration_df, aes(x = ehat_x)) + geom_histogram(bins = 30) + 
      facet_wrap(~ZWrite, ncol = 1) + labs(x = "ehat(X_i), Boost"))
  dev.off()
  
  pdf(file = "Results/Figs/boost_propensity_calibration.pdf", width = 6, height = 4)
  print(
    ggplot(calibration_df, aes(x = ehat_x, y = as.numeric(Z))) + 
      geom_smooth() + geom_abline(slope = 1, lty = 2) +
      labs(x = "ehat(X_i), Boost", y = paste0("Z_i (mean Z_i/ehat(X_i) = ", round(mean(calibration_df$Z / calibration_df$ehat_x), 3), ")"))
  )
  dev.off()
  
  pdf(file = "Results/Figs/boost_muhat1_calibration.pdf", width = 6, height = 4)
  print(
    ggplot(calibration_df[calibration_df$Z==1,], aes(x = muhat_1, y = as.numeric(Y))) + 
      geom_smooth() + geom_abline(slope = 1, lty = 2) + 
      labs(x = "\\hat{\\mu}(X_i, 1)", y = "Y(1)")
  )
  dev.off()
  
  pdf(file = "Results/Figs/boost_muhat0_calibration.pdf", width = 6, height = 4)
  print(
    ggplot(calibration_df[calibration_df$Z==0,], aes(x = muhat_0, y = as.numeric(Y))) + 
      geom_smooth() + geom_abline(slope = 1, lty = 2) + 
      labs(x = "\\hat{\\mu}(X_i, 0)", y = "Y(0)")
  )
  dev.off()
  
}
