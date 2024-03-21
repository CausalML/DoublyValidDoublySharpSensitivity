rerun_firststage  = T
rerun_sensitivity = T 
K                 = 5
B                 = 10000
lambda_values     = unique(sort(c(seq(1, 1.5, by = 0.001), seq(1, 3, by = 0.01))))

source("DVDS.R")
source("simplified_dvds.R")


library(data.table)
library(ggplot2)
library(progress)
library(SuperLearner)
library(gridExtra)
library(e1071)

if ("rhc.csv" %in% list.files("Real_data")) {
  rhc = read.csv("Real_data/rhc.csv")
} else {
  # Downloaded from Github --- devtools::install_github("kolesarm/BWSnooping")
  # Should match Hirano and Imbens, "Estimation of Causal Effects using Propensity Score Weighting: An Application to Data on Right Heart Catheterization"
  library(BWSnooping)
  data(rhc)
  write.csv(rhc, "Real_data/rhc.csv")
}

rhc = as.data.table(rhc)
ncol(rhc) #Only 52 covariates in the dataset

#"Survival" is actually non-survival (matching Connors et al., Table 2)
rhc[, survival := !survival]

#3551 controls (RHC == FALSE), 2184 treated (rhc == TRUE)
print(rhc[, .N, by = rhc])
print(dcast(rhc[, .N, by = .(rhc, survival)], survival ~ rhc, value.var = "N"))


#======================================================================
# Prep Data
set.seed(16)

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
  if (class(rhc[, get(el)]) == "logical") {
    rhc[, c(el)] = as.numeric(rhc[, get(el)])
  }
}

#======================================================================
# Implement Super Learners
data = copy(rhc)
setnames(data, c("rhc", "survival"), c("Z", "Y"))
if (rerun_firststage) {
  cvgroup = make.cvgroup.balanced(data, K, 'Z')
  data$K  = cvgroup 
  #Generate component predictions for SuperLearner 
  preds = list()
  pb <- progress_bar$new(total = 3 * K, format = "  estimating [:bar] :percent eta: :eta")
  for (this_out in c('Z', 'Y1', "Y0")) {
    sub_data = copy(data)
    if (this_out != 'Z') {
      include_row = (sub_data$Z == as.numeric(gsub("Y", "", this_out)))
    } else {
      include_row = rep(TRUE, nrow(sub_data))
    }
    preds[[this_out]] = data.frame('boost' = rep(0, nrow(sub_data)), 'svm' = 0, 'forest' = 0, 'lin' = 0)
    for (k in unique(sub_data$K)) {
      #Get values for dvds argument 
      trainmask = include_row & sub_data$K != k
      testmask  = sub_data$K == k
      form_x    = paste(demean_vars, collapse = " + ")
      form_resp = substr(this_out, 1, 1)
      preds[[this_out]][, "boost"][ testmask] =  boostregn(sub_data, trainmask, testmask, form_x, form_resp,  boostregn_option_bin)
      preds[[this_out]][, "svm"][   testmask] =    svmregn(sub_data, trainmask, testmask, form_x, form_resp,    svmregn_option_bin)
      preds[[this_out]][, "forest"][testmask] = forestregn(sub_data, trainmask, testmask, form_x, form_resp, forestregn_option_bin)
      preds[[this_out]][, "lin"][   testmask] =    linregn(sub_data, trainmask, testmask, form_x, form_resp,    linregn_option_bin)
      pb$tick()
    }
  }
  all_preds = data.frame()
  for (this_name in names(preds)) {
    to_add = preds[[this_name]]
    to_add$Outcome = this_name
    all_preds = rbind(all_preds, to_add)
  }
  write.csv(all_preds, "Real_data/Results/Tabs/RHCExploration_predictions.csv", row.names = FALSE)
}
all_preds = read.csv("Real_data/Results/Tabs/RHCExploration_predictions.csv")


#======================================================================
# Run SuperLearner to maximize log probability 

sl_preds = data.frame()
for (this_out in unique(all_preds$Outcome)) {
  meta_preds = copy(all_preds[all_preds$Outcome == this_out, ])
  #Change to just SL components 
  meta_preds[, c("Outcome")] = NULL
  meta_preds$Zero = 0
  meta_preds$One  = 1
  train_preds = copy(meta_preds)
  train_outs  = as.numeric(data[, get(substr(this_out, 1, 1))])
  if (this_out != 'Z') {
    train_preds = train_preds[data$Z == as.numeric(gsub("Y", "", this_out)),]
    train_outs  = train_outs[ data$Z == as.numeric(gsub("Y", "", this_out))]
  }
  #Maximize logistic error 
  train_preds = as.matrix(train_preds)
  weightScore = function(w) {
    w = w / sum(w)
    preds = train_preds %*% w
    return(-1 * sum(train_outs * log(preds) + (1-train_outs) * log(1-preds)))
  }
  optimal_weights = optim(rep(1/ncol(train_preds), ncol(train_preds)), weightScore, method = "L-BFGS-B", lower = 0)$par
  to_print = optimal_weights 
  names(to_print) <- names(meta_preds)
  print(paste0("The optimal weights for ", this_out, " are:"))
  print(round(100 * to_print / sum(to_print), 2))
  optimal_weights = optimal_weights / sum(optimal_weights) 
  sl_preds = rbind(sl_preds, 
                   data.frame(
                     Outcome = rep(this_out, nrow(data)), 
                     i = 1:nrow(data), 
                     TrueValues = ifelse(rep(this_out, nrow(data)) == "Z", data$Z, ifelse(data$Z == as.numeric(gsub("Y", "", this_out)), data$Y, as.numeric(NA)) ),
                     SuperLearner = as.matrix(meta_preds) %*% optimal_weights)
  )
}
all_preds$i = 1:nrow(data)
all_preds = merge(sl_preds, all_preds, by = c("Outcome", "i"), all.x = TRUE, all.y = TRUE)
all_preds = all_preds[order(all_preds$Outcome, all_preds$i),]

ehat   = all_preds[all_preds$Outcome == "Z",  ("SuperLearner")]
mu1hat = all_preds[all_preds$Outcome == "Y1", ("SuperLearner")]
mu0hat = all_preds[all_preds$Outcome == "Y0", ("SuperLearner")]

#======================================================================
# Sensitivity Analysis



if (rerun_sensitivity) {
  sensitivity_results = data.frame()
  pb <- progress_bar$new(total = length(lambda_values), format = "  sensitivity [:bar] :percent eta: :eta")
  for (lambda_val in lambda_values) {
    
    
    #Run estimates 
    to_add = rbind(
      extrema.aipw(data$Z, data$Y, ehat, mu0hat, mu1hat, Lambda = lambda_val),
      dvds.binary( data$Y, data$Z, mu0hat, mu1hat, ehat, Lambda = lambda_val)
    )
    #Style and add 
    to_add$Method = c("ZSB", "DVDS")
    to_add$Lambda = lambda_val
    sensitivity_results = rbind(sensitivity_results, to_add)
    pb$tick()
  }
  write.csv(sensitivity_results, "Real_data/Results/Tabs/sensitivity_results_binary.csv", row.names = FALSE)
  
  #Find crossing points 
  dvds.binary.wrapper = function(lambda, this_df, e, alpha=0.05, ret_upper = T) {
    if (ret_upper) {
      dvds.binary(this_df$Y, this_df$Z, this_df$mu0hat, this_df$mu1hat, this_df$ehat, lambda, alpha)$upper
    } else {
      dvds.binary(this_df$Y, this_df$Z, this_df$mu0hat, this_df$mu1hat, this_df$ehat, lambda, alpha)$lower
    }
  }
  min_negative = 2 * min(sensitivity_results[sensitivity_results$Method == "DVDS" & sensitivity_results$upper.CI > 0, "Lambda"])
  full_df      = data.frame(ind = 1:nrow(data), Y = data$Y, Z = data$Z, mu0hat, mu1hat, ehat)
  pb           = progress_bar$new(total = B, format = "  bootstrap [:bar] :percent eta: :eta")
  min_lambdas  = c()
  for (b in 1:B) {
    set.seed(1616 + b)
    min_lambdas = c(min_lambdas, uniroot(dvds.binary.wrapper, lower = 1, upper = min_negative, this_df = full_df[sample(full_df$ind, nrow(full_df), replace = T), ])$root)
    pb$tick()
  }
  min_lambdas = data.frame(seed_num = 1616 + 1:B, min_lambda_val = min_lambdas)
  write.csv(min_lambdas, "Real_data/Results/Tabs/sensitivity_results_bootstrap_crossing.csv", row.names = FALSE)
} else {
  sensitivity_results = read.csv("Real_data/Results/Tabs/sensitivity_results_binary.csv")
  min_lambdas         = read.csv("Real_data/Results/Tabs/sensitivity_results_bootstrap_crossing.csv")
}


plotSensitivity = function(lambda_vals = unique(sensitivity_results$Lambda), this_estimand = "ATE", this_estimators = c("DVDS", "ZSB"), 
                           bar_width = 0.02, color_column = "Method", y_add = NULL, this_infer = "Joint", make_as_line = F) {
  
  to_plot = sensitivity_results[round(Lambda, 3) %in% round(lambda_vals, 3) & Method %in% this_estimators]
  to_plot = to_plot[, .(
    Lambda, ColorValue = as.character(get(color_column)), lower, upper,
    lower.CI, upper.CI
  )]
  
  if (make_as_line) {
    return(
      ggplot(to_plot, aes(x = Lambda, color = ColorValue, fill = ColorValue)) +
        geom_line( aes(y = lower)) + geom_line( aes(y = upper)) + 
        geom_point(aes(y = lower)) + geom_point(aes(y = upper)) +
        geom_ribbon(aes(ymin = lower.CI, ymax = lower), alpha = 0.3, color = NA) +
        geom_ribbon(aes(ymin = upper, ymax = upper.CI), alpha = 0.3, color = NA) +
        theme_bw() + theme(panel.grid.major = element_blank()) +
        labs(x = expression(Lambda), y = paste0("Bounds", ifelse(is.null(y_add), "", y_add)), color = "Method", fill = "Method", x = expression(paste(Lambda))) +
        geom_hline(yintercept = 0, alpha = 0.9, lty = 1) + 
        scale_color_manual(values = c("DVDS" = "#F8766D", "ZSB" = "#619CFF", "odds" = "black")) +
        scale_fill_manual(values = c("DVDS" = "#F8766D", "ZSB" = "#619CFF", "odds" = "black")) +
        scale_x_continuous(breaks = seq(0.1 * floor(  min(to_plot$Lambda*10)), 0.1 * ceiling(max(to_plot$Lambda*10)), by = 0.1))
    )
  } else {
    return(
      ggplot(to_plot, aes(x = Lambda, color = ColorValue)) +
        geom_errorbar(aes(ymin = lower, ymax = upper, width = bar_width), position = position_dodge()) +
        geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI, width = bar_width), position = position_dodge(), 
                      alpha = 0.3, linetype = 1) +
        theme_bw() + theme(panel.grid.major = element_blank()) +
        labs(x = expression(Lambda), y = paste0("Bounds", ifelse(is.null(y_add), "", y_add)), color = "Method", x = expression(paste(Lambda))) +
        geom_hline(yintercept = 0, alpha = 0.9, lty = 1) + 
        scale_color_manual(values = c("DVDS" = "#F8766D", "ZSB" = "#619CFF", "odds" = "black")) +
        #geom_hline(yintercept = nominal_estimate, alpha = 0.7, lty = 2) + 
        #geom_ribbon(aes(ymin = nominal_CILB, ymax = nominal_CIUB), alpha = 0.05, color = NA) +
        scale_x_continuous(breaks = seq(0.1 * floor(  min(to_plot$Lambda*10)), 0.1 * ceiling(max(to_plot$Lambda*10)), by = 0.1))
    )
  }
}

getORs = function(col_name, dt = pred_drops) {
  dt[, get(col_name) / (1-get(col_name)) / (None / (1-None))]
}

#Add calibration points to plot 
addCalibration = function(to_plot, curr_max_lambda = 1.5, add_cols = c("ad_neuro", "surv2md1", "pafi1")) {
  for (i in 1:length(add_cols)) {
    this_col = add_cols[i]
    col_ORs = getORs(this_col)
    col_ORs = pmax(col_ORs, 1/col_ORs)
    x_low = min(quantile(col_ORs, 0.25), curr_max_lambda)
    x_hi  = min(quantile(col_ORs, 0.75), curr_max_lambda)
    y_level = -0.15 - 0.05 * i
    arrow_ends = ifelse(x_hi == curr_max_lambda, "first", "both")
    label_name = ifelse(this_col == "sex", "Sex", ifelse(this_col == "dnr1", "DNR Status (Day 1)", ifelse(this_col == "surv2md1", "Estim. surv. prob.", ifelse(this_col == "pafi1", "PaO2/FI02 ratio", ifelse(
      this_col == "ad_neuro", "Neuro. diag.", paste0(this_col, " (Formatting TODO)"))))))
    #Add boxes
    to_plot = to_plot + 
      annotate("rect", xmin = x_low, xmax = x_hi, ymin = y_level - 0.01, ymax = y_level + 0.01, alpha = 0, color = "black") +
      annotate("text", x = (x_low + x_hi) / 2, y = y_level + 0.025, label = label_name, size = 3)
    #Add whiskers 
    if (quantile(col_ORs, 0.1) < min(quantile(col_ORs, 0.25), curr_max_lambda)) {
      to_plot = to_plot +  annotate("segment", x = quantile(col_ORs, 0.1), xend = x_low, y = y_level, yend = y_level, arrow = arrow(ends = "first", angle = 90, length = unit(.2, "cm")))
    }
    if (quantile(col_ORs, 0.75) < min(quantile(col_ORs, 0.9), curr_max_lambda)) {
      #Hiding whisker in the box 
      arrow_ends = ifelse(quantile(col_ORs, 0.9) < curr_max_lambda, "first", "last")
      to_plot = to_plot +  annotate("segment", x = min(quantile(col_ORs, 0.9), curr_max_lambda), xend = x_hi, y = y_level, yend = y_level, arrow = arrow(ends = arrow_ends, angle = 90, length = unit(.2, "cm")))
    }
  }
  return(to_plot)
}


#======================================================================
# Get calibration points for plot 

test_dt = as.data.table(copy(data))

#Don't drop levels of factors 
test_dt[, c(paste0("cat1", 1:8), paste0("cat2", 1:6), paste0("cancer", 1:2)) := NULL]
control_cols = setdiff(names(test_dt), c("Z", "Y"))

#Baseline results 
glm_results = glm(paste0("Z ~ ", paste(control_cols, collapse = " + ")), data = test_dt, family = "binomial")

#Results dropping variables 
pred_drops = data.table(obs_id = 1:nrow(test_dt), Z = test_dt$Z, pred_SL = all_preds[all_preds$Outcome == "Z",]$SuperLearner, None = predict(glm_results, type = "response"))
for (this_name in control_cols) {
  pred_drops[, c(this_name) := predict(glm(paste0("Z ~ ", paste(setdiff(control_cols, this_name), collapse = " + ")), data = test_dt, family = "binomial"), type = "response")]
}


#======================================================================
# Plot sensitivity
sensitivity_results = as.data.table(sensitivity_results)

max_lambda = 1.5
pdf(file = "Real_data/Results/Figs/RHC_InitialSensitivity_boostOnly_aipwOnly_wide.pdf", width = 6, height = 3)
print(plotSensitivity(lambda_vals = unique(sensitivity_results$Lambda[sensitivity_results$Lambda <= max_lambda]), color_column = "Method", y_add = "", this_infer = "Joint", make_as_line = TRUE))
dev.off()

max_lambda = 2
pdf(file = "Real_data/Results/Figs/RHC_InitialSensitivity_boostOnly_aipwOnly_wide_withCalibration.pdf", width = 6, height = 3.5)
print(addCalibration(plotSensitivity(lambda_vals = unique(sensitivity_results$Lambda[sensitivity_results$Lambda <= max_lambda]), color_column = "Method", y_add = "", this_infer = "Joint", make_as_line = TRUE), max_lambda))
dev.off()



#======================================================================
# Statistics for paper 


statistics_for_paper = c(
  # Unconfounded estimate 
  paste0(
    "Using an AIPW estimator, we estimate an ATE of around $",
    sensitivity_results[Lambda == 1 & Method == "DVDS", round(100 * upper, 1)],
    "\\%$ with a 95\\% confidence interval of $[",
    sensitivity_results[Lambda == 1 & Method == "DVDS", round(100 * lower.CI, 1)], 
    "\\%, ",
    sensitivity_results[Lambda == 1 & Method == "DVDS", round(100 * upper.CI, 1)], 
    "\\%]$"
  ),
  
  # Positive point estimate 
  paste0(
    "($\\Lambda = ",
    sensitivity_results[upper >= 0 & Method == "ZSB",  round(min(Lambda), 2)],
    "$ for the \\textit{ZSB} method, $\\Lambda = ",
    sensitivity_results[upper >= 0 & Method == "DVDS", round(min(Lambda), 2)], 
    "$ for the \\textit{DVDS} method)"
  ),
  
  # Positive CI 
  paste0(
    "odds of treatment by a factor of $\\Lambda = ",
    sensitivity_results[upper.CI >= 0 & Method == "DVDS", round(min(Lambda), 2)], 
    "$ could already reverse the original finding"
  ),
  
  #Crossover point 
  paste0(
    "Our lower confidence bound was $",
    round(quantile(min_lambdas$min_lambda_val, 0.05), 3),
    "$, meaning that it is highly unlikely the true ATE would be negative unless ",
    "unobserved confounders increased the odds of treatment in some covariate level by $> ",
    round(100 * (quantile(min_lambdas$min_lambda_val, 0.05) - 1), 1), 
    "\\%$ or decreased the odds of treatment in some covariate level by $> ",
    round(100 * (1 - 1 / quantile(min_lambdas$min_lambda_val, 0.05)), 1),
    "\\%$"
  )
)


fileConn <- file("Real_data/Results/Tabs/statistics_for_paper.tex")
writeLines(statistics_for_paper, fileConn)
close(fileConn)


estimated_sufficient_OR = sensitivity_results[upper >= 0 & Method == "DVDS", round(min(Lambda), 2)]

fileConn <- file("Real_data/Results/Tabs/statistics_sufficient_OR.tex")
writeLines(c(estimated_sufficient_OR, "%"), fileConn)
close(fileConn)


# Statistics on changing the odds of treatment 
glm_changes = melt(pred_drops, id.vars = c("obs_id", "Z", "pred_SL", "None"))
glm_changes[, odds_ratio_change := (None / (1-None)) / (value / (1-value))]
glm_changes[, LargerChange := abs(log(odds_ratio_change)) >= abs(log(estimated_sufficient_OR))]
fwrite(
  glm_changes[, .(PercLargerChanges = 100 * mean(LargerChange)), by = variable][order(-PercLargerChanges)],
  "Real_data/Results/Tabs/statistics_drops_sufficient.csv"
)



for (this_var in c("ad_neuro", "surv2md1", "pafi1")) {
  
  this_var_changes = glm_changes[variable == this_var]
  x_lab = paste0(this_var, " ORs (", this_var_changes[, round(100 * mean(LargerChange), 1)], "% Ratios >= ", estimated_sufficient_OR, ")")
  
  pdf(file = paste0("Real_data/Results/Figs/ORChanges_", this_var, ".pdf"), width = 6, height = 3)
  print(
    ggplot(this_var_changes, aes(x = odds_ratio_change)) + 
      geom_histogram(bins = 30) + scale_x_log10() + theme_bw() + labs(x = x_lab) +
      geom_vline(xintercept = c(1/estimated_sufficient_OR, estimated_sufficient_OR), lty = 2)
  )
  dev.off()
}





