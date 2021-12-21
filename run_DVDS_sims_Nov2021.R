
#Setting to add sims at all (FALSE to skip to tables/graphs)
rerun_any_sims = FALSE      
#Setting for removing existing sims + starting from scratch
rerun_all_sims = FALSE    
#Setting for loading all sims & condensing 
condense_sims  = FALSE   

#Settings for how many/which simulations to run:
#   (1:num_sims)[(1:num_sims %% num_machines) == mchn_remainder]
#   mchn_remainder: one of 0, ..., num_machines-1: 
#     the index of this machine if manually parallelizing the sims
#   R: mchn_remainder = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
num_sims       = 1000
num_machines   = 24
if (rerun_any_sims) {
  mchn_remainder = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) #0 
}

#Settings for estimation/bootstrap
N = 1000
K = 5
quant_boots  = 500 
Lambdas = seq(1, 3, length.out = 13)

#Settings for estimation
lin_covs   = 'X1 + X2 + X3 + X4 + X5'
full_covs  = 'X1 + X2 + X3 + X4 + X5 + abs(X2) + abs(X5) + X1 * X2'
hmskd_covs = 'X1 + X2 + X3 + X4 + X5 + X1 * X2'
bad_covs   = 'X1 + sin(X2) + X3 + X4 + X5'
KS_covs    = 'XK1 + XK2 + XK3 + XK4 + XK5'


#==========================
# Load Functions 
#==========================

source("DVDS.R")

#==========================
# Some Simulation Functions
#==========================


#Proposition 2 from DG: If Z | X ~ Bern(e(X)) and Y | X, Z ~ N(mu(X, Z), sigma^2(X)), then:
#   ATE bounds are E[mu(X, 1) - mu(X, 0)] \pm (L^2 - 1)/L * phi(Phi^{-1}(tau)) * E[sigma(X)]
e.oracle = function(data) {
  1 / (1 + exp(-data$X1 - data$X2 - data$X1*data$X2))
}

mu.oracle = function(data) {
  data$X1 + data$X3 + data$X1 * data$X2
}

sigma.oracle = function(data) {
  1 - abs(data$X2) + abs(data$X5)
}

q.oracle = function(data, tau) {
  mu.oracle(data) + sigma.oracle(data) * qnorm(tau)
}

#Check function TODO 
kappa.oracle = function(data, Lambda, t) {
  mu.oracle(data) + sigma.oracle(data) * (Lambda - Lambda^(-1)) * dnorm(qnorm(Lambda/(Lambda+1))) * (2*t - 1)
}

dgp = function(n) {
  X1 = runif(n, min=-1, max=1) 
  X2 = runif(n, min=-1, max=1) 
  X3 = runif(n, min=-1, max=1) 
  X4 = runif(n, min=-1, max=1) 
  X5 = runif(n, min=-1, max=1) 
  eps = rnorm(n)
  V   = runif(n, 0, 1)
  data = data.frame(X1, X2, X3, X4, X5)
  if (min(sigma.oracle(data)) < 0) {
    stop('Invalid SD formula (must have positive SD)')
  }
  data$Y = mu.oracle(data) + sigma.oracle(data) * eps
  data$Z = (V <= e.oracle(data))
  #Kang and Schafer-type covariates 
  data$XK1 = exp(data$X1/2)
  data$XK2 = 10 + data$X2 / (1 + exp(data$X1))
  data$XK3 = (data$X1 * data$X3 / 25 + 0.6)^3
  data$XK4 = (data$X2 + data$X4 + 20)^2
  data$XK5 = (data$X5)^5
  return(data)
}

#Currently, E[rsd] = 1
mean_sigma = 1
truebounds = ((Lambdas^2-1)/Lambdas) * dnorm(qnorm(Lambdas/(Lambdas+1))) * (mean_sigma)

if (t.test(sigma.oracle(dgp(10^6)) - mean_sigma)$p.value < 0.01) {
  warning("Potentially incorrect truebounds")
}


#======================
# Some special functions for oracle estimates
#======================


bootstrapPropsOracle = function(data) {
  
  start_boot_seed = .Random.seed 
  
  boot_data = mclapply(1:quant_boots, function(b) {
    set.seed(b)
    
    cvboot  = sample(1:5, nrow(data), replace = TRUE)
    realind = sample(1:nrow(data), nrow(data), replace = TRUE)
    cvprop  = e.oracle(data[realind,])
    to_add = merge(data.frame(realind, bootind = 1:nrow(data), prop = cvprop), data.frame(realind = 1:nrow(data), cv = cvboot), by = "realind", all.x = TRUE)
    return(to_add)
  })
  .Random.seed = start_boot_seed 
  
  return(boot_data)
}

dvds.oracle = function(Lambdas, data, cubic_err = F, form_y = 'Y', form_z = 'Z') {
  boot_data = bootstrapPropsOracle(data)
  return(foreach(Lambda=Lambdas, .combine=rbind) %do% {
    eif = list()
    boot_estim = list()
    se_list = list()
    e = e.oracle(data)
    if (cubic_err) {
      e = e + e * ( (e + (1-e) * exp(rnorm(nrow(data), sd=1000^(-1/3))))^(-1) - 1)
      e = pmin(pmax(e, 0.01), 0.99)
      for (i in 1:length(boot_data)) {
        boot_data[[i]]$prop = e[boot_data[[i]]$realind]
      }
      q_err     = rnorm(nrow(data), sd = nrow(data)^(-1/3))
      kappa_err = rnorm(nrow(data), sd = nrow(data)^(-1/3))
    } else {
      q_err   = 0
      kappa_err = 0
    }
    for (z in 0:1) {
      ee = (2*z-1)*e + (1 - z)
      for (t in 0:1) {
        zz  = (2*z-1)*data$Z + (1 - z)
        tau = ((1-t)+t*Lambda)/(Lambda + 1)
        q   = q.oracle(data, tau) + q_err
        kappa = kappa.oracle(data, Lambda, t) + kappa_err 
        eif[[paste(z,t)]] = kappa + zz*(q - kappa)/ee + zz*(data$Y - q)*(ee + (1-ee)*Lambda^((2*t-1)*sign(data$Y - q)))/ee
        
        se_list[[paste(z,t)]] = sd(eif[[paste(z, t)]])/sqrt(nrow(data))
        
        boot_estims_wrongdirection = mclapply(1:quant_boots, function(b) {
          set.seed(b)
          ind_star       = boot_data[[b]]$realind
          data_star      = data[ind_star, ]
          data_star$prop = boot_data[[b]]$prop
          q_star     =     q[ind_star]
          kappa_star = kappa[ind_star]
          e_star   = data_star$prop
          Y_star   = data_star[[form_y]]
          zz_star  = (2*z-1)*data_star[[form_z]] + (1 - z)
          ee_star  = (2*z-1)*e_star + (1 - z)
          c(
            mean(kappa_star + zz_star*(q_star - kappa_star)/ee_star + zz_star*(data_star[[form_y]] - q_star)*(ee_star + (1-ee_star)*Lambda^((2*t-1)*sign(data_star[[form_y]] - q_star)))/ee_star),
            mean(Y_star),
            mean(zz_star)
          )
        })
        
        
        for (b in 1:quant_boots) {
          boot_estim[[paste(z,t)]][b]       = boot_estims_wrongdirection[[b]][1]
          boot_estim[[paste(z,'Ymean')]][b] = boot_estims_wrongdirection[[b]][2]
          boot_estim[[paste(z,'Zmean')]][b] = boot_estims_wrongdirection[[b]][3]
        }
      }
    }
    
    #ATE standard errors
    se_list[['1 1 0 0']] = sd(eif[['1 1']]-eif[['0 0']])/sqrt(nrow(data))
    se_list[['1 0 0 1']] = sd(eif[['1 0']]-eif[['0 1']])/sqrt(nrow(data))
    
    boot_facts = summarizeBoots(boot_estim)
    
    summarizeResults(Lambda, data, form_y, form_z, eif, se_list, boot_facts) 
  })
}

qb.oracle = function(Lambdas, data, cubic_err = F, form_y = 'Y', form_z = 'Z') {
  data = data%>%arrange(!! sym('Y'))
  e    = e.oracle(data)
  boot_data = bootstrapPropsOracle(data)
  if (cubic_err) {
    e = e + e * ( (e + (1-e) * exp(rnorm(nrow(data), sd=1000^(-1/3))))^(-1) - 1)
    e = pmin(pmax(e, 0.01), 0.99)
    for (i in 1:length(boot_data)) {
      boot_data[[i]]$prop = e[boot_data[[i]]$realind]
    }
    q_err   = rnorm(nrow(data), sd = nrow(data)^(-1/3))
  } else {
    q_err   = 0
    kappa_err = 0
  }
  return(foreach(Lambda=Lambdas, .combine=rbind) %do% {
    iif = list()
    se_list = list()
    boot_estim = list()
    for (z in 0:1) {
      zz = (2*z-1)*data$Z + (1 - z)
      ee = (2*z-1)*e + (1 - z)
      for (t in 0:1) {
        #First-stage quantiles 
        tau = ((1-t)+t*Lambda)/(Lambda + 1)
        q = q.oracle(data, tau) + q_err
        Q = q
        #Get quantiles for quantile balancing 
        quant_df = data.frame(Y = data$Y, Q = q)
        if (sd(quant_df$Q) != 0) {
          q = predict(rq(Y ~ Q, tau=tau, weights=((1-ee)/ee)[zz==1], data = quant_df[zz==1,]), quant_df)
        } else {
          q = predict(rq(Y ~ 1, tau=tau, weights=((1-ee)/ee)[zz==1], data = quant_df[zz==1,]), quant_df)
        }
        #Return i's "influence function," the quantile balancing term 
        # (TODO: does this produce valid standard errors?)
        iif[[paste(z,t)]] = (zz*q/ee + zz*(data$Y - q)*(ee + (1-ee)*Lambda^((2*t-1)*sign(data$Y - q)))/ee) / mean(zz/ee)
        
        se_list[[paste(z,t)]] = sd(iif[[paste(z, t)]])/sqrt(nrow(data))
        
        boot_estims_wrongdirection = mclapply(1:quant_boots, function(b) {
          set.seed(b)
          ind_star       = boot_data[[b]]$realind
          data_star      = data[ind_star, ]
          data_star$prop = boot_data[[b]]$prop
          q_star   = Q[ind_star]
          e_star   = data_star$prop
          Y_star   = data_star[[form_y]]
          zz_star  = (2*z-1)*data_star[[form_z]] + (1 - z)
          ee_star  = (2*z-1)*e_star + (1 - z)
          
          quant_df_star  = data.frame(Y = data_star[[form_y]], Q = q_star)
          if (sd(q_star) == 0) {
            q_balance_star = predict(rq(Y ~ 1, tau=tau, weights=((1-ee_star)/ee_star)[zz_star==1], data = quant_df_star[zz_star==1,]), quant_df_star)
          } else {
            q_balance_star = predict(rq(Y ~ Q, tau=tau, weights=((1-ee_star)/ee_star)[zz_star==1], data = quant_df_star[zz_star==1,]), quant_df_star)
          }
          c(
            mean((zz_star*q_balance_star/ee_star + zz_star*(data_star[[form_y]] - q_balance_star)*(ee_star + (1-ee_star)*Lambda^((2*t-1)*sign(data_star[[form_y]] - q_balance_star)))/ee_star) / mean(zz_star/ee_star)),
            mean(Y_star),
            mean(zz_star)
          )
        })
        
        for (b in 1:quant_boots) {
          boot_estim[[paste(z,t)]][b]       = boot_estims_wrongdirection[[b]][1]
          boot_estim[[paste(z,'Ymean')]][b] = boot_estims_wrongdirection[[b]][2]
          boot_estim[[paste(z,'Zmean')]][b] = boot_estims_wrongdirection[[b]][3]
        }
        
      }
    }
    
    #ATE standard errors
    se_list[['1 1 0 0']] = sd(iif[['1 1']]-iif[['0 0']])/sqrt(nrow(data))
    se_list[['1 0 0 1']] = sd(iif[['1 0']]-iif[['0 1']])/sqrt(nrow(data))
    
    boot_facts = summarizeBoots(boot_estim)
    
    summarizeResults(Lambda, data, form_y, form_z, iif, se_list, boot_facts)
  })
}

zsb.oracle = function(Lambdas, data, cubic_err = F, form_y = 'Y', form_z = 'Z', center_mu = F) {
  data = data%>%arrange(!! sym('Y'))
  e    = e.oracle(data)
  boot_data = bootstrapPropsOracle(data)
  mu = list("1" = rep(0, nrow(data)), "0" = rep(0, nrow(data)))
  if (center_mu) {
    for (z in 0:1) {
      mu[[as.character(z)]] = mu.oracle(data)
    }
    rm(z)
  }
  if (cubic_err) {
    e = e + e * ( (e + (1-e) * exp(rnorm(nrow(data), sd=1000^(-1/3))))^(-1) - 1)
    e = pmin(pmax(e, 0.01), 0.99)
    for (i in 1:length(boot_data)) {
      boot_data[[i]]$prop = e[boot_data[[i]]$realind]
    }
  } 
  return(foreach(Lambda=Lambdas, .combine=rbind) %do% {
    iif = list()
    boot_estim = list()
    se_list = list()
    for (z in 0:1) {
      data_mu    = mu[[as.character(z)]]
      data_resid = data$Y - data_mu
      zz = (2*z-1)*data$Z + (1 - z)
      ee = (2*z-1)*e + (1 - z)
      #Inverse of lower propensity (higher weight)
      L  = ee / (ee + (1-ee)*Lambda^(+1))
      U  = ee / (ee + (1-ee)*Lambda^(-1))
      for (t in 0:1) {
        zsb_thresh  = optimCutHajek(data_resid, zz, L, U, run_min = as.logical(1-t))
        opt_invprop = (ee + (1-ee)*Lambda^((2*t-1)*sign(data_resid - zsb_thresh)))/ee
        iif[[paste(z,t)]] = data_mu + (zz * data_resid * opt_invprop) / mean(zz * opt_invprop)
        
        se_list[[paste(z,t)]] = sd(iif[[paste(z, t)]])/sqrt(nrow(data))
        
        boot_estims_wrongdirection = mclapply(1:quant_boots, function(b) {
          set.seed(b)
          ind_star       = boot_data[[b]]$realind
          data_star      = data[ind_star, ]
          data_star$prop = boot_data[[b]]$prop
          e_star   = data_star$prop
          Y_star   = data_star[[form_y]]
          mu_star  = data_mu[ind_star]
          resid_star = Y_star - mu_star
          zz_star  = (2*z-1)*data_star[[form_z]] + (1 - z)
          ee_star  = (2*z-1)*e_star + (1 - z)
          L_star  = ee_star / (ee_star + (1-ee_star)*Lambda^(+1))
          U_star  = ee_star / (ee_star + (1-ee_star)*Lambda^(-1))
          zsb_thresh_star  = optimCutHajek(resid_star, zz_star, L_star, U_star, run_min = as.logical(1-t))
          opt_invprop_star = (ee_star + (1-ee_star)*Lambda^((2*t-1)*sign(resid_star - zsb_thresh_star)))/ee_star
          
          c(
            mean(Y_star + (zz_star * resid_star * opt_invprop_star) / mean(zz_star * opt_invprop_star)),
            mean(Y_star),
            mean(zz_star)
          )
        })
        
        for (b in 1:quant_boots) {
          boot_estim[[paste(z,t)]][b]       = boot_estims_wrongdirection[[b]][1]
          boot_estim[[paste(z,'Ymean')]][b] = boot_estims_wrongdirection[[b]][2]
          boot_estim[[paste(z,'Zmean')]][b] = boot_estims_wrongdirection[[b]][3]
        }
        
      }
    }
    #ATE standard errors
    se_list[['1 1 0 0']] = sd(iif[['1 1']]-iif[['0 0']])/sqrt(nrow(data))
    se_list[['1 0 0 1']] = sd(iif[['1 0']]-iif[['0 1']])/sqrt(nrow(data))
    
    boot_facts = summarizeBoots(boot_estim)
    summarizeResults(Lambda, data, form_y, form_z, iif, se_list, boot_facts)
    
  })
  
}




#==========================
# Prep Files/Solve For Sims to Run
#   Sims are saved in Results/estimates_[number].csv
#==========================


if (rerun_any_sims) {
  potential_sims = (1:num_sims)[(1:num_sims %% num_machines) == mchn_remainder]
  if (rerun_all_sims) {
    files_to_delete = as.numeric(gsub(".csv", "", gsub("estimates_", "", list.files("Results/Sims/")[grepl("estimates_", list.files("Results/Sims/"))], fixed = T), fixed = T))
    if (length(files_to_delete) > 0) {
      files_to_delete = paste0("estimates_", files_to_delete[files_to_delete %% num_machines == mchn_remainder], ".csv")
      for (el in files_to_delete) {
        file.remove(paste0("Results/Sims/", el))
      }
    }
    rm(el, files_to_delete)
  }
  
  run_sims = setdiff(potential_sims, as.numeric(gsub(".csv", "", gsub("estimates_", "", list.files("Results/Sims/")[grepl("estimates_", list.files("Results/Sims/"))], fixed = T), fixed = T)))
} else {
  run_sims = c()
}

#==========================
# Specs to Run
#==========================

#What specs to run? 
spec_list = list()

for (this_formula in c("zsba", "zsb", "qb", "dvds")) {
  
  spec_list[[length(spec_list) + 1]] = paste0("results.", this_formula, ".cubicErr", paste(rep(" ", 9 - nchar(this_formula)), collapse = ""), "= ", gsub("zsba", "zsb", this_formula), ".oracle(", paste(rep(" ", 7 - nchar(this_formula)), collapse = ""), "Lambdas, data", ifelse(this_formula == "zsba", ", center_mu = T", ""), ", cubic_err = T)")
  spec_list[[length(spec_list) + 1]] = paste0("results.", this_formula, ".oracle  ", paste(rep(" ", 9 - nchar(this_formula)), collapse = ""), "= ", gsub("zsba", "zsb", this_formula), ".oracle(", paste(rep(" ", 7 - nchar(this_formula)), collapse = ""), "Lambdas, data", ifelse(this_formula == "zsba", ", center_mu = T", ""), ")")
  
  #for (this_component in c("includeabs", "svm", "boost", "justlinear", "justhomosked", "KSmisspec", "KSbstkappa")) {
  for (this_component in c("includeabs", "svm", "boost", "justlinear", "justhomosked")) {
    
    #Pass covariates: Including absolute value only for includeabs, include sin(x2) only for wrongx2
    these_covs   = ifelse(this_component == "includeabs", "full_covs", ifelse(this_component == "justhomosked", "hmskd_covs", ifelse(this_component == "wrongx2", " bad_covs", ifelse(this_component == "KSmisspec" | this_component == "KSbstkappa", "  KS_covs", " lin_covs"))))
    add_function = paste0(gsub("zsba", "zsb", this_formula), "(", paste(rep(" ", 14 - nchar(this_formula)), collapse = ""), "Lambdas, data, ", these_covs, ", ")
    
    #Propensity model: Boost, correct spec (true and partial linear), or misspecified with missing absolute value 
    add_function = paste0(add_function, "'Z', 'Y', ", ifelse(this_component == "boost", "boostregn, boostregn_option_bin, ", ifelse(this_component == "svm", "  svmregn,   svmregn_option_bin, ", "  linregn,   linregn_option_bin, ")))
    
    #muhat model for ZSB (AIPW)
    if (this_formula == "zsba") {
      add_function = paste0(add_function, ifelse(this_component %in% c("boost", "boostokprop", "ptlin", "KSbstkappa"), "boostregn,  boostregn_option_cts, ", ifelse(this_component == "svm", "  svmregn,    svmregn_option_cts, ", "  linregn,    linregn_option_cts, ")))
      add_function = paste0(add_function, "semiadaptive=T, ")
    }
    
    if (this_formula != "zsb" & this_formula != "zsba") {
      #Quantiles: linear outside of boost specifications 
      add_function = paste0(add_function, 
                            ifelse(this_component %in% c("boost", "svm", "KSbstkappa", "boostokprop"),
                                   "forestquant, list(reuse_forests=T), ",
                                   ifelse(this_component == "linconst", 
                                          "constquant, constquant_option,      ", 
                                          "  linquant,   linquant_option,      "
                                                        )))
      
      #Kappa-hat estimation 
      if (this_formula == "dvds") {
        add_function = paste0(add_function, ifelse(this_component %in% c("boost", "boostokprop", "ptlin", "KSbstkappa"), "boostregn,  boostregn_option_cts, ", ifelse(this_component == "svm", "  svmregn,    svmregn_option_cts, ", "  linregn,    linregn_option_cts, ")))
        add_function = paste0(add_function, "semiadaptive=T, ")
        if (this_component %in% c("boost", "boostokprop", "ptlin", "KSbstkappa", "svm")) {
          add_function = paste0(add_function, "form_x_kappa = paste0('Q + ', ", ifelse(this_component == "KSbstkappa", these_covs, "lin_covs"), "), ")
        } else {
          add_function = paste0(add_function, paste0("form_x_kappa = paste0('Q + ', ", these_covs, "), "))
        }
      }
      
      add_function = paste0(add_function, "form_x_quant = ", these_covs, ", ")
      
    }
    
    #Add bootstrap SEs 
    add_function = paste0(
      add_function, 
      "boot_infer = T, boot_settings = list(reuse_boot = T, refit_propensities = ",
      ifelse(this_component %in% c("boost", "svm"), "F", "T"), 
      ", B = ", quant_boots, ", return_seed = TRUE), "
    )
    
    spec_list[[length(spec_list) + 1]] = paste0(
      "results.", this_formula, ".", this_component, paste(rep(" ", 12 - nchar(this_component) + 4 - nchar(gsub("zsba", "zsb", this_formula))), collapse = ""), " = ", 
      add_function, "trim.type='clip', normalize = F)"
    )
  }
}

print("The following are the specifications to run:")
for (el in spec_list) {
  print(el)
}





#==========================
# Run+Save Sims
#==========================

for (this_sim in run_sims) {
  
  #print(paste0("Starting simulation ", this_sim))
  
  set.seed(this_sim)
  
  #Simulate data
  data = dgp(N)
  
  #Simulate bootstrap data + pre-compute overall and bootstrapped propensities 
  start_time = proc.time()[3]
  
  #Ensure same seed in all estimates to avoid randomness in cross-fitting
  estim_seed = .Random.seed
  
  #Run all specifications 
  for (el in spec_list) {
    
    if (this_sim == run_sims[1]) {
      print(paste0("Starting ", strsplit(el, "=")[[1]][1]))
    }
    
    .Random.seed = estim_seed
    eval(parse(text=el))
    
    if (this_sim == run_sims[1]) {
      if (length(warnings()) > 0) {
        print(warnings())
        assign("last.warning", NULL, envir = baseenv())
      }
    }
  }
  end_time = proc.time()[3]
  print(paste0("All done with sim ", this_sim, ". Time elapsed: ", round((end_time - start_time)/60, 2), " minutes."))
  
  
  #Combine results for saving 
  all_results_combined = data.frame()
  for (these_results in ls()[grepl("results\\.", ls())]) {
    to_add = get(these_results) %>% mutate(method = gsub("results.", "", these_results, fixed = T))
    all_results_combined = rbind(all_results_combined, to_add)
  }
  
  write.csv(all_results_combined, paste0("Results/Sims/estimates_", this_sim, ".csv"))
  
}

print("All done with simulations.")




#==========================
# APO Bounds by Simulation
#==========================

#APO bounds (from Dorn and Guo, Corollary 3):
# E[Y(1)]: E[\mu(X)] \pm APO_lambda_factor * E[(1-Z) \sigma(X)]
# E[Y(0)]: E[\mu(X)] \pm APO_lambda_factor * E[ (Z)  \sigma(X)]

#TODO: make this analytic 
if (condense_sims) {
  temp_df = dgp(10^7)
  APO_0_factor = mean(  e.oracle(temp_df)   * sigma.oracle(temp_df))
  APO_1_factor = mean((1-e.oracle(temp_df)) * sigma.oracle(temp_df))
  rm(temp_df)
  mu_mean = 0
  APO_lambda_factor = ((Lambdas - 1/Lambdas) * dnorm(qnorm(Lambdas/(Lambdas + 1))))
  
  #Assuming E[Z] = 1/2 (TODO: Add checks)
  true_bound_df = rbind(
    data.frame(Lambda=Lambdas, estimand='ATE', upper=T, truth= truebounds), 
    data.frame(Lambda=Lambdas, estimand='ATE', upper=F, truth=-truebounds),
    data.frame(Lambda=Lambdas, estimand='1',   upper=T, truth=mu_mean+APO_lambda_factor * APO_1_factor), 
    data.frame(Lambda=Lambdas, estimand='1',   upper=F, truth=mu_mean-APO_lambda_factor * APO_1_factor),
    data.frame(Lambda=Lambdas, estimand='0',   upper=T, truth=mu_mean+APO_lambda_factor * APO_0_factor), 
    data.frame(Lambda=Lambdas, estimand='0',   upper=F, truth=mu_mean-APO_lambda_factor * APO_0_factor),
    data.frame(Lambda=Lambdas, estimand='ATT', upper=T, truth=+APO_lambda_factor * APO_0_factor / (1/2)), 
    data.frame(Lambda=Lambdas, estimand='ATT', upper=F, truth=-APO_lambda_factor * APO_0_factor / (1/2)),
    data.frame(Lambda=Lambdas, estimand='ATC', upper=T, truth=+APO_lambda_factor * APO_1_factor / (1/2)), 
    data.frame(Lambda=Lambdas, estimand='ATC', upper=F, truth=-APO_lambda_factor * APO_1_factor / (1/2))
  )
  
}
#==========================
# Load Sims
#==========================

#Load table of results (using data.table for faster reading)
library(data.table)

if (condense_sims) {
  all_results = list()
  for (spec in list.files("Results/Sims/")[list.files("Results/Sims/") %like% "estimates"]) {
    all_results[[length(all_results) + 1]] = fread(paste0("Results/Sims/", spec))[, Sim := gsub("estimates_", "", strsplit(spec, "\\.")[[1]][1])]
  }
  all_results = rbindlist(all_results)[, ':='(V1 = NULL, nobs = N)]
  all_results$Lambda2   = round(  all_results$Lambda, 4)
  true_bound_df$Lambda2 = round(true_bound_df$Lambda, 4)
  all_results = merge(
    all_results,
    true_bound_df[, c("Lambda2", "estimand", "upper", "truth")],
    by = c('Lambda2', 'estimand', 'upper'), all.x = TRUE
  )
  all_results[, Lambda2 := NULL]
  
  #Coverage
  all_results[upper == TRUE,  AlphaIIF  := 1 - pnorm((truth - estimate) / sterr_iif)]
  all_results[upper == FALSE, AlphaIIF  :=     pnorm((truth - estimate) / sterr_iif)]
  all_results[upper == TRUE,  AlphaBoot := 1 - pnorm((truth - estimate) / sterr_boot)]
  all_results[upper == FALSE, AlphaBoot :=     pnorm((truth - estimate) / sterr_boot)]
  
  setkeyv(all_results, c("estimand", "method", "nobs", "Lambda", "upper", "AlphaIIF"))
  
  fwrite(all_results[, .(estimand, method, nobs, Lambda, upper, Sim, estimate, AlphaIIF, AlphaBoot, truth, quantiles90, quantiles95)], "Results/Sims/combined_table.csv")
  
  #Pre-process for visualization: Point estimates 
  fwrite(all_results[!is.na(truth), .(estimand, method, Lambda, upper, Sim, estimate, truth)], "combined_table_pointEstims.csv")
  
  #Pre-process for visualization: Coverage 
  fwrite(all_results[, .(estimand, method, Lambda, upper, Sim, AlphaIIF, AlphaBoot)], "combined_table_coverageStats.csv")
} else {
  all_results = fread("Results/Sims/combined_table.csv")
}


setkeyv(all_results, c("estimand", "method", "nobs", "Lambda", "upper", "AlphaIIF"))
all_results[, CoverageIIF  := (.N:1) / .N, by = c("estimand", "method", "nobs", "Lambda", "upper")]

setkeyv(all_results, c("estimand", "method", "nobs", "Lambda", "upper", "AlphaBoot"))
all_results[, CoverageBoot := (.N:1) / .N, by = c("estimand", "method", "nobs", "Lambda", "upper")]

#TODO: Add true value percentile bootstrap (need to calculate in sims)

all_results[, ':='(
  Coverage95CIBelowIIF  = qbinom(0.025, .N, 1-AlphaIIF) / .N,
  Coverage95CIAboveIIF  = qbinom(0.975, .N, 1-AlphaIIF) / .N,
  Coverage95CIBelowBoot = qbinom(0.025, .N, 1-AlphaBoot) / .N,
  Coverage95CIAboveBoot = qbinom(0.975, .N, 1-AlphaBoot) / .N
), by = c("estimand", "method", "nobs", "Lambda", "upper")]


all_results[, CoverageTestIIF   := ifelse(CoverageIIF   < Coverage95CIBelowIIF   | CoverageIIF   > Coverage95CIAboveIIF,  "Reject (Z-Test)", "Accept (Z-Test)")]
all_results[, CoverageTestBoot  := ifelse(CoverageBoot  < Coverage95CIBelowBoot  | CoverageBoot  > Coverage95CIAboveBoot, "Reject (Z-Test)", "Accept (Z-Test)")]
#TODO: Add this once we keep the percentile bootstrap p-value on the truth in the sims
#all_results[, CoverageTestPBoot := ifelse(CoveragePBoot < Coverage95CIBelowPBoot | CoveragePBoot > Coverage95CIAbovePBoot, "Reject (Z-Test)", "Accept (Z-Test)")]

#For comparison to oracle-only sims 
#sims_oracleOnly = fread("../EfficiencyAlgorithms/DVDS_JustPointEstimates/Results/Tabs/oracle_only_sims.csv")[HetVar==TRUE & HetMean==TRUE & HetProp==TRUE & AbsProp == FALSE, .(Lambda, upper, Truth, Seed = as.character(Seed), estimate_oo = estimate, sterr_oo = sterr)]
#sims_thisSim    = all_results[estimand == "ATE" & method == "dvds.oracle" & Lambda %in% 1:4, .(Lambda, upper, Seed = Sim, estimate_dvds = estimate, sterr_dvds = sterr)]
#sims_comparison = merge(sims_oracleOnly, sims_thisSim, by = c("Lambda", "upper", "Seed"))[, .(Lambda, upper, Truth, Seed, estimate_oo, estimate_dvds, sterr_oo, sterr_dvds)][order(as.numeric(Seed))]


#Rename inputs (from app)
list_mapping = list(
  "Oracle" = "oracle",
  "Oracle + n^(-1/3) Error" = "cubicErr",
  "MLE (Well-Specified)" = "includeabs",
  "MLE (Misspecified Quantiles)" = "justhomosked",
  "MLE (No-Interaction Quantiles + Propensities)" = "justlinear",
  "MLE (Misspecified Quantiles + Propensities)" = "wrongx2",
  "Boost (Forest Quantiles)" = "boost",
  "Kang and Schafer" = "KSmisspec",
  "Kang and Schafer (Boost Kappa + Forest Quantiles)" = "KSbstkappa",
  "SVM (Forest Quantiles)" = "svm"
)
method_names = list()
for (i in 1:length(list_mapping)) {
  method_names[[i]] = data.table(Inputs = list_mapping[[i]], InputsPrint = names(list_mapping)[i])
}
method_names = rbindlist(method_names)



true_values = dcast(all_results[, .(truth = truth[1]), by = .(estimand, Lambda, upper)], estimand + Lambda ~ upper, value.var = "truth")
setnames(true_values, c("FALSE", "TRUE"), c("Lower", "Upper"))

#Get long table of results
all_results_long = dcast(all_results, estimand + Lambda + method + Sim ~ upper, value.var = "estimate")
setnames(all_results_long, c("FALSE", "TRUE"), c("Lower", "Upper"))

#Add nuisance names 
formatTab = function(dt) {
  dt = copy(dt)
  dt[, ':='(Estimator = strsplit(method, "\\.")[[1]][1], Inputs = strsplit(method, "\\.")[[1]][2]), by = method]
  dt = merge(dt, method_names, by = "Inputs", all.x = TRUE)
  
  #Format order of methods for plotting
  dt[, MethodWrite := paste0(InputsPrint, " (", Estimator, ")")]
  #Ensures order of choices is order of plotted point estimates 
  dt[, InputsPrint := factor(InputsPrint, levels = method_names$InputsPrint)]
  setkeyv(dt, c("InputsPrint", "Estimator"))
  dt[, method := factor(method, levels = unique(dt$method))]
  
  #For printing methods
  dt[, InputsPaper := ifelse(Inputs == "oracle", "Oracle", ifelse(Inputs == "cubicErr", "Oracle + n^(-1/3) Error", ifelse(Inputs == "boost", "Boost", ifelse(Inputs == "includeabs", "Parametric", ifelse(Inputs == "svm", "Machine Learning", "Not in Paper Currently")))))]
  dt[, InputsPaper := factor(InputsPaper, levels = c("Oracle", "Oracle + n^(-1/3) Error", "Parametric", "NP", "Boost", "Machine Learning", "Not in Paper Currently"))]
  
  return(dt)
}

all_results_long = formatTab(all_results_long)



coveragePlot = function(this_lam, this_estimand, these_inputs = c("oracle", "svm", "includeabs"), these_estimators = c("zsb", "zsba", "dvds", "qb")) {
  
  this_dt     = copy(all_results_long)[estimand == this_estimand & round(Lambda, 3) == round(as.numeric(this_lam), 3) & Inputs %in% these_inputs & Estimator %in% these_estimators]
  true_values = true_values[estimand == this_estimand & round(Lambda, 3) == round(as.numeric(this_lam), 3)]
  
  this_dt = this_dt[InputsPrint %in% c("Oracle", "Oracle + n^(-1/3) Error", "MLE (Well-Specified)", "SVM (Forest Quantiles)")]
  
  this_dt[, Estimator := ifelse(Estimator == 'zsb', 'ZSB-IPW', ifelse(Estimator == 'zsba', 'ZSB-AIPW', toupper(Estimator)))]
  
  return(
    ggplot(this_dt) +
      geom_boxplot(aes(x = Estimator, y = Lower, col = InputsPrint)) +
      geom_boxplot(aes(x = Estimator, y = Upper, col = InputsPrint)) +
      facet_wrap(~InputsPaper, ncol = length(unique(this_dt$InputsPaper))) +
      theme_bw() +
      theme(legend.position = "None") +
      geom_hline(yintercept = true_values$Lower, lty=2) +
      geom_hline(yintercept = true_values$Upper, lty=2) +
      xlab("") +
      ylab("Bounds") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  )
}
pdf(file = "Results/Figs/DVDS_Simulation_PointEstimates.pdf", width = 8, height = 5)
print(coveragePlot(2, "ATE"))
dev.off()


#=========================
# Simulated Coverage
#=========================

library(xtable)

#TODO: Make this be bootstrap true p-values (need to recover from the sims)
coverage_rates = all_results[, .(
  AlphaTwoSidedIIF  = min(AlphaIIF),
  AlphaTwoSidedBoot = min(AlphaBoot),
  Cover90PBoot = min( (2 * upper - 1) * (quantiles90 - truth) > 0),
  Cover95PBoot = min( (2 * upper - 1) * (quantiles95 - truth) > 0)
), by = .(estimand, Lambda, nobs, method, Sim)]


coverage_rates = formatTab(coverage_rates[, .(
  Cover95IIF   = mean(AlphaTwoSidedIIF   >= (1-.95)/2),
  Cover90IIF   = mean(AlphaTwoSidedIIF   >= (1-.90)/2),
  Cover95Boot  = mean(AlphaTwoSidedBoot  >= (1-.95)/2),
  Cover90Boot  = mean(AlphaTwoSidedBoot  >= (1-.90)/2),
  Cover95PBoot = mean(Cover95PBoot),
  Cover90PBoot = mean(Cover90PBoot)
), by = .(estimand, Lambda, nobs, method)][, ':='(
  Cover95Joint = ifelse(method %like% "dvds", Cover95IIF, Cover95PBoot),
  Cover90Joint = ifelse(method %like% "dvds", Cover90IIF, Cover90PBoot)
)])

if (condense_sims) {
  coverage_rates[, nobs := all_results$nobs[1]]
  fwrite(coverage_rates, "coverage_rates.csv")
}


getCoverage = function(dt, this_estimand, these_lambdas, these_inputs = c("oracle", "svm", "includeabs"), these_estimators = c("zsb", "zsba", "dvds", "qb"), infer_type) {
  to_ret = dt[estimand == this_estimand & round(Lambda, 3) %in% round(these_lambdas, 3) & Inputs %in% these_inputs & Estimator %in% these_estimators]
  
  #Flip Wide + Long
  to_ret = melt(to_ret[, .(Lambda, Nuisances = InputsPaper, Estimator = toupper(Estimator), Cover95 = get(paste0("Cover95", infer_type)), Cover90 = get(paste0("Cover90", infer_type)))], id.vars = c("Lambda", "Nuisances", "Estimator"))
  to_ret[, ':='(variable = paste0(gsub("Cover", "", variable), "\\%"), value = paste0(round(100*value, 1), "\\%"))]
  to_ret = dcast(to_ret, Lambda + Nuisances + variable ~ Estimator, value.var = "value")
  setnames(to_ret, "variable", "CI")
  keep_cols = intersect(c("Lambda", "Nuisances", "CI", toupper(these_estimators)), colnames(to_ret))
  to_ret = to_ret[, mget(keep_cols)]
  if (length(these_lambdas) == 1) {
    to_ret[, Lambda := NULL]
  }
  if ('ZSB' %in% names(to_ret)) {
    setnames(to_ret, 'ZSB', 'ZSB-IPW')
  }
  if ('ZSBA' %in% names(to_ret)) {
    setnames(to_ret, 'ZSBA', 'ZSB-AIPW')
  }
  return(print(xtable(to_ret), include.rownames = F, sanitize.text.function = function(x) gsub("n^(-1/3) ", "", gsub("Error", "$n^\\{-1/3}$ Error", x), fixed = T)))
}

#Print to latex: ATE for 2
for (this_infer in c("IIF", "Boot", "PBoot", "Joint")) {
  coverage_results = strsplit(getCoverage(coverage_rates, "ATE", 2, infer_type = this_infer), "\n")[[1]]
  coverage_start   = (1:length(coverage_results))[coverage_results %like% "\\{tabular\\}"][1]
  coverage_end     = (1:length(coverage_results))[coverage_results %like% "\\{tabular\\}"][2]
  coverage_results = coverage_results[coverage_start:coverage_end]
  
  fileConn <- file(paste0("Results/Tabs/coverage_ATE_Lambda2_", this_infer, ".tex"))
  writeLines(coverage_results, fileConn)
  close(fileConn)
  
  if (this_infer == "Joint") {
    fileConn <- file("Results/Tabs/coverage_ATE_Lambda2.tex")
    writeLines(coverage_results, fileConn)
    close(fileConn)
  }
  
}

