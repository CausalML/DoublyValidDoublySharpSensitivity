
#===========================
# Load Necessary Packages
#===========================
packages.needed = c('foreach','tidyverse','quantreg','randomForest', 'quantreg','gbm', 'parallel', 'grf', 'e1071');
lapply(packages.needed, library, character.only = TRUE);

#===========================
# Quick Helper Functions
#===========================

# split [n] into K folds; return n-vector of fold membership id
make.cvgroup = function(n, K, right = TRUE) {
  split     = runif(n)
  return(as.numeric(cut(split, quantile(split, probs = seq(0, 1, 1/K)), include.lowest = TRUE, right = right)))
}

# as above but conditions on fraction of 1s to 0s being similar across folds
make.cvgroup.balanced = function(data, K, form_z) {
  cvgroup = numeric(nrow(data))
  cvgroup[data[[form_z]]==1] = make.cvgroup(sum(data[[form_z]]==1), K, right = TRUE)
  cvgroup[data[[form_z]]==0] = make.cvgroup(sum(data[[form_z]]==0), K, right = FALSE)
  return(cvgroup)
}

#===========================
# Medium Wrapper Functions
#===========================

# returns list with attributes prop and keep
# prop is the cross-fitted propensities
# when trim.type='drop', keep is which units are within the trim interval, otherwise keep is all T
# when trim.type='clip', clip propensities to the trim interval
# normalize = "0" to multiply odds to get E_n[(1-Z) / (1-\hat{e}(X))] = 1 
#   normalize = "1" to multiply odds to get E_n[Z / \hat{e}(X)] = 1 
#   I.e. imposes mean( (Z == normalize) / \hat{E}[Z==normalize|X]) = 1
cross.fit.propensities = function(data, cvgroup, form_x, form_z, method_prop, option_prop, trim=c(0.01,0.99), trim.type='clip', normalize="NA") {
  K = max(cvgroup)
  prop = numeric(nrow(data))
  for (k in 1:K) {
    prop[cvgroup==k] = method_prop(data, cvgroup!=k, cvgroup==k, form_x, form_z, option_prop)
  }
  if(trim.type == 'drop') {
    keep = (prop>trim[1] & prop<trim[2])
    prop[!keep] = 0.5
  } else {
    keep = rep(T, nrow(data))
  }
  if(trim.type == 'clip') {
    prop[prop<trim[1]] = trim[1]
    prop[prop>trim[2]] = trim[2]
  }
  if (normalize=="1") {
    mean_value = mean(data[[form_z]][keep] / prop[keep])
    mean_z     = mean(data[[form_z]][keep])
    prop[keep] = (1 + (1-mean_z)/(mean_value-mean_z) * (1-prop[keep])/prop[keep])^(-1)
  } else if (normalize=="0") {
    mean_value = mean((1-data[[form_z]][keep]) / (1-prop[keep]))
    mean_z     = mean(data[[form_z]][keep])
    prop[keep] = 1 - (1 + mean_z / (mean_value + mean_z - 1) * prop[keep] / (1-prop[keep]))^(-1)
  }
  return(list(prop=prop, keep=keep))
}



# Sub-function to estimate ZSB (2019) for E[Y(1)]
#   L and U: MSM bounds on propensities with L \leq U for every observation 
optimCutHajek = function(Y, Z, Li, Ui, run_min = FALSE) {
  if (run_min) {
    return(-1 * optimCutHajek(-1 * Y, Z, Li, Ui))
  }
  
  #TODO: Add better error code 
  if (min(Li - Ui) > 0 | min(Li) <= 0) {
    stop("Impossible Propensities")
  }
  
  #Subset to treated observations
  Z = as.logical(Z)
  Y  =  Y[Z==TRUE]
  Li = Li[Z==TRUE]
  Ui = Ui[Z==TRUE]
  rm(Z)
  
  #Order observations
  new_order = order(Y)
  Y  =  Y[new_order]
  Li = Li[new_order]
  Ui = Ui[new_order]
  rm(new_order)
  
  #Will not upweight below naive mean 
  start_val = mean(Y / Li) / mean(1 / Li)
  curr_ind  = max(which.max(Y > start_val) - 1, 1)
  curr_nmr  = sum(Y[1:curr_ind] / Ui[1:curr_ind]) + sum(Y[(1+curr_ind):length(Y)] / Li[(1+curr_ind):length(Y)])
  curr_dnm  = sum(     1        / Ui[1:curr_ind]) + sum(             1            / Li[(1+curr_ind):length(Y)])
  
  #Increase cutoff index until optimum is reached
  old_opt = -Inf
  while(curr_nmr / curr_dnm > old_opt) {
    old_opt  = curr_nmr / curr_dnm
    curr_ind = curr_ind + 1
    curr_nmr = curr_nmr + (1 / Ui[curr_ind] - 1 / Li[curr_ind]) * Y[curr_ind]
    curr_dnm = curr_dnm + (1 / Ui[curr_ind] - 1 / Li[curr_ind])
  }
  
  return(old_opt)
}



#===========================
# Regression Functions
#===========================

constregn_option = list()
constregn = function(data, trainmask, testmask, form_x, form_resp, option) {
  rep(mean(if(form_resp%in%colnames(data)) data[[form_resp]] else model.matrix(as.formula(paste('~',form_resp,'-1')), data=data[trainmask,])[,2]), sum(testmask))
}

boostregn_option_bin = list(distribution = 'bernoulli', bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=5, verbose = FALSE)
boostregn_option_cts = list(distribution = 'gaussian', bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=5, verbose = FALSE)
boostregn = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(paste(form_resp, "~", form_x));
  #tryCatch({
  #Sink to suppress gbm's undesirable "CV:" cat expressions 
  sink("temp")
  fit = do.call(gbm, append(list(formula=form, data=data[trainmask,]), option));
  best = if('cv.folds' %in% names(option) && option[['cv.folds']]>0) gbm.perf(fit,plot.it=FALSE,method="cv") else gbm.perf(fit,plot.it=FALSE,method="OOB");
  sink()
  return(predict(fit, n.trees=best, newdata=data[testmask,],  type="response"))
  #}, error = function(err) {warning('failed to fit conditional model; reverting to marginal model'); constregn(data, trainmask, testmask, form_x, form_resp, constregn_option)})
}

svmregn_option_bin = list(kernel="radial", probability=TRUE)
svmregn_option_cts = list(kernel="radial", probability=FALSE, type = "eps-regression")
svmregn = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(paste(ifelse(option$probability, "factor(", ""), form_resp, ifelse(option$probability, ")", ""), "~", form_x));
  fit = do.call(svm, append(list(formula=form, data=data[trainmask,]), option))
  preds = predict(fit, newdata=data[testmask,], probability=option$probability)
  if (option$probability) {
    if (length(unique(data[, form_resp])) > 2) {
      error("Non-binary treatment")
    }
    z_values = unique(data[, form_resp])
    z_treat  = as.character(z_values[as.logical(z_values) == TRUE][1])
    return(attr(preds, "probabilities")[, z_treat])
  } else{
    return(preds)
  }
}

forestregn_option_bin = list(nodesize=1, ntree=1000, na.action=na.omit, replace=TRUE, binary=T)
forestregn_option_cts = list(nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE, binary=F)
forestregn = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(if(option$binary) paste("as.factor(",form_resp, ") ~", form_x) else paste(form_resp, "~", form_x));
  #tryCatch({
  fit = do.call(randomForest, append(list(formula=form, data=data[trainmask,]), option))
  return(if(option$binary) predict(fit, newdata=data[testmask,],  type="prob")[,2] else predict(fit, newdata=data[testmask,],  type="response"))
  #}, error = function(err) {warning('failed to fit conditional model; reverting to marginal model'); constregn(data, trainmask, testmask, form_x, form_resp, constregn_option)})
}

linregn_option_bin = list(family = "binomial")
linregn_option_cts = list(family = "gaussian")
linregn = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(paste(form_resp, "~", form_x));
  #tryCatch({
  fit = do.call(glm, append(list(formula=form, data=data[trainmask,]), option))
  return(predict(fit, newdata=data[testmask,],  type="response"))
  #}, error = function(err) {warning('failed to fit conditional model; reverting to marginal model'); constregn(data, trainmask, testmask, form_x, form_resp, constregn_option)})
}

constquant_option = list()
constquant = function(data, trainmask, testmask, form_x, form_resp, tau, option) {
  rep(quantile(if(form_resp%in%colnames(data)) data[[form_resp]] else model.matrix(as.formula(paste('~',form_resp,'-1')), data=data[trainmask,])[,2], tau), sum(testmask))
}

linquant_option = list()
linquant = function(data, trainmask, testmask, form_x, form_resp, tau, option) {
  form = as.formula(paste(form_resp, "~", form_x));
  tryCatch({
    fit = do.call(rq, append(list(formula=form, data=data[trainmask,], tau=tau), option))
    return(predict(fit, newdata=data[testmask,],  type="response"))
  }, error = function(err) {warning('failed to fit conditional model; reverting to marginal model'); constquant(data, trainmask, testmask, form_x, form_resp, tau, constquant_option)})
}

# Can pass a quantile forest object as option$pretrained_forest 
forestquant_option = list()
forestquant = function(data, trainmask, testmask, form_x, form_resp, tau, option) {
  form = as.formula(paste(form_resp, "~", form_x))
  train_df = model.matrix(lm(paste(form_resp, "~ 0 + ", form_x), data = data)) 
  if (is.null(option$pretrained_forest)) {
    tryCatch({
        fit = do.call(quantile_forest, append(list(formula=form, data=train_df[trainmask,], quantiles=tau), option))
        return(predict(fit, newdata=train_df[testmask,], quantiles=tau)[,1])
      }, error = function(err) {warning('failed to fit forest model; reverting to linear model'); linquant(data, trainmask, testmask, form_x, form_resp, tau, linquant_option)})
  } else {
    tryCatch({
        return(predict(option$pretrained_forest, newdata = train_df[testmask,], quantiles=tau)[,1])
      }, error = function(err) {warning('failed to fit pre-trained forest model; reverting to re-trained forest model'); forestquant(data, trainmask, testmask, form_x, form_resp, tau, forestquant_option)})
  }
}


forestconddist_option = list(nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
forestconddist = function(data, trainmask, testmask, form_x, form_y, option) {
  form = as.formula(paste(form_y, "~", form_x));
  rf = do.call(randomForest, append(list(formula=form, data=data[trainmask,]), option))
  trainids = attr(predict(rf,data[trainmask,],nodes=T),"nodes")
  testids  = attr(predict(rf,data[testmask,],nodes=T),"nodes")
  w = matrix(0L, sum(testmask), sum(trainmask))
  for(i in 1:sum(testmask)) {
    P = t(trainids)==testids[i,]
    w[i,] = colSums(P/rowSums(P))/(dim(P)[1])
  }
  return(w)
}




kernelconddist_option = list(h=.1)
kernelconddist = function(data, trainmask, testmask, form_x, form_y, option) {
  X = model.matrix(as.formula(paste('~',form_x)), data=data)
  X = X[,2:ncol(X)] # remove intercept
  dec = chol( cov(X) )
  tmp = forwardsolve(t(dec), t(X) )
  d2  = as.matrix(dist(t(tmp)))^2
  diag(d2) = Inf
  return(1.*(d2[testmask, trainmask] <= option$h))
}

#Wrapper to use pre-calculated propensities (for simulations, not for production code)
getprop_option = list(col_name = 'nameMissing') 
getprop = function(data, trainmask, testmask, form_x, form_y, option) {
  data$realind = 1:nrow(data)
  data_merge = merge(data, orig_propensities, by = "realind", all.x = TRUE)
  get_col    = gsub("prop_", "prop_Norefit", option$col_name)
  return(data_merge[testmask, c(get_col)])
}

#Functions for binary outcomes (tau = Lambda / (Lambda + 1) for t = 1 and tau = 1/(Lambda + 1) for t = 0)
binaryExtrapolation_option_linear = list(MuhatFunction =   linregn, options =   linregn_option_bin)
binaryExtrapolation_option_boost  = list(MuhatFunction = boostregn, options = boostregn_option_bin)
# Recall that we are only estimating E[Y(1)]; E[Y(0)] happens through subtraction 
# Option: List(muhat = muhat)
binaryExtrapolationquant = function(data, trainmask, testmask, form_x, form_resp, tau, option) {
  as.numeric(option$muhat >= 1-tau)
  #ifelse(rep(tau, nrow(data)) >= 0.5, as.numeric(option$muhat >= 1-tau), as.numeric(option$muhat >= tau))
}

#option: muhat = option$muhat, tau = tau
binaryExtrapolationKappa = function(data, trainmask, testmask, form_x, form_resp, option) {
  tau   = option$tau 
  muhat = option$muhat
  Qhat = binaryExtrapolationquant(data, trainmask, testmask, form_x, form_resp, tau, list(muhat = muhat))
  Lambda = max(tau/(1-tau), (1-tau)/tau)
  if (tau >= 0.5) {
    return(ifelse(Qhat == 1, 1-Lambda^(-1)*(1-muhat), Lambda^( 1)*muhat))
  } else {
    return(ifelse(Qhat == 1, 1-Lambda^( 1)*(1-muhat), Lambda^(-1)*muhat))
  }
}


#===========================
# Estimation Functions
#===========================

## Helper function to bootstrap propensities
bootstrapProps = function(data, form_x, form_z, K, method_prop, option_prop, trim, trim.type, normalize, prop, B, reuse_boot, return_seed, refit_propensities) {
  
  start_boot_seed = .Random.seed 
  
  boot_data = mclapply(1:B, function(b) {
    set.seed(b)
    if (reuse_boot & identical(method_prop, getprop)) {
      prop_col = gsub("prop_", paste0("prop_", ifelse(refit_propensities, "Yes", "No"), "refit"), option_prop$col_name)
      to_add = boot_predata[[b]][, c("realind", "bootind", prop_col, "cv")]
      names(to_add)  = c("realind", "bootind", "prop", "cv")
    } else {
      cvboot  = make.cvgroup.balanced(data, K, form_z)
      realind = sample(1:nrow(data), nrow(data), replace = TRUE)
      if (refit_propensities == TRUE) {
        cvprop = cross.fit.propensities(data[realind,], cvboot[realind], form_x, form_z, method_prop, option_prop, trim=trim, trim.type=trim.type, normalize=normalize)$prop
      } else {
        cvprop = prop$prop[realind]
      }
      to_add = merge(data.frame(realind, bootind = 1:nrow(data), prop = cvprop), data.frame(realind = 1:nrow(data), cv = cvboot), by = "realind", all.x = TRUE)
    } 
    return(to_add)
  })
  if (return_seed) {
    .Random.seed = start_boot_seed 
  }
  return(boot_data)
}

## Helper function to summarize bootstraps
#   No bootstrap estimates ==> everything is returned as NA
summarizeBoots = function(boot_estim) {
  
  #Percentile bootstrap bounds (NA for non-bootstrapped data)
  ci_mean = list()
  ci_50   = list()
  ci_90   = list()
  ci_95   = list()
  ci_sd   = list()
  for (t in 0:1) {
    for (z in 0:1) {
      ci_mean[[paste(z, t)]] = mean(as.numeric(boot_estim[[paste(z, t)]]))
      ci_50[[  paste(z, t)]] = quantile(       boot_estim[[paste(z, t)]], ifelse(t == 1, 1-1.00/2, 1.00/2))
      ci_90[[  paste(z, t)]] = quantile(       boot_estim[[paste(z, t)]], ifelse(t == 1, 1-0.10/2, 0.10/2))
      ci_95[[  paste(z, t)]] = quantile(       boot_estim[[paste(z, t)]], ifelse(t == 1, 1-0.05/2, 0.05/2))
      ci_sd[[  paste(z, t)]] =       sd(       boot_estim[[paste(z, t)]])
      
      #ATT & ATC 
      ci_mean[[paste('AT', z, t)]] = mean(as.numeric((boot_estim[[paste(z, 'Ymean')]] - boot_estim[[paste(1-z, 1-t)]]) / boot_estim[[paste(z, 'Zmean')]]))
      ci_50[[  paste('AT', z, t)]] = quantile((       boot_estim[[paste(z, 'Ymean')]] - boot_estim[[paste(1-z, 1-t)]]) / boot_estim[[paste(z, 'Zmean')]], ifelse(t == 1, 1-1.00/2, 1.00/2))
      ci_90[[  paste('AT', z, t)]] = quantile((       boot_estim[[paste(z, 'Ymean')]] - boot_estim[[paste(1-z, 1-t)]]) / boot_estim[[paste(z, 'Zmean')]], ifelse(t == 1, 1-0.10/2, 0.10/2))
      ci_95[[  paste('AT', z, t)]] = quantile((       boot_estim[[paste(z, 'Ymean')]] - boot_estim[[paste(1-z, 1-t)]]) / boot_estim[[paste(z, 'Zmean')]], ifelse(t == 1, 1-0.05/2, 0.05/2))
      ci_sd[[  paste('AT', z, t)]] =       sd((       boot_estim[[paste(z, 'Ymean')]] - boot_estim[[paste(1-z, 1-t)]]) / boot_estim[[paste(z, 'Zmean')]])
    }
  }
  
  #ATE
  ci_mean[[paste('1 1 0 0')]] = mean(as.numeric(boot_estim[['1 1']] - boot_estim[['0 0']]))
  ci_50[[  paste('1 1 0 0')]] = quantile(       boot_estim[['1 1']] - boot_estim[['0 0']], 1-1.00/2)
  ci_90[[  paste('1 1 0 0')]] = quantile(       boot_estim[['1 1']] - boot_estim[['0 0']], 1-0.10/2)
  ci_95[[  paste('1 1 0 0')]] = quantile(       boot_estim[['1 1']] - boot_estim[['0 0']], 1-0.05/2)
  ci_sd[[  paste('1 1 0 0')]] =       sd(       boot_estim[['1 1']] - boot_estim[['0 0']])
  
  ci_mean[[paste('1 0 0 1')]] = mean(as.numeric(boot_estim[['1 0']] - boot_estim[['0 1']]))
  ci_50[[  paste('1 0 0 1')]] = quantile(       boot_estim[['1 0']] - boot_estim[['0 1']], 1.00/2)
  ci_90[[  paste('1 0 0 1')]] = quantile(       boot_estim[['1 0']] - boot_estim[['0 1']], 0.10/2)
  ci_95[[  paste('1 0 0 1')]] = quantile(       boot_estim[['1 0']] - boot_estim[['0 1']], 0.05/2)
  ci_sd[[  paste('1 0 0 1')]] =       sd(       boot_estim[['1 0']] - boot_estim[['0 1']])
  
  return(list(ci_mean, ci_50, ci_90, ci_95, ci_sd))
}



## Helper function to turn estimates into a summary table 
summarizeResults = function(Lambda, data, form_y, form_z, iif, se_list, boot_facts) {
  
  ci_mean = boot_facts[[1]]
  ci_50   = boot_facts[[2]]
  ci_90   = boot_facts[[3]]
  ci_95   = boot_facts[[4]]
  ci_sd   = boot_facts[[5]]
  
  ATT_plus  = mean(data[[form_y]] - iif[['0 0']]) / mean(data[[form_z]])
  ATT_minus = mean(data[[form_y]] - iif[['0 1']]) / mean(data[[form_z]])
  
  ATC_plus  = mean(data[[form_y]] - iif[['1 0']]) / mean(1-data[[form_z]])
  ATC_minus = mean(data[[form_y]] - iif[['1 1']]) / mean(1-data[[form_z]])
  
  
  data.frame(
    Lambda = Lambda,
    estimand = c('1','1','0','0','ATE','ATE', 'ATT', 'ATT', 'ATC', 'ATC'),
    upper = c(T,F,T,F,T,F,T,F,T,F),
    estimate = c(
      mean(iif[['1 1']]),
      mean(iif[['1 0']]),
      mean(iif[['0 1']]),
      mean(iif[['0 0']]),
      mean(iif[['1 1']]-iif[['0 0']]),
      mean(iif[['1 0']]-iif[['0 1']]),
      ATT_plus,
      ATT_minus,
      ATC_plus,
      ATC_minus
    ),
    sterr_iif = c(
      se_list[['1 1']],
      se_list[['1 0']],
      se_list[['0 1']],
      se_list[['0 0']],
      se_list[['1 1 0 0']],
      se_list[['1 0 0 1']],
      sd((data[[form_y]] - iif[['0 0']] - data[[form_z]] * ATT_plus )/mean(data[[form_z]])) / sqrt(nrow(data)), #(Y - psi_i - Z*psihat)/Zbar where psi_i = IF for E[Y(0)] and psihat is ATT estimate
      sd((data[[form_y]] - iif[['0 1']] - data[[form_z]] * ATT_minus)/mean(data[[form_z]])) / sqrt(nrow(data)), 
      sd((data[[form_y]] - iif[['1 0']] - (1-data[[form_z]]) * ATC_plus )/mean(1-data[[form_z]])) / sqrt(nrow(data)),
      sd((data[[form_y]] - iif[['1 1']] - (1-data[[form_z]]) * ATC_minus)/mean(1-data[[form_z]])) / sqrt(nrow(data))
    ),
    sterr_boot = c(
      ci_sd[['1 1']],
      ci_sd[['1 0']],
      ci_sd[['0 1']],
      ci_sd[['0 0']],
      ci_sd[['1 1 0 0']],
      ci_sd[['1 0 0 1']],
      ci_sd[['AT 1 1']],
      ci_sd[['AT 1 0']],
      ci_sd[['AT 0 1']],
      ci_sd[['AT 0 0']]
    ),
    bootmean = c(
      ci_mean[['1 1']],
      ci_mean[['1 0']],
      ci_mean[['0 1']],
      ci_mean[['0 0']],
      ci_mean[['1 1 0 0']],
      ci_mean[['1 0 0 1']],
      ci_mean[['AT 1 1']],
      ci_mean[['AT 1 0']],
      ci_mean[['AT 0 1']],
      ci_mean[['AT 0 0']]
    ),
    quantiles50 = c(
      ci_50[['1 1']],
      ci_50[['1 0']],
      ci_50[['0 1']],
      ci_50[['0 0']],
      ci_50[['1 1 0 0']],
      ci_50[['1 0 0 1']],
      ci_50[['AT 1 1']],
      ci_50[['AT 1 0']],
      ci_50[['AT 0 1']],
      ci_50[['AT 0 0']]
    ),
    quantiles90 = c(
      ci_90[['1 1']],
      ci_90[['1 0']],
      ci_90[['0 1']],
      ci_90[['0 0']],
      ci_90[['1 1 0 0']],
      ci_90[['1 0 0 1']],
      ci_90[['AT 1 1']],
      ci_90[['AT 1 0']],
      ci_90[['AT 0 1']],
      ci_90[['AT 0 0']]
    ),
    quantiles95 = c(
      ci_95[['1 1']],
      ci_95[['1 0']],
      ci_95[['0 1']],
      ci_95[['0 0']],
      ci_95[['1 1 0 0']],
      ci_95[['1 0 0 1']],
      ci_95[['AT 1 1']],
      ci_95[['AT 1 0']],
      ci_95[['AT 0 1']],
      ci_95[['AT 0 0']]
    )
  )
}


## Lambdas: list of MSM parameters to consider
## method_prop, option_prop: binary regression method for learning propensities
## method_quant, option_quant: quantile regression method for learning q
## method_regn, option_regn: regression method for learning kappa
## method_conddist, option_conddist: method for learning conditional distribution
## conddist_quant: TRUE means use method_conddist to fit quantiles, FALSE means use method_quant
## conddist_kappa: TRUE means use method_conddist to fit kappa, FALSE means use method_regn
## semiadaptive: TRUE means use the same out-of-fold data to fit q and kappa, FALSE means split the out-of-fold into two (requires K>=5 odd)
## bootstrap_settings: return_seed means return the seed at the end of the bootstrap estimates to the seed at start for comparability 
dvds = function(Lambdas, data, form_x, form_z, form_y, method_prop, option_prop, method_quant=NULL, option_quant=NULL, method_regn=NULL, option_regn=NULL, method_conddist=NULL, option_conddist=NULL, conddist_quant=F, conddist_kappa=F, K=5, semiadaptive=FALSE, trim=c(0.01,0.99), trim.type='clip', normalize="F", form_x_quant = NULL, form_x_kappa = NULL, boot_infer = F, boot_settings = list(reuse_boot = T, refit_propensities = F, B = 500, return_seed = TRUE)) {
  data = data%>%arrange(!! sym(form_y))
  if (identical(method_prop, getprop)) {
    cvgroup = cvgroup_orig
  } else {
    cvgroup = make.cvgroup.balanced(data, K, form_z)
  }
  if (is.null(form_x_quant)) {
    form_x_quant = form_x
  }
  if (is.null(form_x_kappa)) {
    form_x_kappa = paste0('Q + ', form_x)
  }
  prop = cross.fit.propensities(data, cvgroup, form_x, form_z, method_prop, option_prop, trim=trim, trim.type=trim.type, normalize=normalize)
  data = data[prop$keep,]
  e = prop$prop[prop$keep]
  
  #For binary outcomes: Estimate kappa with option_regn$MuhatFunction, with options from option_regn$options
  if (identical(method_quant, binaryExtrapolationquant) | identical(method_regn, binaryExtrapolationKappa)) {
    muhat_list = list()
    for (z in 0:1) {
      muhat_list[[as.character(z)]] = rep(as.numeric(NA), nrow(data))
      for (k in unique(cvgroup)) {
        muhat_list[[as.character(z)]][cvgroup==k] = option_quant$MuhatFunction(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_z]]==z, cvgroup==k, form_x_quant, form_y, option_quant$options)
      }
    }
  }
  #Initialize list of standard errors
  se_list = list()
  if (boot_infer) {
    boot_data = bootstrapProps(data, form_x, form_z, K, method_prop, option_prop, trim, trim.type, normalize, prop, boot_settings$B, boot_settings$reuse_boot, boot_settings$return_seed, boot_settings$refit_propensities)
  }
  #print("Done with bootstrap propensity calculation")
  
  
  if(conddist_quant || conddist_kappa) {
    w = list()
    wc = list()
    for (z in 0:1) {
      w[[paste(z)]] = matrix(0L,nrow=nrow(data),ncol=nrow(data))
      for (k in 1:K) {
        w[[paste(z)]][cvgroup==k, cvgroup!=k & data[[form_z]]==z] = method_conddist(data, cvgroup!=k & data[[form_z]]==z, cvgroup==k, form_x, form_y, option_conddist)
      }
      w[[paste(z)]] = w[[paste(z)]]/apply(w[[paste(z)]], 1, sum)
      wc[[paste(z)]] = t(apply(w[[paste(z)]], 1, cumsum))
    }
  }
  
  #Train quantile forests (if applicable)
  if (identical(method_quant, forestquant) & !is.null(option_quant$reuse_forests)) {
    if (option_quant$reuse_forests == T) {
      pretrained_forests = list()
      all_taus = sort(unique(c(1 / (1 + Lambdas), Lambdas / (1+Lambdas))))
      train_df = model.matrix(lm(paste0(form_y, " ~ 0 + ", form_x_quant), data = data)) 
      for (k in 1:K) {
        for (z in 0:1) {
          #data table to allow functions in form_x 
          pretrained_forests[[paste(k, z)]] = quantile_forest(train_df[cvgroup!=k & data[[form_z]]==z,], data[[form_y]][cvgroup!=k & data[[form_z]]==z], quantiles = all_taus)
        }
      }
      option_quant$pretrained_forests = pretrained_forests 
      rm(train_df, all_taus, k, z)
    }
  }
  
  return(foreach(Lambda=Lambdas, .combine=rbind) %do% {
    eif        = list()
    boot_estim = list()
    
    for (z in 0:1) {
      for (t in 0:1) {
        q = numeric(nrow(data));
        tau = ((1-t)+t*Lambda)/(Lambda + 1)
        if (conddist_quant) {
          q = data[[form_y]][apply(wc[[paste(z)]]>=tau,1,which.max)]
        } else if (identical(method_quant, binaryExtrapolationquant)) {
          option_quant = list(muhat = muhat_list[[as.character(z)]])
          q = method_quant(data, NA, NA, form_x_quant, form_y, tau, option_quant)
        } else {
          for (k in 1:K) {
            if (!is.null(option_quant$pretrained_forests[[paste(k, z)]])) {
              option_quant$pretrained_forest = option_quant$pretrained_forests[[paste(k, z)]]
            }
            q[cvgroup==k] = method_quant(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==0) & data[[form_z]]==z, cvgroup==k, form_x_quant, form_y, tau, option_quant)
          }
        }
        kappa = numeric(nrow(data));
        data[['.hinge']] = data[[form_y]]/Lambda + (1. - 1. / Lambda) * (q + (data[[form_y]]-q)*Lambda^sign((2*t-1)*data[[form_y]] - q))
        if(conddist_kappa) {
          kappa = w[[paste(z)]] %*% data[['.hinge']]
        } else if (identical(method_regn, binaryExtrapolationKappa)) {
          kappa = method_regn(data, NA, NA, NA, NA, list(muhat = muhat_list[[as.character(z)]], tau = tau))
        } else {
          data$Q = q
          for (k in 1:K) {
            kappa[cvgroup==k] = method_regn(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_z]]==z, cvgroup==k, form_x_kappa, '.hinge', option_regn)
          }
          data$Q = NULL
        }
        zz = (2*z-1)*data[[form_z]] + (1 - z)
        ee = (2*z-1)*e + (1 - z)
        eif[[paste(z,t)]] = kappa + zz*(q - kappa)/ee + zz*(data[[form_y]] - q)*(ee + (1-ee)*Lambda^((2*t-1)*sign(data[[form_y]] - q)))/ee
        
        if (boot_infer) {
          start_boot_seed = .Random.seed
          
          boot_estims_wrongdirection = mclapply(1:boot_settings$B, function(b) {
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
          
          for (b in 1:boot_settings$B) {
            boot_estim[[paste(z,t)]][b]       = boot_estims_wrongdirection[[b]][1]
            boot_estim[[paste(z,'Ymean')]][b] = boot_estims_wrongdirection[[b]][2]
            boot_estim[[paste(z,'Zmean')]][b] = boot_estims_wrongdirection[[b]][3]
          }
          
          if (boot_settings$return_seed) {
            .Random.seed = start_boot_seed 
          }
        } 
        se_list[[paste(z,t)]] = sd(eif[[paste(z, t)]])/sqrt(nrow(data))
      }
    }
    
    #ATE standard errors
    se_list[['1 1 0 0']] = sd(eif[['1 1']]-eif[['0 0']])/sqrt(nrow(data))
    se_list[['1 0 0 1']] = sd(eif[['1 0']]-eif[['0 1']])/sqrt(nrow(data))
    
    #No bootstrap estimates ==> returns NA or NaN
    boot_facts = summarizeBoots(boot_estim)
    summarizeResults(Lambda, data, form_y, form_z, eif, se_list, boot_facts)
  })
}


## Lambdas: list of MSM parameters to consider
## method_prop, option_prop: binary regression method for learning propensities
zsb = function(Lambdas, data, form_x, form_z, form_y, method_prop, option_prop, method_regn=NULL, option_regn=NULL, K=5, semiadaptive=FALSE, trim=c(0.01,0.99), trim.type='clip', normalize="F", boot_infer = F, boot_settings = list(reuse_boot = T, refit_propensities = F, B = 500, return_seed = TRUE)) {
  #Now main processing
  data = data%>%arrange(!! sym(form_y))
  if (identical(method_prop, getprop)) {
    cvgroup = cvgroup_orig
  } else {
    cvgroup = make.cvgroup.balanced(data, K, form_z)
  }
  prop = cross.fit.propensities(data, cvgroup, form_x, form_z, method_prop, option_prop, trim=trim, trim.type=trim.type, normalize=normalize)
  data  = data[     prop$keep,]
  e     = prop$prop[prop$keep]
  cvgroup = cvgroup[prop$keep]
  #Regression centering (zero for IPW)
  muhat = list("1" = rep(0, nrow(data)), "0" = rep(0, nrow(data)))
  if (!is.null(method_regn)) {
    for (z in 0:1) {
      for (k in unique(cvgroup)) {
        muhat[[as.character(z)]][cvgroup == k] = method_regn(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_z]]==z, cvgroup==k, form_x, form_y, option_regn)
      }
    }
    rm(z, k)
  }
  #Initialize list of standard errors
  se_list = list()
  if (boot_infer) {
    boot_data = bootstrapProps(data, form_x, form_z, K, method_prop, option_prop, trim, trim.type, normalize, prop, boot_settings$B, boot_settings$reuse_boot, boot_settings$return_seed, boot_settings$refit_propensities)
  }
  #print("Done with bootstrap propensity calculation")
  return(foreach(Lambda=Lambdas, .combine=rbind) %do% {
    
    iif        = list()
    boot_estim = list()
    for (z in 0:1) {
      data_muhat = muhat[[as.character(z)]]
      data_resid = data[[form_y]] - data_muhat
      zz = (2*z-1)*data[[form_z]] + (1 - z)
      ee = (2*z-1)*e + (1 - z)
      #Inverse of lower propensity (higher weight)
      L  = ee / (ee + (1-ee)*Lambda^(+1))
      U  = ee / (ee + (1-ee)*Lambda^(-1))
      for (t in 0:1) {
        zsb_thresh  = optimCutHajek(data_resid, zz, L, U, run_min = as.logical(1-t))
        opt_invprop = (ee + (1-ee)*Lambda^((2*t-1)*sign(data_resid - zsb_thresh)))/ee
        iif[[paste(z,t)]] = data_muhat + (zz * data_resid * opt_invprop) / mean(zz * opt_invprop)
        
        if (boot_infer) {
          start_boot_seed = .Random.seed
          
          boot_estims_wrongdirection = mclapply(1:boot_settings$B, function(b) {
            set.seed(b)
            ind_star       = boot_data[[b]]$realind
            data_star      = data[ind_star, ]
            data_star$prop = boot_data[[b]]$prop
            e_star   = data_star$prop
            Y_star   = data_star[[form_y]]
            muhat_star = muhat[[as.character(z)]][ind_star]
            resid_star = Y_star - muhat_star
            zz_star  = (2*z-1)*data_star[[form_z]] + (1 - z)
            ee_star  = (2*z-1)*e_star + (1 - z)
            L_star  = ee_star / (ee_star + (1-ee_star)*Lambda^(+1))
            U_star  = ee_star / (ee_star + (1-ee_star)*Lambda^(-1))
            zsb_thresh_star  = optimCutHajek(resid_star, zz_star, L_star, U_star, run_min = as.logical(1-t))
            opt_invprop_star = (ee_star + (1-ee_star)*Lambda^((2*t-1)*sign(resid_star - zsb_thresh_star)))/ee_star
            
            c(
              mean(muhat_star + (zz_star * resid_star * opt_invprop_star) / mean(zz_star * opt_invprop_star)),
              mean(Y_star),
              mean(zz_star)
            )
          })
          
          for (b in 1:boot_settings$B) {
            boot_estim[[paste(z,t)]][b]       = boot_estims_wrongdirection[[b]][1]
            boot_estim[[paste(z,'Ymean')]][b] = boot_estims_wrongdirection[[b]][2]
            boot_estim[[paste(z,'Zmean')]][b] = boot_estims_wrongdirection[[b]][3]
          }
          
          if (boot_settings$return_seed) {
            .Random.seed = start_boot_seed 
          }
        } 
        se_list[[paste(z,t)]] = sd(iif[[paste(z, t)]])/sqrt(nrow(data))
        
      }
    }
    #ATE standard errors
    se_list[['1 1 0 0']] = sd(iif[['1 1']]-iif[['0 0']])/sqrt(nrow(data))
    se_list[['1 0 0 1']] = sd(iif[['1 0']]-iif[['0 1']])/sqrt(nrow(data))
    
    
    #No bootstrap estimates ==> returns NA or NaN
    boot_facts = summarizeBoots(boot_estim)
    
    summarizeResults(Lambda, data, form_y, form_z, iif, se_list, boot_facts)
  })
}

## Lambdas: list of MSM parameters to consider
## method_prop, option_prop: binary regression method for learning propensities
## method_quant, option_quant: quantile regression method for learning q
## semiadaptive: TRUE means use the same out-of-fold data to fit q, FALSE means split the out-of-fold into two (requires K>=5 odd)
qb = function(Lambdas, data, form_x, form_z, form_y, method_prop, option_prop, method_quant=NULL, option_quant=NULL, K=5, semiadaptive=FALSE, trim=c(0.01,0.99), trim.type='clip', normalize="F", form_x_quant = NULL, boot_infer = F, boot_settings = list(reuse_boot = T, refit_propensities = F, B = 500, return_seed = TRUE)) {
  
  if (is.null(form_x_quant)) {
    form_x_quant = form_x
  }
  data = data%>%arrange(!! sym(form_y))
  if (identical(method_prop, getprop)) {
    cvgroup = cvgroup_orig
  } else {
    cvgroup = make.cvgroup.balanced(data, K, form_z)
  }
  prop = cross.fit.propensities(data, cvgroup, form_x, form_z, method_prop, option_prop, trim=trim, trim.type=trim.type, normalize=normalize)
  data = data[prop$keep,]
  e = prop$prop[prop$keep]
  #For binary outcomes: can use muhat for qhat 
  if (identical(method_quant, binaryExtrapolationquant)) {
    muhat_list = list()
    for (z in 0:1) {
      muhat_list[[as.character(z)]] = rep(as.numeric(NA), nrow(data))
      for (k in unique(cvgroup)) {
        muhat_list[[as.character(z)]][cvgroup==k] = option_quant$MuhatFunction(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_z]]==z, cvgroup==k, form_x_quant, form_y, option_quant$options)
      }
    }
  }
  #Initialize list of standard errors
  se_list = list()
  boot_estim = list()
  if (boot_infer) {
    boot_data  = bootstrapProps(data, form_x, form_z, K, method_prop, option_prop, trim, trim.type, normalize, prop, boot_settings$B, boot_settings$reuse_boot, boot_settings$return_seed, boot_settings$refit_propensities)
  }
  
  
  
  #Train quantile forests (if applicable)
  if (identical(method_quant, forestquant) & !is.null(option_quant$reuse_forests)) {
    if (option_quant$reuse_forests == T) {
      pretrained_forests = list()
      all_taus = sort(unique(c(1 / (1 + Lambdas), Lambdas / (1+Lambdas))))
      train_df = model.matrix(lm(paste0(form_y, " ~ 0 + ", form_x_quant), data = data)) 
      for (k in 1:K) {
        for (z in 0:1) {
          #data table to allow functions in form_x 
          pretrained_forests[[paste(k, z)]] = quantile_forest(train_df[cvgroup!=k & data[[form_z]]==z,], data[[form_y]][cvgroup!=k & data[[form_z]]==z], quantiles = all_taus)
        }
      }
      option_quant$pretrained_forests = pretrained_forests 
      rm(train_df, all_taus, k, z)
    }
  }
  return(foreach(Lambda=Lambdas, .combine=rbind) %do% {
    iif = list()
    for (z in 0:1) {
      zz = (2*z-1)*data[[form_z]] + (1 - z)
      ee = (2*z-1)*e + (1 - z)
      for (t in 0:1) {
        #First-stage quantiles 
        q = numeric(nrow(data));
        tau = ((1-t)+t*Lambda)/(Lambda + 1)
        if (identical(method_quant, binaryExtrapolationquant)) {
          option_quant = list(muhat = muhat_list[[as.character(z)]])
          q = method_quant(data, NA, NA, form_x_quant, form_y, tau, option_quant)
        } else {
          for (k in 1:K) {
            if (!is.null(option_quant$pretrained_forests[[paste(k, z)]])) {
              option_quant$pretrained_forest = option_quant$pretrained_forests[[paste(k, z)]]
            }
            q[cvgroup==k] = method_quant(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==0) & data[[form_z]]==z, cvgroup==k, form_x_quant, form_y, tau, option_quant)
          }
        }
        #Get quantiles for quantile balancing 
        quant_df = data.frame(Y = data[[form_y]], Q = q)
        if (sd(q) == 0) {
          #warning("Constant quantile estimates. Quantile balancing will only impose one constraint.")
          q_balance = predict(rq(Y ~ 1, tau=tau, weights=((1-ee)/ee)[zz==1], data = quant_df[zz==1,]), quant_df)
        } else {
          q_balance = predict(rq(Y ~ Q, tau=tau, weights=((1-ee)/ee)[zz==1], data = quant_df[zz==1,]), quant_df)
        }
        
        #Return i's "influence function," the quantile balancing term 
        # (TODO: does this produce valid standard errors?)
        iif[[paste(z,t)]] = (zz*q_balance/ee + zz*(data[[form_y]] - q_balance)*(ee + (1-ee)*Lambda^((2*t-1)*sign(data[[form_y]] - q_balance)))/ee) / mean(zz/ee)
        
        if (boot_infer) {
          start_boot_seed = .Random.seed
          
          boot_estims_wrongdirection = mclapply(1:boot_settings$B, function(b) {
            set.seed(b)
            ind_star       = boot_data[[b]]$realind
            data_star      = data[ind_star, ]
            data_star$prop = boot_data[[b]]$prop
            q_star   = q[ind_star]
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
          
          for (b in 1:boot_settings$B) {
            boot_estim[[paste(z,t)]][b]       = boot_estims_wrongdirection[[b]][1]
            boot_estim[[paste(z,'Ymean')]][b] = boot_estims_wrongdirection[[b]][2]
            boot_estim[[paste(z,'Zmean')]][b] = boot_estims_wrongdirection[[b]][3]
          }
          
          if (boot_settings$return_seed) {
            .Random.seed = start_boot_seed 
          }
        } 
        se_list[[paste(z,t)]] = sd(iif[[paste(z, t)]])/sqrt(nrow(data))
      }
    }
    #ATE standard errors
    se_list[['1 1 0 0']] = sd(iif[['1 1']]-iif[['0 0']])/sqrt(nrow(data))
    se_list[['1 0 0 1']] = sd(iif[['1 0']]-iif[['0 1']])/sqrt(nrow(data))
    
    boot_facts = summarizeBoots(boot_estim)
    summarizeResults(Lambda, data, form_y, form_z, iif, se_list, boot_facts)
  })
}


