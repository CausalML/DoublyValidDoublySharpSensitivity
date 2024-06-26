
library(dplyr)
library(progress)
library(ggplot2)
library(grf)

# ==================================================================
# Helper functions

phi.continuous = function(X, Y, Z, mu, e, Lambdas) {
  # DVDS generic influence function for continuous outcomes
  n = length(Y)
  m = length(Lambdas)
  taus = Lambdas/(Lambdas+1)
  
  # Estimate quantiles using quantile forests
  Q = matrix(numeric(m*n), nrow = n)
  qfit = quantile_forest(X[Z==1,], Y[Z==1], quantiles = taus)
  Q[Z==1,] = predict(qfit)$predictions; Q[Z==0,] = predict(qfit, X[Z==0,])$predictions
  
  # Estimate CV@R using random forests.  Uses the GRF weights from the median tau for all taus
  cvar = matrix(numeric(m*n), nrow = n)
  k = median(1:m)
  check = sapply(1:m, function(i) { Q[,i] + pmax(0, Y-Q[,i])/(1-taus[i]) })
  cvarfit = regression_forest(X[Z==1,], check[Z==1,k])
  cvar[Z==1,] = as.matrix(get_forest_weights(cvarfit) %*% check[Z==1,])
  cvar[Z==0,] = as.matrix(get_forest_weights(cvarfit, X[Z==0,]) %*% check[Z==1,])
  
  # Form the influence function
  eif = sapply(1:m, function(i) {
    kappa = (1/Lambdas[i])*mu + (1-1/Lambdas[i])*cvar[,i]
    pseudo= (1/Lambdas[i])*Y  + (1-1/Lambdas[i])*check[,i]
    Z*Y + (1-Z)*kappa + (1-e)/e*Z*(pseudo - kappa)
  })
  return(eif)
}


# ==================================================================
# Main DVDS functions

# X is only used for continuous, but included here to analogize 
dvds.binary = function(Y, Z, mu0, mu1, e, Lambda=1, alpha=0.05) {
  # DVDS sensitivity analysis for binary outcomes
  tau = Lambda/(Lambda+1)
  Q1plus = (mu1 > 1-tau)
  Q0plus = (mu0 > 1-tau)
  Q1mins = (mu1 > tau)
  Q0mins = (mu0 > tau)
  kappa1plus = pmin(1 - 1/Lambda + mu1/Lambda, mu1*Lambda)
  kappa0plus = pmin(1 - 1/Lambda + mu0/Lambda, mu0*Lambda)
  kappa1mins = pmax(1 - Lambda + mu1*Lambda, mu1/Lambda)
  kappa0mins = pmax(1 - Lambda + mu0*Lambda, mu0/Lambda)
  ATEplus1 = (Z*Y + (1-Z)*kappa1plus + (1-e)/e*Z*(Q1plus + Lambda^sign(Y-Q1plus)*(Y-Q1plus) - kappa1plus))
  ATEplus0 = ((1-Z)*Y + Z*kappa0mins + e/(1-e)*(1-Z)*(Q0mins + Lambda^sign(Q0mins-Y)*(Y-Q0mins) - kappa0mins))
  ATEmins1 = (Z*Y + (1-Z)*kappa1mins + (1-e)/e*Z*(Q1mins + Lambda^sign(Q1mins-Y)*(Y-Q1mins) - kappa1mins))
  ATEmins0 = ((1-Z)*Y + Z*kappa0plus + e/(1-e)*(1-Z)*(Q0plus + Lambda^sign(Y-Q0plus)*(Y-Q0plus) - kappa0plus))
  ATEplus = 
    (Z*Y + (1-Z)*kappa1plus + (1-e)/e*Z*(Q1plus + Lambda^sign(Y-Q1plus)*(Y-Q1plus) - kappa1plus)) -
    ((1-Z)*Y + Z*kappa0mins + e/(1-e)*(1-Z)*(Q0mins + Lambda^sign(Q0mins-Y)*(Y-Q0mins) - kappa0mins))
  ATEmins = 
    (Z*Y + (1-Z)*kappa1mins + (1-e)/e*Z*(Q1mins + Lambda^sign(Q1mins-Y)*(Y-Q1mins) - kappa1mins)) -
    ((1-Z)*Y + Z*kappa0plus + e/(1-e)*(1-Z)*(Q0plus + Lambda^sign(Y-Q0plus)*(Y-Q0plus) - kappa0plus))
  c = qnorm(1 - alpha/2)/sqrt(length(Y))
  return(data.frame("Lambda" = Lambda,
                    "upper" = mean(ATEplus),
                    "lower" = mean(ATEmins),
                    "upper.CI" = mean(ATEplus) + c*sd(ATEplus),
                    "lower.CI" = mean(ATEmins) - c*sd(ATEmins),
                    "upper.1" = mean(ATEplus1),
                    "lower.1" = mean(ATEmins1),
                    "upper.0" = mean(ATEmins0),
                    "lower.0" = mean(ATEplus0)))
}


dvds.continuous = function(X, Y, Z, mu0, mu1, e, Lambdas=1, alpha=0.05) {
  # DVDS sensitivity analysis for continuous outcomes
  # Uses random forests for all quantile and CV@R nuisances
  n = length(Y)
  phi1plus = phi.continuous(X, Y, Z, mu1, e, Lambdas)
  phi1mins =-phi.continuous(X,-Y, Z,-mu1, e, Lambdas)
  phi0plus = phi.continuous(X, Y, 1-Z, mu0, 1-e, Lambdas)
  phi0mins =-phi.continuous(X,-Y, 1-Z,-mu0, 1-e, Lambdas)
  ATEplus = phi1plus - phi0mins
  ATEmins = phi1mins - phi0plus
  c = qnorm(1 - alpha/2)/sqrt(n)
  lapply(1:length(Lambdas), function(i) {
    data.frame("Lambda" = Lambdas[i],
               "lower" = mean(ATEmins[,i]),
               "upper" = mean(ATEplus[,i]),
               "lower.CI" = mean(ATEmins[,i]) - c*sd(ATEmins[,i]),
               "upper.CI" = mean(ATEplus[,i]) + c*sd(ATEplus[,i]),
               "upper.1" = mean(phi1plus),
               "lower.1" = mean(phi1mins),
               "upper.0" = mean(phi0plus),
               "lower.0" = mean(phi0mins)
    )
  }) %>% bind_rows() %>% return()
}


# ==================================================================
# ZSB comparison
extrema = function(A, Y, gamma, fitted.prob) {
  # Helper function for ZSB sensitivity analysis
  fitted.logit <- qlogis(fitted.prob)
  eg <- exp(-fitted.logit)
  Y <- Y[A == 1]
  eg <- eg[A == 1]
  eg <- eg[order(-Y)]
  Y <- Y[order(-Y)]
  num.each.low <- Y * (1 + exp(-gamma) * eg)
  num.each.up <- Y * (1 + exp(gamma) * eg)
  num <- c(0, cumsum(num.each.up)) + c(rev(cumsum(rev(num.each.low))), 0)
  den.each.low <- (1 + exp(-gamma) * eg)
  den.each.up <- (1 + exp(gamma) * eg)
  den <- c(0, cumsum(den.each.up)) + c(rev(cumsum(rev(den.each.low))), 0)
  maximum <- max(num/den)
  num <- c(0, cumsum(num.each.low)) + c(rev(cumsum(rev(num.each.up))), 0)
  den <- c(0, cumsum(den.each.low)) + c(rev(cumsum(rev(den.each.up))), 0)
  minimum <- min(num/den)
  c(minimum, maximum)
}

extrema.os = function(Z, Y, Lambda, e) {
  extrema(Z, Y, log(Lambda), e) - rev(extrema(1-Z, Y, log(Lambda), 1-e))
}

extrema.aipw = function(Z, Y, Lambda, e, mu0, mu1) {
  # ZSB sensitivity analysis for AIPW
  eps = Y - ifelse(Z==1, mu1, mu0)
  bounds = mean(mu1 - mu0) + extrema.os(Z, eps, Lambda, e)
  return(data.frame("Lambda" = Lambda,
                    "upper" = bounds[2],
                    "lower" = bounds[1],
                    "upper.CI" = NA,
                    "lower.CI" = NA,
                    "upper.1" = NA,
                    "lower.1" = NA,
                    "upper.0" = NA,
                    "lower.0" = NA))
}




