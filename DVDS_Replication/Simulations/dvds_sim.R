library(dplyr)
library(progress)
library(grf)
source("simplified_dvds.R")
set.seed(1)

Nsims = 1000
Lambdas = seq(1, 3, length.out = 20)
pb = progress_bar$new(":eta [:bar] :percent", total=Nsims)
for (i in 1:Nsims) {
  pb$tick() 
  
  # Generate data
  n = 1000
  d = 5
  X = matrix(runif(n*d, min=-1, max=1), nrow=n)
  Z = rbinom(n, 1, 1/(1 + exp(X[,1] + 0.5*(X[,2]>0) + 0.5*X[,2]*X[,3])))
  Ybinary = rbinom(n, 1, 1/(1 + exp(0.5*X[,1] + X[,2] + 0.25*X[,2]*X[,3])))
  Ycontin = 2*sign(X[,1]) + X[,2] + X[,2]*X[,3] + (1 + X[,4]^2)*rnorm(n)
  
  # Estimate shared nuisances
  e = predict(probability_forest(X=X, Y=factor(Z)))$prediction[,"1"]
  mu0binary = mu1binary = mu0contin = mu1contin = numeric(n)
  fit0binary = probability_forest(X = X[Z==0,], Y = factor(Ybinary[Z==0]))
  fit1binary = probability_forest(X = X[Z==1,], Y = factor(Ybinary[Z==1]))
  fit0contin = regression_forest(X = X[Z==0,], Y = Ycontin[Z==0])
  fit1contin = regression_forest(X = X[Z==1,], Y = Ycontin[Z==1])
  mu0binary[Z==0] = predict(fit0binary)$predictions[,"1"]
  mu0binary[Z==1] = predict(fit0binary, X[Z==1,])$predictions[,"1"]
  mu1binary[Z==0] = predict(fit1binary, X[Z==0,])$predictions[,"1"]
  mu1binary[Z==1] = predict(fit1binary)$predictions[,"1"]
  mu0contin[Z==0] = predict(fit0contin)$predictions
  mu0contin[Z==1] = predict(fit0contin, X[Z==1,])$predictions
  mu1contin[Z==0] = predict(fit1contin, X[Z==0,])$predictions
  mu1contin[Z==1] = predict(fit1contin)$predictions

  # Fit bounds
  write.table(
    x = bind_rows(lapply(Lambdas, function(l) {
      dvds.binary(X=X, Y=Ybinary, Z=Z, mu0=mu0binary, mu1=mu1binary, e=e, Lambda=l)
    })),
    file = "dvds_output/binary_dvds.csv",
    append = (i != 1),
    row.names = F,
    col.names = (i==1)
  )
  write.table(
    x = dvds.continuous(X=X, Y=Ycontin, Z=Z, mu0=mu0contin, mu1=mu1contin, e=e, Lambdas=Lambdas),
    file = "dvds_output/continuous_dvds.csv",
    append = (i != 1),
    row.names = F,
    col.names = (i==1)
  )
  write.table(
    x = bind_rows(lapply(Lambdas, function(l) {
      extrema.aipw(Z=Z, X=X, Y=Ybinary, Lambda=l, e=e, mu0=mu0binary, mu1=mu1binary)
    })),
    file = "dvds_output/binary_zsb.csv",
    append = (i != 1),
    row.names = F,
    col.names = (i==1)
  )
  write.table(
    x = bind_rows(lapply(Lambdas, function(l) {
      extrema.aipw(Z=Z, X=X, Y=Ycontin, Lambda=l, e=e, mu0=mu0contin, mu1=mu1contin)
    })),
    file = "dvds_output/continuous_zsb.csv",
    append = (i != 1),
    row.names = F,
    col.names = (i==1)
  )
}
