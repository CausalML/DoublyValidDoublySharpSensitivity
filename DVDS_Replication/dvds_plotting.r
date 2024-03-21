library(dplyr)
library(ggplot2)
set.seed(1)

# Plots DVDS simulation results

# =============================================
# Reading data
binary_dvds = read.table("Simulations/dvds_output/binary_dvds.csv", header=T) %>% 
  mutate(outcome = "Binary", method = "DVDS")
binary_zsb = read.table("Simulations/dvds_output/binary_zsb.csv", header=T) %>% 
  mutate(outcome = "Binary", method = "ZSB")
continuous_dvds = read.table("Simulations/dvds_output/continuous_dvds.csv", header=T) %>%
  mutate(outcome = "Continuous", method = "DVDS")
continuous_zsb = read.table("Simulations/dvds_output/continuous_zsb.csv", header=T) %>%
  mutate(outcome = "Continuous", method = "ZSB")
df = bind_rows(binary_dvds, binary_zsb, continuous_dvds, continuous_zsb)


# =============================================
# Computing sharp bounds

# Binary bounds via simulation
Lambdas = as.numeric(unique(df$Lambda))
n = 500000
d = 5
X = matrix(runif(n*d, -1, 1), nrow = n)
e = 1/(1 + exp(X[,1] + 0.5*(X[,2]>0) + 0.5*X[,2]*X[,3]))
mu = 1/(1 + exp(0.5*X[,1] + X[,2] + 0.25*X[,2]*X[,3]))
true_upper = sapply(Lambdas, function(l) { 
  psi1plus = mean(e*mu + (1-e)*pmin(1 - 1/l + mu/l, mu*l)) 
  psi0mins = mean((1-e)*mu + e*pmax(1 - l + mu*l, mu/l))
  psi1plus - psi0mins
})
true_lower = sapply(Lambdas, function(l) {
  psi1mins = mean(e*mu + (1-e)*pmax(1 - l + mu*l, mu/l))
  psi0plus = mean((1-e)*mu + e*pmin(1 - 1/l + mu/l, mu*l))
  psi1mins - psi0plus
})
binary_truth = data.frame("outcome" = "Binary", "Lambda" = Lambdas,
                          "true_upper" = true_upper,
                          "true_lower" = true_lower)

# Continuous bounds via explicit formulas
continuous_truth = data.frame("outcome" = "Continuous", "Lambda" = Lambdas,
                              "true_upper" = (1+1/3)*(Lambdas^2 - 1)/Lambdas*dnorm(qnorm(Lambdas/(Lambdas+1))),
                              "true_lower" =-(1+1/3)*(Lambdas^2 - 1)/Lambdas*dnorm(qnorm(Lambdas/(Lambdas+1))))
truebounds = bind_rows(binary_truth, continuous_truth)
df = left_join(df, truebounds, by=c("Lambda", "outcome"))


# =============================================
# Plotting results
pdf(file = "Simulations/dvds_output/art/binary_comparison.pdf", width=6, height=4)
df %>%
  filter(outcome == "Binary") %>%
  mutate(method = ifelse(method == "DVDS", "Binary DVDS", "Binary ZSB")) %>%
  mutate(Lambda = round(Lambda, 2)) %>%
  ggplot() + 
  geom_boxplot(aes(x = factor(Lambda), y = upper, col = factor(method), fill = factor(method)), alpha = 0.25) +
  geom_boxplot(aes(x = factor(Lambda), y = lower, col = factor(method), fill = factor(method)), alpha = 0.25) +
  geom_line(aes(x = as.numeric(factor(Lambda)), y = true_lower), size = 1) +
  geom_line(aes(x = as.numeric(factor(Lambda)), y = true_upper), size = 1) +
  facet_wrap(.~ method) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "None"
  ) +
  scale_x_discrete(breaks = levels(factor(round(df$Lambda, 2)))[c(1, 6, 11, 16, 20)]) +
  xlab(expression(Lambda)) +
  ylab("")
dev.off()

pdf(file = "Simulations/dvds_output/art/continuous_comparison.pdf", width=6, height=4)
df %>%
  filter(outcome == "Continuous") %>%
  mutate(method = ifelse(method == "DVDS", "Continuous DVDS", "Continuous ZSB")) %>%
  mutate(Lambda = round(Lambda, 2)) %>%
  ggplot() + 
  geom_boxplot(aes(x = factor(Lambda), y = upper, col = factor(method), fill = factor(method)), alpha = 0.25) +
  geom_boxplot(aes(x = factor(Lambda), y = lower, col = factor(method), fill = factor(method)), alpha = 0.25) +
  geom_line(aes(x = as.numeric(factor(Lambda)), y = true_lower), size = 1) +
  geom_line(aes(x = as.numeric(factor(Lambda)), y = true_upper), size = 1) +
  facet_wrap(.~ method) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "None"
  ) +
  scale_x_discrete(breaks = levels(factor(round(df$Lambda, 2)))[c(1, 6, 11, 16, 20)]) +
  xlab(expression(Lambda)) +
  ylab("") 
dev.off()

# =============================================
# Compute coverage
coverage = df %>%
  filter(method == "DVDS") %>%
  group_by(Lambda, outcome) %>%
  mutate(covered = (lower.CI < true_lower & true_upper < upper.CI)) %>%
  summarise(coverage = mean(covered))

binary.coverage     = coverage %>% filter(outcome == "Binary") %>% pull(coverage)
continuous.coverage = coverage %>% filter(outcome == "Continuous") %>% pull(coverage)

names(binary.coverage)     = round(unique(df$Lambda), 2)
names(continuous.coverage) = round(unique(df$Lambda), 2)


coveragePrint = function(cover.rates) {
  return(paste0(
    "coverage ranging from $",
    round(100 * min(cover.rates), 1),
    "\\%$ (when $\\Lambda = ",
    names(cover.rates)[which.min(cover.rates)], 
    "$) to $", 
    round(100 * max(cover.rates), 1),
    "\\%$ (when $\\Lambda = ", 
    names(cover.rates)[which.max(cover.rates)],
    "$) with average coverage of $", 
    round(100 * mean(cover.rates), 1), 
    "\\%"
  ))
}


fileConn <- file("Simulations/dvds_output/binary_coverageRates.tex")
writeLines(coveragePrint(binary.coverage), fileConn)
close(fileConn)


fileConn <- file("Simulations/dvds_output/continuous_coverageRates.tex")
writeLines(coveragePrint(continuous.coverage), fileConn)
close(fileConn)





