# Testing the phase-transition of community alpha GNAR models
# Load required packages and datasets
library(GNAR)
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_gnar_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/communal_gnar_simulation_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/global_gnar_fitting.R")
# PNACF level plots
source("/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/nacf_level_plots.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_interactions_fitting_methods.R")
# fiveNet basic network time series blocks
W = weights_matrix(fiveNet, 3)
fiveNet_stages <- get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(fiveNet)), 3)

# Testing out with global-alpha model

simulated_vts <- GNARsim(50, fiveNet, alphaParams = list(rep(0.37, 5)), betaParams = list(c(0.24, 0.18)))
sd(simulated_vts[, 2])
linear_model_test = build_linear_model_edt(fiveNet_stages, 1, c(2), simulated_vts, W)

eigen(t(linear_model_test[, 2:4]) %*% linear_model_test[, 2:4])

summary(lm(linear_model_test[, 1] ~ linear_model_test[, 2:4] + 0))
sqrt(sum((lm(linear_model_test[, 1] ~ linear_model_test[, 2:4] + 0)$coefficients - c(0.37, 0.24, 0.18))^2))

n_sample = nrow(linear_model_test)
tau_test = min(eigen(t(linear_model_test[, 2:4]) %*% linear_model_test[, 2:4])$values) / n_sample
dummy_bound = 2 * sqrt(3) / (tau_test * n_sample) * max(abs(t(linear_model_test[, 2:4]) %*% rnorm(n_sample)))
print(dummy_bound)
delta = 0.10
up_bound = 2 * sqrt(3) / tau_test * sqrt(2 * log(3) / n_sample  + delta)
print(up_bound)
bound_prob = 1 - 2 * exp(- delta^2 * n_sample * tau_test / 2)
print(bound_prob)
sample_sizes = 1000
bound_points = 1 - 2 * exp(- delta^2 * c(1:sample_sizes) * tau_test / 2)
bound_points = vapply(seq(1:sample_sizes), function(x) {ifelse(bound_points[x] >= 0, bound_points[x], 0)}, 0.0)
plot(1:sample_sizes, bound_points, "l")

