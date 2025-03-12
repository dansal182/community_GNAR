# Simulation study presented in the community-alpha GNAR paper
# We illustrate model selection, fitting methods and interpretation for a community-alpha with interactions GNAR model using the USA 
# state-wise network as underlying structure, further we illustrate the performance of the NACF, PNACF and our estimators with missing
# data. 
library(GNAR)
library(igraph)
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/community_GNAR_simulation")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_gnar_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/communal_gnar_simulation_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/global_gnar_fitting.R")
# PNACF level plots
source("/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/nacf_level_plots.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_interactions_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/com_GNARfit_tv_weights.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/weight_adjustment_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/nacf_pnacf_edited.r")
load("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/data_sets/USA_presidential_GNAR.RData")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/one-stepcommunity_GNAR_prediction.r", encoding = "UTF-8")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/corbit_plot_edit.r")
# The simulation study is a synthetic longitudinal study of public policies that have a period of four-years.
# We compute the time varying weights for the USA state-wise network

red_weight_function <- function(d_ij, t) {
    return ((1 + cos(t * pi / 2) + 0.1) * 2^(-d_ij / 2))
}

blue_weight_function <- function(d_ij, t) {
  return ((1 + sin(t * pi / 2) + 0.1) / 2^(d_ij))
}

swing_weight_function <- function(d_ij, t) {
  return ((1 + sin(t * pi / 4) * cos(t * pi / 4) + 0.1) / 2^(d_ij))
}

weight_function <- function(d_ij, j_community, t) {
  if (!is.finite(d_ij) || d_ij == 0) {
    out = 0.0
  } else if (j_community == 'red') {
    out = red_weight_function(d_ij, t)
  } else if (j_community == "blue") {
    out = blue_weight_function(d_ij, t)
  } else {
    out = swing_weight_function(d_ij, t)
  }
  return (out)
}

time_step_matrix <- function(t, distance_matrix) {
  d = ncol(distance_matrix)
  out = vapply(seq(1:d), function(x) {row_matrix(t, distance_matrix, x)}, rep(0, d))
  # out = out - diag(d)
  return (out)
}

row_matrix <- function(t, distance_matrix, ith_row) {
  d = ncol(distance_matrix)
  row = vapply(seq(1:d), function(x) {weight_function(distance_matrix[ith_row, x], state_communities[x], t)}, 0.0)
  return (row)
}

# Plot the USA state-wise network
plot(usa_net)

# Compute the node distance matrix
USA_net_distances = distances(usa_net)


# Compute community indicator 
state_communities = vapply(seq(1:51), function(x) {ifelse(red_states[x] == 1, "red", ifelse(blue_states[x] == 1, "blue", "swing"))}, '')

# Compute the weight matrix at each time step from 1976 to 2020, i.e., 44 time steps (years)
time_steps = 100
# Specify model order for the simulation
# Different lag orders for each community
alpha_order = c(3, 2, 1)
# r-stage neighbourhood and interaction term regression order
beta_orders = list(c(2, 2, 1), c(3, 2), c(3))
# interaction orders
interaction_order = list(c(3), c(3), c(1, 2))

# Set seed for our simulation
sample_sizes = c(c(10:1000), 10000)
l2_errors = matrix(rep(0, length(sample_sizes) * 4), nrow = 4, ncol = length(sample_sizes))
colnames(l2_errors) <- vapply(c(1:length(sample_sizes)), function(x) {paste0("T=", as.character(sample_sizes[x]))}, "")
rownames(l2_errors) <- c("K1", "K2", "K3", "K")
i = 1
for (t in sample_sizes) {
  set.seed(2024)
  W_time_varying_list = lapply(seq(1:t), function(x) {time_step_matrix(x, USA_net_distances)})
# simulate a staitonary idealised vector of parameters
  theta0 = stationary_theta_sim(alpha_order, beta_orders, interaction_order,  c(0.25, 0.22, 0.28))

  alpha_params0 = list(c(theta0[1], theta0[6], theta0[11]), c(theta0[14], theta0[21]), c(theta0[26]))

  beta_params0 = list(list(c(theta0[2], theta0[3]), c(theta0[7], theta0[8]), c(theta0[12])), 
                    list(c(theta0[15], theta0[16], theta0[17]), c(theta0[22], theta0[23])),
                    list(c(theta0[27:29])))

  gamma_params0 = list(list(c(theta0[4:5]), c(theta0[9:10]), c(theta0[13])),
                     list(c(theta0[18:20]), c(theta0[24:25])),
                     list(c(theta0[30:35])))

# Simulate a stationary community-alpha GNAR model with interactions
  simulate_longitudinal_data <- community_GNARsim(t, usa_net_GNAR, alpha_params0, beta_params0, gamma_params0, W_time_varying_list,
                                                list(red_states, blue_states, swing_states), interaction_order)
  colnames(simulate_longitudinal_data) <- usa_net_data$Tag
  

# Fitting the model with oracle knowledge of model order


  longitudinal_fitted_model <- com_GNARfit_tv_weights(simulate_longitudinal_data, usa_net_GNAR, alpha_order, beta_orders, W_time_varying_list,
                                                    list(red_states, blue_states, swing_states), interaction_order)
  
  # red states community l2-error bound
  l2_errors[1, i] = sqrt(sum((longitudinal_fitted_model$coefficients[1:13] - theta0[1:13])^2))
  # blue states l2-error bound
  l2_errors[2, i] = sqrt(sum((longitudinal_fitted_model$coefficients[14:25] - theta0[14:25])^2))
  # swing states l2-error bound
  l2_errors[3, i] = sqrt(sum((longitudinal_fitted_model$coefficients[26:35] - theta0[26:35])^2))
  # complete model l2-error bound
  l2_errors[4, i] = sqrt(sum((longitudinal_fitted_model$coefficients - theta0)^2))
  if (i %% 10 == 0) {
    cat(paste0("done: ", as.character(i), "\n"))
  }
  # update counter
  i = i + 1
}


write.csv(data.frame(l2_errors, check.names = FALSE), "/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/code_scripts/l2-error_norm.csv")


# random parameter estimator

random_errors = vapply(seq(1:10000), function(x) {
  return(sqrt(sum((stationary_theta_sim(alpha_order, beta_orders, interaction_order,  runif(n=3, min=0, max=1)) - theta0)^2)))},
                      0.0)
print(mean(random_errors))
print(sd(random_errors))
print(hist(random_errors))

l2_errors_data <- read.csv("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/code_scripts/l2-error_norm.csv", 
                           header=TRUE)
plot(c(1:991), l2_errors_data[1, 2:992], 'l', col='orange', ylim=c(0, 1), main = expression("Plot of ||"*hat(theta) - theta[0] * "||"[2]), xlab = "realisation length: T", 
     ylab="error-norm")+
  lines(c(1:991), l2_errors_data[2, 2:992], col='blue')+
  lines(c(1:991), l2_errors_data[3, 2:992], col='darkgreen')+
  lines(c(1:991), l2_errors_data[4, 2:992], col='black')
