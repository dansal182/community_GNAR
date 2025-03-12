# asymptotic normality simulations
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
# simulate a base vector of parameters

# Different lag orders for each community
alpha_order = c(3, 2, 1)
# r-stage neighbourhood and interaction term regression order
beta_orders = list(c(2, 2, 1), c(3, 2), c(3))
# interaction orders
interaction_order = list(c(3), c(3), c(1, 2))

theta0 = stationary_theta_sim(alpha_order, beta_orders, interaction_order,  c(0.25, 0.22, 0.28))

alpha_params0 = list(c(theta0[1], theta0[6], theta0[11]), c(theta0[14], theta0[21]), c(theta0[26]))

beta_params0 = list(list(c(theta0[2], theta0[3]), c(theta0[7], theta0[8]), c(theta0[12])), 
                    list(c(theta0[15], theta0[16], theta0[17]), c(theta0[22], theta0[23])),
                    list(c(theta0[27:29])))

gamma_params0 = list(list(c(theta0[4:5]), c(theta0[9:10]), c(theta0[13])),
                     list(c(theta0[18:20]), c(theta0[24:25])),
                     list(c(theta0[30:35])))


sample_sizes = c(10, 20, 60, 100, 1000, 10000)

nsims = 1000

q = length(theta0)
W_basic = weights_matrix(usa_net_GNAR, 11)
W_list = lapply(seq(1:sample_sizes[2]), function(x) {return(W_basic)})


fit_simulated_model <- function(time_steps, x) {
  sim_data = community_GNARsim(time_steps, usa_net_GNAR, alpha_params0, beta_params0, gamma_params0, W_list,
                               list(red_states, blue_states, swing_states), interaction_order)
  
  longitudinal_fitted_model <- com_GNARfit_tv_weights(sim_data, usa_net_GNAR, alpha_order, beta_orders, W_list,
                                                      list(red_states, blue_states, swing_states), interaction_order)
  # l2-errors
  l2_errors = longitudinal_fitted_model$coefficients - theta0
  print(paste("done:", as.character(x)))
  return (l2_errors)
}


simulated_l2_errors <- vapply(c(1:nsims), function(x) {fit_simulated_model(sample_sizes[2], x)}, rep(0, q))
max(colMeans(simulated_l2_errors))
min(colMeans(simulated_l2_errors))
mean(colMeans(simulated_l2_errors))
sd(simulated_l2_errors)

hist(simulated_l2_errors[, 33])

hist(l2_errors[, 33])



