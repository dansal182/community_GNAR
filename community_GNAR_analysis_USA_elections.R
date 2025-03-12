# Presidential election cycles analysis by using community-alpha GNAR models
library(GNAR) # Cobrit and Wagner plots
library(fastCorbit)
library(glmnet)
# community-alpha GNAR fitting code
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_gnar_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/communal_gnar_simulation_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/global_gnar_fitting.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_interactions_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/one-stepcommunity_GNAR_prediction.r")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/com_GNARfit_tv_weights.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/community_GNAR_simulation")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/community_GNARfit.r")
# PNACF level plots
source("/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/nacf_level_plots.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/latex_output_and_utilities/latex_outputs_utilities.R")
# Load USA network and R data
load("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/data_sets/USA_presidential_GNAR.RData")

# Community-wise decomposition for each vector time series
red_vts <- t(vapply(seq(1:12), function(x) {return(red_states * votes_vts_matrix[x, ])}, rep(0, 51)))
blue_vts <- t(vapply(seq(1:12), function(x) {return(blue_states * votes_vts_matrix[x, ])}, rep(0, 51)))
swing_vts <- t(vapply(seq(1:12), function(x) {return(swing_states * votes_vts_matrix[x, ])}, rep(0, 51)))

votes_vts_matrix_diff <- votes_vts_matrix[2:12, ] - votes_vts_matrix[1:11, ]
colnames(votes_vts_matrix_diff) <- usa_net_data$Tag
red_vts_diff <- t(vapply(seq(1:11), function(x) {return(red_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))
blue_vts_diff <- t(vapply(seq(1:11), function(x) {return(blue_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))
swing_vts_diff <- t(vapply(seq(1:11), function(x) {return(swing_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))


# Wagner plots
# NACF
# PNACF
r_corbit_plot_cpp(list(red_vts, blue_vts, swing_vts), list(usa_net_GNAR, usa_net_GNAR, usa_net_GNAR), 4, 3, list(W_red, W_blue, W_all), 
            c("Red States", "Blue States", "Swing States"), 
            same_net = "no", partial = "yes")

r_corbit_plot_cpp(list(red_vts_diff, blue_vts_diff, swing_vts_diff), list(usa_net_GNAR, usa_net_GNAR, usa_net_GNAR), 8, 6,  list(W_all, W_all, W_all), 
              c("Red States", "Blue States", "Swing States"), 
              same_net = "no", partial = "yes")
# PNACF Level plots
pnacf_level_plot(red_vts, usa_net_GNAR, 2, 11, W_red + W_blue)
pnacf_level_plot(blue_vts, usa_net_GNAR, 2, 11, W_blue + W_red)
pnacf_level_plot(swing_vts, usa_net_GNAR, 2, 3, W_all)

corbit_plot(vts = swing_vts, net = usa_net_GNAR, max_lag = 5, max_stage = 11)

# Standardise the vector time series
votes_vts_matrix_std <- vapply(seq(1:51), function(x) {return((votes_vts_matrix[, x] - mean(votes_vts_matrix[, x])) / sd(votes_vts_matrix[, x]))}, rep(0, 12))

colnames(votes_vts_matrix_std) <- usa_net_data$Tag



# Community GNAR fits
# Fit a community-alpha GNAR([2, 2, 2], {(1, 0), (1, 0), (1, 0)}, 3) to three state-wise communities
com_gnar_model <- community_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, c(2, 2, 2), list(c(0, 1), c(0, 1), c(1, 1)), 
                                     W_all, list(red_states, blue_states, swing_states))
summary(com_gnar_model)
# Fit a community-alpha GNAR([2, 2, 2], {(1, 0), (1, 0), (1, 0)}, 3) to three state-wise communities with interactions
com_gnar_model_with_interactions_std <- community_interaction_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, c(2, 2, 2), list(c(1, 1), c(1, 1), c(1, 1)), 
                                                                   W_all, list(red_states, blue_states, swing_states), list(c(2, 3), c(1, 3), c(1, 2)))
summary(com_gnar_model_with_interactions_std)
cat(sprintf(latex_table(cbind(com_gnar_model_with_interactions_std$coefficients, com_gnar_model_with_interactions_std$coefficients), 
                        names(com_gnar_model_with_interactions_std$coefficients))))
# non-standardised
com_gnar_model_with_interactions <- community_interaction_gnar_fit(votes_vts_matrix, usa_net_GNAR, c(2, 2, 2), list(c(1, 1), c(1, 1), c(2, 2)), 
                                                                   W_all, list(red_states, blue_states, swing_states), list(c(2, 3), c(1, 3), c(1, 2)))
summary(com_gnar_model_with_interactions)
# global-alpha for comparison
summary(global_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, 2, c(0, 1), W_all))
# Analysis for one-lag differenced data
votes_vts_matrix_diff <- votes_vts_matrix[2:12, ] - votes_vts_matrix[1:11, ]
colnames(votes_vts_matrix_diff) <- usa_net_data$Tag
red_vts_diff <- t(vapply(seq(1:11), function(x) {return(red_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))
blue_vts_diff <- t(vapply(seq(1:11), function(x) {return(blue_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))
swing_vts_diff <- t(vapply(seq(1:11), function(x) {return(swing_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))


# Standardise series
votes_vts_matrix_diff_std <- vapply(seq(1:51), function(x) {return((votes_vts_matrix_diff[, x] - mean(votes_vts_matrix_diff[, x])) / sd(votes_vts_matrix_diff[, x]))}, 
                                    rep(0, 11))
colnames(votes_vts_matrix_diff_std) <- usa_net_data$Tag
red_vts_std <- t(vapply(seq(1:11), function(x) {return(red_states * votes_vts_matrix_diff_std[x, ])}, rep(0, 51)))
blue_vts_std <- t(vapply(seq(1:11), function(x) {return(blue_states * votes_vts_matrix_diff_std[x, ])}, rep(0, 51)))
swing_vts_std <- t(vapply(seq(1:11), function(x) {return(swing_states * votes_vts_matrix_diff_std[x, ])}, rep(0, 51)))

r_corbit_plot_cpp(list(red_vts_diff, blue_vts_std, swing_vts_std), list(usa_net_GNAR), 8, 6, list(W_red, W_blue, W_swing), c("Red States", "Blue States", "Swing States"), 
            same_net = "yes")
r_corbit_plot_cpp(list(red_vts_diff, blue_vts_diff, swing_vts_diff), list(usa_net_GNAR, usa_net_GNAR, usa_net_GNAR), 8, 6, list(W_all, W_all, W_all), c("Red States", "Blue States", "Swing States"), 
            same_net = "no", partial = "yes", viridis_color_option = "cividis")
pnacf_level_plot(red_vts_diff, usa_net_GNAR, 2, 11, W_red)
pnacf_level_plot(blue_vts_diff, usa_net_GNAR, 2,11, W_blue)
pnacf_level_plot(swing_vts_diff, usa_net_GNAR, 2, 11, W_all)



# Fitting Communal-alpha GNAR models with different orders
com_gnar_model_diff <-  community_gnar_fit(votes_vts_matrix_diff_std, usa_net_GNAR, c(3, 3, 3), list(c(1, 0, 0), c(1, 0, 0), c(0, 1, 0)), W_all, 
                                           list(red_states, blue_states, swing_states))
summary(com_gnar_model_diff)

cat(sprintf(latex_table(cbind(com_gnar_model_diff$coefficients, com_gnar_model_diff$coefficients), 
                        names(com_gnar_model_diff$coefficients))))

com_gnar_model_with_interactions_diff <- community_interaction_gnar_fit(votes_vts_matrix_diff_std, usa_net_GNAR, c(3, 4, 2), list(c(0, 0, 0), c(1, 0, 0, 0), c(1, 1)), 
                                                                        W_all, list(red_states, blue_states, swing_states), list(c(0), c(0), c(1, 2)))
summary(com_gnar_model_with_interactions_diff)

summary(global_gnar_fit(votes_vts_matrix_diff_std, usa_net_GNAR, 3, c(0, 1, 0), W_all))
sum(abs(global_gnar_fit(votes_vts_matrix_diff_std, usa_net_GNAR, 3, c(0, 1, 0), W_all)$coefficients))

# Fitting ridge-like models 

com_gnar_model_with_interactions_std <- community_interaction_GNAR_RIDGEfit(votes_vts_matrix_std, usa_net_GNAR, c(5, 5, 5), list(rep(2, 5), rep(2, 5), rep(2, 5)), 
                                                                       W_all, list(red_states, blue_states, swing_states), list(c(2, 3), c(1, 3), c(1, 2)))
coef(com_gnar_model_with_interactions_std)
