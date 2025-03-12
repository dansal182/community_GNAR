# load the required scripts

library(GNAR)
library(igraph)
library(sparsevar)
library(CARBayesST)
# We centre the data before fitting all models, i.e., subtract the mean to remove some trend.
votes_vts_matrix_std <- vapply(seq(1:51), function(x) {return((votes_vts_matrix[1:11, x] - mean(votes_vts_matrix[1:11, x])) / sd(votes_vts_matrix[1:11, x]))}, rep(0, 11))
colnames(votes_vts_matrix_std) <- usa_net_data$Tag
corbit_plot_cpp(vts = votes_vts_matrix, net = usa_net_GNAR, max_lag = 7, max_stage = 3, partial = "yes", viridis_color_option = "cividis")
#############################################################################
# community-alpha GNAR model
com_model_two <- community_interaction_gnar_fit(votes_vts_matrix_std[1:11, ], usa_net_GNAR, c(2, 2, 2), list(c(2, 2), c(2, 2), c(2, 2)), 
                                                W_all, list(red_states, blue_states, swing_states), list(c(0), c(0), c(0)))
summary(com_model_two)
hist(com_model_two$residuals)
one_step_ahead = community_GNARpred(usa_net_GNAR, W_all, list(c(2, 2, 2), list(c(2, 2), c(2, 2), c(2, 2))), com_model_two$coefficients, 
                                    votes_vts_matrix_std[10:11, ],
                                    list(red_states, blue_states, swing_states), list(c(0), c(0), c(0)))
names(one_step_ahead) <- usa_net_data$Tag
one_step_ahead
#round(sqrt(sum((one_step_ahead - votes_vts_matrix_std[12, ])^2) / 51), digits = 3)
xhat <- vapply(seq(1:51), function(x) {return(sd(votes_vts_matrix[1:11, x]) * one_step_ahead[x] + mean(votes_vts_matrix[1:11, x]))}, 0.0)
round(sqrt(sum((xhat - votes_vts_matrix[12, ])^2) / 51), digits = 3)
##################################################################################

# community-alpha with interactions GNAR model
com_model_two_noncentre_interactions <- community_interaction_gnar_fit(votes_vts_matrix_std[1:11, ], usa_net_GNAR, c(2, 2, 2), list(c(2, 2), c(2, 2), c(2, 2)), 
                                                                       W_all, list(red_states, blue_states, swing_states), list(c(2, 3), c(1, 3), c(1, 2)))
summary(com_model_two_noncentre_interactions)
hist(com_model_two_noncentre_interactions$residuals)
one_step_ahead = community_GNARpred(usa_net_GNAR, W_all, list(c(2, 2, 2), list(c(2, 2), c(2, 2), c(2, 2))), com_model_two_noncentre_interactions$coefficients, 
                                    votes_vts_matrix_std[10:11, ],
                                    list(red_states, blue_states, swing_states), list(c(2, 3), c(1, 3), c(1, 2)))
names(one_step_ahead) <- usa_net_data$Tag
one_step_ahead
#round(sqrt(sum((one_step_ahead - votes_vts_matrix_std[12, ])^2) / 51), digits = 3)
xhat <- vapply(seq(1:51), function(x) {return(sd(votes_vts_matrix[1:11, x]) * one_step_ahead[x] + mean(votes_vts_matrix[1:11, x]))}, 0.0)
round(sqrt(sum((xhat - votes_vts_matrix[12, ])^2) / 51), digits = 3)
################################################################################

# Sparse VAR model
set.seed(2024)
sparse_var_forecast <- fitVAR(votes_vts_matrix_std[1:11, ] - colMeans(votes_vts_matrix_std[1:11, ]), p = 2)
xhat =  sparse_var_forecast$A[[1]] %*% votes_vts_matrix_std[11, ] + 
  sparse_var_forecast$A[[2]] %*% votes_vts_matrix_std[10, ]
#round(sqrt(sum((xhat - votes_vts_matrix_std[12, ])^2) / 51), digits = 3)
b = sparse_var_forecast$A[[1]]
b[abs(b) > 0.0] = 1
sum(rowSums(b))
b = sparse_var_forecast$A[[2]]
b[abs(b) > 0.0] = 1
sum(rowSums(b))
yhat <- vapply(seq(1:51), function(x) {return(sd(votes_vts_matrix[1:11, x]) * xhat[x] + mean(votes_vts_matrix[1:11, x]))}, 0.0)
round(sqrt(sum((yhat - votes_vts_matrix[12, ])^2) / 51), digits = 3)
##################################################################################

# global-alpha GNAR Model
global_alpha <- GNARfit(votes_vts_matrix_std[1:11, ], usa_net_GNAR, alphaOrder = 2, betaOrder = c(2, 2))
summary(global_alpha)
#one_step_global_alpha = global_alpha$mod$coefficients[1] * votes_vts_matrix_std[11, ] +
#                                                                                                    global_alpha$mod$coefficients[2] * (as.matrix(S1) * W_all) %*% votes_vts_matrix_std[11, ] +
#  global_alpha$mod$coefficients[3] * votes_vts_matrix_std[11, ]                                 
#round(sqrt(sum((predict(global_alpha)- votes_vts_matrix_std[12, ])^2) / 51), digits = 3)
one_step_ahead = predict(global_alpha)
xhat <- vapply(seq(1:51), function(x) {return(sd(votes_vts_matrix[1:11, x]) * one_step_ahead[x] + mean(votes_vts_matrix[1:11, x]))}, 0.0)
round(sqrt(sum((xhat - votes_vts_matrix[12, ])^2) / 51), digits = 3)
# local-alpha GNAR model
local_alpha <- GNARfit(votes_vts_matrix_std[1:11, ] - colMeans(votes_vts_matrix_std[1:11, ]), usa_net_GNAR, alphaOrder = 2, betaOrder = c(2, 2), globalalpha = FALSE)
summary(local_alpha)
one_step_ahead = predict(local_alpha)
#  colMeans(votes_vts_matrix_std[1:11, ]) + local_alpha$mod$coefficients[1:51] * (votes_vts_matrix_std[11, ] - colMeans(votes_vts_matrix_std[1:11, ])) +
#  local_alpha$mod$coefficients[52] * (as.matrix(S1) * W_all) %*% (votes_vts_matrix_std[11, ] - colMeans(votes_vts_matrix_std[1:11, ])) +
#  local_alpha$mod$coefficients[53:103] * (votes_vts_matrix_std[11, ] - colMeans(votes_vts_matrix_std[1:11, ]))  
#round(sqrt(sum((one_step_ahead - votes_vts_matrix_std[12, ])^2) / 51), digits = 3)
xhat <- vapply(seq(1:51), function(x) {return(sd(votes_vts_matrix[1:11, x]) * one_step_ahead[x] + mean(votes_vts_matrix[1:11, x]))}, 0.0)
round(sqrt(sum((xhat - votes_vts_matrix[12, ])^2) / 51), digits = 3)
###########################################

# spatio-temporal CAR.ar model

votes_vector_form <- rep(0, 11 * 51)

for (i in 1:11) {
  top = 1 + 51 * (i - 1)
  bottom = 51 + 51 * (i - 1)
  votes_vector_form[top:bottom] = votes_vts_matrix_std[i, ]
}
names(votes_vector_form) <- rep(usa_net_data$Tag, 11)
votes_vector_frame <- data.frame(votes_vector_form)
S1 = as.matrix(usa_net)
# synthetically connect Alaska and Hawaii
S1[2, 11] = 1
S1[11, 2] = 1
set.seed(2024)
car_ar_model <- ST.CARar(formula = votes_vector_form ~ 1, family = "gaussian", data = votes_vector_frame,
                         W = as.matrix(S1), AR = 2, burnin = 20000, n.sample = 220000, thin = 10) 

car_ar_model$summary.results[1, 1]
car_ar_model$summary.results[5, 1]
car_ar_model$summary.results[6, 1]

car_ar_model$summary.results

one_step_car_ar = car_ar_model$summary.results[1, 1] + car_ar_model$summary.results[5, 1] * (votes_vts_matrix_std[11, ] - car_ar_model$summary.results[1, 1]) +
  car_ar_model$summary.results[6, 1] * (votes_vts_matrix_std[10, ] - car_ar_model$summary.results[1, 1])
# round(sqrt(sum((one_step_car_ar - votes_vts_matrix_std[12, ])^2) / 51), digits = 3)
xhat <- vapply(seq(1:51), function(x) {return(sd(votes_vts_matrix[1:11, x]) * one_step_car_ar[x] + mean(votes_vts_matrix[1:11, x]))}, 0.0)
round(sqrt(sum((xhat - votes_vts_matrix[12, ])^2) / 51), digits = 3)

##########################
# naive, i.e., use the previous observation as forecast value

#round(sqrt(sum((votes_vts_matrix_std[11, ] - votes_vts_matrix_std[12, ])^2) / 51), digits = 3)
xhat <- vapply(seq(1:51), function(x) {return(sd(votes_vts_matrix[1:11, x]) * votes_vts_matrix_std[11, x] + mean(votes_vts_matrix[1:11, x]))}, 0.0)
round(sqrt(sum((xhat - votes_vts_matrix[12, ])^2) / 51), digits = 3)
