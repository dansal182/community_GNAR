# load the required scripts
library(GNAR)
library(igraph)
library(sparsevar)
library(CARBayesST)
# We centre the data before fitting all models, i.e., subtract the mean to remove some trend.
corbit_plot(vts = votes_vts_matrix, net = usa_net_GNAR, max_lag = 8, max_stage = 5, viridis_color_option = "cividis", partial="yes")
#######
# community-alpha GNAR model
com_model_two_noncentre <- com_GNARfit_tv_weights(votes_vts_matrix[1:11, ], usa_net_GNAR, c(2, 2, 2), list(c(2, 2), c(2, 2), c(2, 2)), 
                                        lapply(seq(1:12), function(x) {return(W_all)}), list(red_states, blue_states, swing_states), 
                                        list(c(0), c(0), c(0)))
summary(com_model_two_noncentre)
hist(com_model_two_noncentre$residuals)
one_step_ahead = community_GNARpred(usa_net_GNAR, W_all, list(c(2, 2, 2), list(c(2, 2), c(2, 2), c(2, 2))), com_model_two_noncentre$coefficients, 
                                    votes_vts_matrix[10:11, ],
                                    list(red_states, blue_states, swing_states), list(c(0), c(0), c(0)))
names(one_step_ahead) <- usa_net_data$Tag
one_step_ahead
round(sqrt(sum((one_step_ahead - votes_vts_matrix[12, ])^2) / 51), digits = 3)

#########################################################################################################################

# community-alpha with interactions GNAR model
com_model_two_noncentre_interactions <- community_interaction_gnar_fit(votes_vts_matrix[1:11, ], usa_net_GNAR, c(2, 2, 2), list(c(2, 2), c(2, 2), c(2, 2)), 
                                                                       W_all, list(red_states, blue_states, swing_states), list(c(2, 3), c(1, 3), c(1, 2)))
summary(com_model_two_noncentre_interactions)
hist(com_model_two_noncentre_interactions$residuals)
one_step_ahead = community_GNARpred(usa_net_GNAR, W_all, list(c(2, 2, 2), list(c(2, 2), c(2, 2), c(2, 2))), com_model_two_noncentre_interactions$coefficients, 
                                    votes_vts_matrix[10:11, ],
                                    list(red_states, blue_states, swing_states), list(c(2, 3), c(1, 3), c(1, 2)))
names(one_step_ahead) <- usa_net_data$Tag
one_step_ahead
round(sqrt(sum((one_step_ahead - votes_vts_matrix[12, ])^2) / 51), digits = 3)

##################################################################################
# Sparse VAR model
set.seed(2024)
sparse_var_forecast <- fitVAR(votes_vts_matrix[1:11, ], p = 2)
xhat = sparse_var_forecast$A[[1]] %*% votes_vts_matrix[11, ] + sparse_var_forecast$A[[2]] %*% votes_vts_matrix[10, ]
round(sqrt(sum((xhat - votes_vts_matrix[12, ])^2) / 51), digits = 3)
b = sparse_var_forecast$A[[1]]
b[abs(b) > 0.0] = 1
sum(rowSums(b))
b = sparse_var_forecast$A[[2]]
b[abs(b) > 0.0] = 1
sum(rowSums(b))

##################################################################################

# global-alpha GNAR Model
global_alpha <- GNARfit(votes_vts_matrix[1:11, ], usa_net_GNAR, alphaOrder = 2, betaOrder = c(1, 0))
summary(global_alpha)
one_step_global_alpha = predict(global_alpha)
round(sqrt(sum((one_step_global_alpha - votes_vts_matrix[12, ])^2) / 51), digits = 3)

# local-alpha GNAR model
local_alpha <- GNARfit(votes_vts_matrix[1:11, ], usa_net_GNAR, alphaOrder = 2, betaOrder = c(1, 0), globalalpha = FALSE)
summary(local_alpha)
one_step_local_alpha = predict(local_alpha)
round(sqrt(sum((one_step_local_alpha - votes_vts_matrix[12, ])^2) / 51), digits = 3)

###########################################

# spatio-temporal CAR.ar model

votes_vector_form <- rep(0, 11 * 51)

for (i in 1:11) {
  top = 1 + 51 * (i - 1)
  bottom = 51 + 51 * (i - 1)
  votes_vector_form[top:bottom] = votes_vts_matrix[i, ]
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
car_ar_model

one_step_car_ar = 0.8439 * votes_vts_matrix[11, ] + 0.1121 * votes_vts_matrix[10, ]

round(sqrt(sum((one_step_car_ar - votes_vts_matrix[12, ])^2) / 51), digits = 3)

one_step_car_ar = car_ar_model$summary.results[1, 1] + car_ar_model$summary.results[5, 1] * (votes_vts_matrix[11, ] - car_ar_model$summary.results[1, 1]) +
  car_ar_model$summary.results[6, 1] * (votes_vts_matrix[10, ] - car_ar_model$summary.results[1, 1])
round(sqrt(sum((one_step_car_ar - votes_vts_matrix[12, ])^2) / 51), digits = 3)

##########################
# naive, i.e., use the previous observation as forecast value

round(sqrt(sum((votes_vts_matrix[11, ] - votes_vts_matrix[12, ])^2) / 51), digits = 3)


save(usa_net_data, usa_net, usa_net_GNAR, rep_data, republican_votes_vts, votes_vts_matrix, red_states, blue_states,
     swing_states, W_all, W_red, W_blue, W_swing, one_step_ahead_noncentre, one_step_ahead,
     file = "/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/data_sets/Daniel_Salnikov_484.RData")
load("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/data_sets/Daniel_Salnikov_486.RData")

colMeans(votes_vts_matrix)
hist(red_vts[red_vts > 0])
hist(blue_vts[blue_vts > 0])
hist(swing_vts[swing_vts > 0])
mean(red_vts[red_vts > 0])
mean(blue_vts[blue_vts > 0])
mean(swing_vts[swing_vts > 0])
