# Analysis of presidential elections in the USA from 1976 to 2020
# We begin by fitting a community-alpha GNAR([2, 2, 2], {[1, 1], [1, 1], [1, 1]}, {(2, 3), (1, 3), (1, 2)}, [3])

votes_vts_matrix_std <- vapply(seq(1:51), function(x) {return((votes_vts_matrix[, x] - mean(votes_vts_matrix[, x])) / sd(votes_vts_matrix[, x]))}, rep(0, 12))
votes_vts_matrix = votes_vts_matrix / 100
colnames(votes_vts_matrix_std) <- usa_net_data$Tag

# Community GNAR fits
# Fit a community-alpha GNAR([2, 2, 2], {(1, 0), (1, 0), (1, 0)}, 3) to three state-wise communities
com_gnar_model <- community_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, c(2, 2, 2), list(c(0, 1), c(0, 1), c(1, 1)), 
                                     W_all, list(red_states, blue_states, swing_states))
summary(com_gnar_model)
# Fit a community-alpha GNAR([2, 2, 2], {(1, 0), (1, 0), (1, 0)}, 3) to three state-wise communities with interactions
com_gnar_model_with_interactions <- com_GNARfit_tv_weights(votes_vts_matrix[1:11, ], usa_net_GNAR, c(2), list(c(1, 1)), 
                                                                       lapply(seq(1:11), function(x) {return(W_all)}), list(rep(1, 51)), 
                                                               list(c(0)))
summary(com_gnar_model_with_interactions)
hist(com_gnar_model_with_interactions$residuals)

one_step_ahead = 55.48 + (0.861174 * diag(51) + 0.007583 * W_all * as.matrix(S1)) %*% (votes_vts_matrix[11, ] - 55.48) + 
  (0.106524 * diag(51) - 0.008305 * W_all * as.matrix(S1)) %*% (votes_vts_matrix[10, ] - 55.48)
one_step_ahead = community_GNARpred(usa_net_GNAR, W_all, list(c(2), list(c(1, 1))), com_gnar_model_with_interactions$coefficients, votes_vts_matrix[1:11, ],
                   list(rep(1, 51)), list(c(0)))
names(one_step_ahead) <- usa_net_data$Tag
one_step_ahead
round(sqrt(sum((one_step_ahead - votes_vts_matrix[12, ])^2) / 51), digits = 3)
# global-alpha for comparison
set.seed(2024)
sparse_var_forecast <- fitVAR(votes_vts_matrix[1:11, ] - colMeans(votes_vts_matrix[1:11, ]), p = 2)
xhat =  colMeans(votes_vts_matrix[1:11, ]) + sparse_var_forecast$A[[1]] %*% (votes_vts_matrix[11, ] - colMeans(votes_vts_matrix[1:11, ])) + 
  sparse_var_forecast$A[[2]] %*% (votes_vts_matrix[10, ] - colMeans(votes_vts_matrix[1:11, ]))
round(sqrt(sum((xhat - votes_vts_matrix[12, ])^2 / 51)), digits = 3)
sparse_var_pe[i] = output[4, 3]
b = sparse_var_forecast$A[[1]]
b[abs(b) > 0.0] = 1
sum(rowSums(b))
b = sparse_var_forecast$A[[2]]
b[abs(b) > 0.0] = 1
sum(rowSums(b))
output[4, 1] = "Sparse VAR(1)"
summary(sparse_var_forecast)


global_alpha <- GNARfit(votes_vts_matrix - colMeans(votes_vts_matrix[1:11, ]), usa_net_GNAR, alphaOrder = 2, betaOrder = c(1, 0))
summary(global_alpha)
one_step_global_alpha = predict(global_alpha) + colMeans(votes_vts_matrix[1:11, ])
round(sqrt(sum((one_step_global_alpha - votes_vts_matrix[12, ])^2) / 51), digits = 3)
community_GNAR_AIC(c(residuals(global_alpha)), 3)
community_GNAR_BIC(c(residuals(global_alpha)), 4)


local_alpha <- GNARfit(votes_vts_matrix - colMeans(votes_vts_matrix[1:11, ]), usa_net_GNAR, alphaOrder = 2, betaOrder = c(1, 0), globalalpha = FALSE)
summary(local_alpha)
one_step_local_alpha = predict(local_alpha) + colMeans(votes_vts_matrix[1:11, ])
round(sqrt(sum((one_step_local_alpha - votes_vts_matrix[12, ])^2) / 51), digits = 3)
community_GNAR_AIC(c(residuals(local_alpha)), 104)
community_GNAR_BIC(c(residuals(local_alpha)), 104)
#############################################################
# Comparing the models after standardising the data
com_gnar_model <- community_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, c(2, 2, 2), list(c(0, 1), c(0, 1), c(1, 1)), 
                                     W_all, list(red_states, blue_states, swing_states))
summary(com_gnar_model)
# Fit a community-alpha GNAR([2, 2, 2], {(1, 0), (1, 0), (1, 0)}, 3) to three state-wise communities with interactions
com_gnar_model_with_interactions_std <- com_GNARfit_tv_weights(votes_vts_matrix_std[1:11, ], usa_net_GNAR, c(2, 2, 2), list(c(1, 0), c(1, 0), c(1, 1)), 
                                                               lapply(seq(1:11), function(x) {return(W_all)}), list(red_states, blue_states, swing_states), 
                                                               list(c(2, 3), c(1, 3), c(1, 2)))
summary(com_gnar_model_with_interactions_std)

one_step_ahead = community_GNARpred(usa_net_GNAR, W_all,  list(c(2, 2, 2), list(c(1, 0), c(1, 0), c(1, 1))), com_gnar_model_with_interactions_std$coefficients,
                                    votes_vts_matrix_std[1:11, ], list(red_states, blue_states, swing_states), list(c(2, 3), c(1, 3), c(1, 2)))
names(one_step_ahead) <- usa_net_data$Tag
round(sqrt(sum((one_step_ahead - votes_vts_matrix_std[12, ])^2)), digits = 3)
hist(one_step_ahead - votes_vts_matrix_std[12, ])
hist(com_gnar_model_with_interactions_std$residuals)
# global-alpha for comparison
summary(global_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, 2, c(1, 1), W_all))

round(sqrt(sum((votes_vts_matrix_std[12, ])^2)), digits = 3)


set.seed(2024)
sparse_var_forecast <- fitVAR(votes_vts_matrix[1:11, ], p = 2)
xhat = sparse_var_forecast$A[[1]] %*% votes_vts_matrix[11, ] + sparse_var_forecast$A[[2]] %*% votes_vts_matrix[10, ]
round(sqrt(sum((xhat - votes_vts_matrix[12, ])^2) / 51), digits = 3)
sparse_var_pe[i] = output[4, 3]
b = sparse_var_forecast$A[[2]]
b[abs(b) > 0.0] = 1
sum(rowSums(b))
output[4, 1] = "Sparse VAR(1)"
summary(sparse_var_forecast)
hist(sparse_var_forecast$residuals)
hist(xhat - votes_vts_matrix_std[12, ])



global_alpha <- GNARfit(votes_vts_matrix[1:11, ], usa_net_GNAR, alphaOrder = 2, betaOrder = c(0, 0))
summary(global_alpha)
round(sqrt(sum(predict(global_alpha) - votes_vts_matrix[12, ])^2), digits = 3)

local_alpha <- GNARfit(votes_vts_matrix_std[1:11, ], usa_net_GNAR, alphaOrder = 2, betaOrder = c(1, 1), globalalpha = FALSE)
summary(local_alpha)
round(sqrt(sum(predict(local_alpha) - votes_vts_matrix[12, ])^2), digits = 3)
hist(predict(local_alpha) - votes_vts_matrix_std[12, ])

#############################################
votes_vector_form <- rep(0, 12 * 51)

for (i in 1:12) {
  top = 1 + 51 * (i - 1)
  bottom = 51 + 51 * (i - 1)
  votes_vector_form[top:bottom] = votes_vts_matrix_std[i, ]
}
names(votes_vector_form) <- rep(usa_net_data$Tag, 12)
votes_vector_frame <- data.frame(votes_vector_form)
S1 = as.matrix(usa_net)
# synthetically connect Alaska and Hawaii
S1[2, 11] = 1
S1[11, 2] = 1
set.seed(2024)
car_ar_model <- ST.CARar(formula = votes_vector_form ~ 1, family = "gaussian", data = votes_vector_frame,
                         W = as.matrix(S1), AR = 2, burnin = 20000, n.sample = 220000, thin = 10) 
car_ar_model
