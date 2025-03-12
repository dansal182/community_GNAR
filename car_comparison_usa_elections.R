library(CARBayesST)

votes_vector_form <- rep(0, 11 * 51)

for (i in 1:11) {
  top = 1 + 51 * (i - 1)
  bottom = 51 + 51 * (i - 1)
  votes_vector_form[top:bottom] = votes_vts_matrix[i, ] - colMeans(votes_vts_matrix[1:11, ])
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

one_step_car_ar = colMeans(votes_vts_matrix[1:11, ]) +  0.7069 * (votes_vts_matrix[11, ] - colMeans(votes_vts_matrix[1:11, ])) - 0.1077 * (votes_vts_matrix[10, ] - 
                                                                                                                          colMeans(votes_vts_matrix[1:11, ]))
round(sqrt(sum((one_step_car_ar - votes_vts_matrix[12, ])^2 / 51)), digits = 3)
corbit_plot(vts = votes_vts_matrix, net = usa_net_GNAR, max_lag = 8, max_stage = 11, viridis_color_option = "cividis")
hist(residuals(car_ar_model), breaks = 20)

round(sqrt(sum((48.8839 - votes_vts_matrix[12, ])^2) / 51), digits = 3)
round(sqrt(sum((votes_vts_matrix[11, ] - votes_vts_matrix[12, ])^2) / 51), digits = 3)
mean(votes_vts_matrix[11, ] - votes_vts_matrix[12, ])
mean(one_step_car_ar - votes_vts_matrix[12, ])
cor(one_step_car_ar, votes_vts_matrix[12, ])
cor(votes_vts_matrix[11, ], votes_vts_matrix[12, ])

#########################

com_model_two <- com_GNARfit_tv_weights(votes_vts_matrix[1:11, ], usa_net_GNAR, c(2, 2, 2), list(c(1, 1), c(1, 1), c(1, 1)), 
                       lapply(seq(1:11), function(x) {return(W_all)}), list(red_states, blue_states, swing_states), 
                       list(c(2), c(0), c(1, 2)))
summary(com_model_two)
hist(com_model_two$residuals)
one_step_ahead = community_GNARpred(usa_net_GNAR, W_all, list(c(2, 2, 2), list(c(1, 1), c(1, 1), c(1, 1))), com_model_two$coefficients, votes_vts_matrix[1:11, ],
                                    list(red_states, blue_states, swing_states), list(c(2), c(0), c(1, 2)))
names(one_step_ahead) <- usa_net_data$Tag
one_step_ahead
round(sqrt(sum((one_step_ahead - votes_vts_matrix[12, ])^2 / 51)), digits = 3)
pred_residuals = one_step_ahead - votes_vts_matrix[12, ]
mean(pred_residuals)
cor(one_step_ahead,  votes_vts_matrix[12, ])


##################################
com_model_two <- com_GNARfit_tv_weights(votes_vts_matrix[1:11, ] - colMeans(votes_vts_matrix[1:11, ]), usa_net_GNAR, c(2, 2, 2), list(c(1, 1), c(1, 1), c(1, 1)), 
                                        lapply(seq(1:12), function(x) {return(W_all)}), list(red_states, blue_states, swing_states), 
                                        list(c(2, 3), c(1, 3), c(1, 2)))
summary(com_model_two)
hist(com_model_two$residuals)
one_step_ahead = community_GNARpred(usa_net_GNAR, W_all, list(c(2, 2, 2), list(c(1, 1), c(1, 1), c(1, 1))), com_model_two$coefficients, 
                                    votes_vts_matrix[1:11, ] - colMeans(votes_vts_matrix[1:11, ]),
                                    list(red_states, blue_states, swing_states), list(c(2, 3), c(1, 3), c(1, 2))) + colMeans(votes_vts_matrix[1:11, ])
names(one_step_ahead) <- usa_net_data$Tag
one_step_ahead
round(sqrt(sum((one_step_ahead - votes_vts_matrix[12, ])^2) / 51), digits = 3)
pred_residuals = one_step_ahead - votes_vts_matrix[12, ]
mean(pred_residuals)
cor(one_step_ahead,  votes_vts_matrix[12, ])
hist(com_model_two$residuals, breaks = 20)

round(sqrt(sum((colMeans(votes_vts_matrix[1:11, ]) - votes_vts_matrix[12, ])^2) / 51), digits = 3)
round(sqrt(sum((votes_vts_matrix[11, ] - votes_vts_matrix[12, ])^2) / 51), digits = 3)
sqrt(sum((red_states * (votes_vts_matrix[11, ] - votes_vts_matrix[12, ])^2)))
sqrt(sum((blue_states * (votes_vts_matrix[11, ] - votes_vts_matrix[12, ])^2)))
sqrt(sum((swing_states * (votes_vts_matrix[11, ] - votes_vts_matrix[12, ])^2)))




com_model_two <- com_GNARfit_tv_weights(votes_vts_matrix - colMeans(votes_vts_matrix), usa_net_GNAR, c(2, 2, 2), list(c(1, 1), c(1, 1), c(1, 1)), 
                                        lapply(seq(1:12), function(x) {return(W_all)}), list(red_states, blue_states, swing_states), 
                                        list(c(2, 3), c(1, 3), c(1, 2)))
summary(com_model_two)
hist(com_model_two$residuals)
      
com_model_number_parameters = length(com_model_two$coefficients)
community_GNAR_AIC(com_model_two$residuals, com_model_number_parameters)
community_GNAR_BIC(com_model_two$residuals, com_model_number_parameters)
