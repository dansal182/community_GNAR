# Difference the network time series
votes_vts_matrix_diff <- votes_vts_matrix[2:12, ] - votes_vts_matrix[1:11, ]
colnames(votes_vts_matrix_diff) <- usa_net_data$Tag
red_vts_diff <- t(vapply(seq(1:11), function(x) {return(red_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))
blue_vts_diff <- t(vapply(seq(1:11), function(x) {return(blue_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))
swing_vts_diff <- t(vapply(seq(1:11), function(x) {return(swing_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))
# standardised series
votes_vts_matrix_diff_std <- vapply(seq(1:51), function(x) {return((votes_vts_matrix_diff[, x] - mean(votes_vts_matrix_diff[, x])) / sd(votes_vts_matrix_diff[, x]))}, 
                               rep(0, 11))
colnames(votes_vts_matrix_diff_std) <- usa_net_data$Tag
red_vts_std <- t(vapply(seq(1:11), function(x) {return(red_states * votes_vts_matrix_diff_std[x, ])}, rep(0, 51)))
blue_vts_std <- t(vapply(seq(1:11), function(x) {return(blue_states * votes_vts_matrix_diff_std[x, ])}, rep(0, 51)))
swing_vts_std <- t(vapply(seq(1:11), function(x) {return(swing_states * votes_vts_matrix_diff_std[x, ])}, rep(0, 51)))
# Corbit plots by block
# Red States
corbit_plot(red_vts_std, usa_net_GNAR, 5, 5, W_red)
corbit_plot(red_vts_diff, usa_net_GNAR, 5, 5, W_red, partial = "yes")
# Blue States
corbit_plot(blue_vts_std, usa_net_GNAR, 5, 5, W_blue)
corbit_plot(blue_vts_diff, usa_net_GNAR, 5, 5, W_blue, partial = "yes")
# Swing States
corbit_plot(swing_vts_std, usa_net_GNAR, 5, 5, W_swing)
corbit_plot(swing_vts_diff, usa_net_GNAR, 5, 5, W_swing, partial = "yes")

# Wagner Plots
wagner_plot(list(red_vts_std, blue_vts_std, swing_vts_std), list(usa_net_GNAR), 8, 6, list(W_red, W_blue, W_swing), c("Red States", "Blue States", "Swing States"), 
            same_net = "yes")
wagner_plot(list(red_vts_std, blue_vts_std, swing_vts_std), list(usa_net_GNAR), 5, 3, list(W_red, W_blue, W_swing), c("Red States", "Blue States", "Swing States"), 
            same_net = "yes", partial = "yes")
# Wagner Plots
wagner_plot(list(red_vts_diff, blue_vts_std, swing_vts_std), list(usa_net_GNAR), 8, 6, list(W_red, W_blue, W_swing), c("Red States", "Blue States", "Swing States"), 
            same_net = "yes")
wagner_plot(list(red_vts_diff, blue_vts_diff, swing_vts_diff), list(usa_net_GNAR, usa_net_GNAR, usa_net_GNAR), 5, 5, list(W_red, W_blue, W_swing), c("Red States", "Blue States", "Swing States"), 
            same_net = "no", partial = "yes")
pnacf_level_plot(red_vts_diff, usa_net_GNAR, 2, 11, W_red)
pnacf_level_plot(blue_vts_diff, usa_net_GNAR, 2,11, W_blue)
pnacf_level_plot(swing_vts_diff, usa_net_GNAR, 2, 11, W_swing)
# Fitting Communal-alpha GNAR models with different orders
com_gnar_model_diff <-  community_gnar_fit(votes_vts_matrix_diff_std, usa_net_GNAR, c(3, 3, 3), list(c(0, 2, 0), c(0, 3, 0), c(0, 2, 0)), W_all, 
                                              list(red_states, blue_states, swing_states))
summary(com_gnar_model_diff)
# Prints the model summary in Latex format
cat(latex_table(round(rbind(com_gnar_model_diff$coefficients, c(0.08584, 0.06390, 0.08241, 0.07572, 0.05457, 0.06899,  0.07993, 0.05383, 0.07229)), 3),
                c("$\\hat{\\alpha}_{21}$", "\\hat{\\alpha}_{21}$", "$\\hat{\\alpha}_{31}$", "$\\hat{\\alpha}_{12}$",
                                                                        "$\\hat{\\alpha}_{22}$", "$\\hat{\\alpha}_{32}$",
                                                             "$\\hat{\\alpha}_{13}$", "$\\hat{\\alpha}_{23}$", "$\\hat{\\alpha}_{33}$")))

hist(com_gnar_model_diff$residuals, breaks = 20,  xlab = expression(hat(u[t])), 
     main = "Residual histogram of GNAR(3, {[0, 0, 0]}, 3) fit")
mean(com_gnar_model_diff$residuals)
qqnorm(com_gnar_model_diff$residuals)
quantile((com_gnar_model_diff$residuals + 0.02) / 0.76 , probs = c(0.025, 0.975))
sd(com_gnar_model_diff$residuals[1:200])
sd(com_gnar_model_diff$residuals[200:400])

# Parallel estimator (stacks the block-models)
model_two <- community_gnar_fit_list(votes_vts_matrix_diff_std, usa_net_GNAR, c(1, 3, 2), list(c(1), c(0, 0, 0), c(1, 0)), W_all, 
                                     list(red_states, blue_states, swing_states))
summary(model_two[[1]])
summary(model_two[[2]])
summary(model_two[[3]])
# Verify the estimated model is stationary and some model validaiton aids
sum(vapply(seq(1:3), function(x) {ifelse(sum(abs(model_two[[x]]$coefficients)) < 1.0, 1, 0)}, 0))


hist(model_two[[1]]$residuals + model_two[[2]]$residuals + model_two[[3]]$residuals)
global_gnar_model <- global_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, 2, c(1, 0), W_all)
summary(global_gnar_model)
hist(global_gnar_model$residuals)
mean(global_gnar_model$residuals)
quantile(global_gnar_model$residuals / 0.78, probs = c(0.025, 0.975))

# Compare with a global-alpha GNAR with model order (1, [2])
summary(community_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, c(1), list(c(2)), W_all, list(rep(1, 51))))

com_gnar_model_std <-  community_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, c(2, 3, 2), list(c(0, 1), c(0, 0, 0), c(1, 0)), W_all, 
                                              list(red_states, blue_states, swing_states))
summary(com_gnar_model_std)
# Verify the estimated model is stationary and some model validaiton aids
sum(vapply(seq(1:3), function(x) {ifelse(sum(abs(model_two[[x]]$coefficients)) < 1.0, 1, 0)}, 0))
hist(com_gnar_model_std$residuals)
mean(com_gnar_model_std_vts$residuals)
qqnorm(com_gnar_model_std_vts$residuals)
quantile((com_gnar_model_std_vts$residuals + 0.02) / 0.76 , probs = c(0.025, 0.975))

####################################################
# Plot the time series by blocks
# Red States
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/red_states_std_vts.pdf")
par(mfrow = c(3, 5))
for(i in 1:51) {
  if (sum(red_vts[, i]) != 0) {
    plot(seq(1980, 2020, 4), red_vts_std[, i], xlab = "year", ylab = expression(X["i, t"]), main = usa_net_data$Tag[i], type = "l")
  }
}
par(mfrow = c(1, 1))
dev.off()
# Blue States
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/blue_states_std_vts.pdf")
par(mfrow = c(4, 5))
for(i in 1:51) {
  if (sum(blue_vts[, i]) != 0) {
    plot(seq(1980, 2020, 4), blue_vts_std[, i], xlab = "year", ylab = expression(X["i, t"]), main = usa_net_data$Tag[i], type = "l")
  }
}
# Swing States
par(mfrow = c(1, 1))
dev.off()
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/swing_states_std_vts.pdf")
par(mfrow = c(4, 5))
for(i in 1:51) {
  if (sum(swing_vts[, i]) != 0) {
    plot(seq(1980, 2020, 4), swing_vts_std[, i], xlab = "year", ylab = expression(X["i, t"]), main = usa_net_data$Tag[i], type = "l")
  }
}
par(mfrow = c(1, 1))
dev.off()
# All in one plot
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/all_states_std_one_diff_vts.pdf")
par(mfrow = c(3, 5))
for(i in 1:51) {
  if (sum(red_vts[, i]) != 0) {
    plot(seq(1980, 2020, 4), red_vts_std[, i], xlab = "year", ylab = expression(X["i, t"]), main = usa_net_data$Tag[i], type = "l", col = "orange")
  } else if(sum(blue_vts[, i]) != 0) {
    plot(seq(1980, 2020, 4), blue_vts_std[, i], xlab = "year", ylab = expression(X["i, t"]), main = usa_net_data$Tag[i], type = "l", col = "lightblue")
  } else {
    plot(seq(1980, 2020, 4), swing_vts_std[, i], xlab = "year", ylab = expression(X["i, t"]), main = usa_net_data$Tag[i], type = "l", col = "black")
  }
}
par(mfrow = c(1, 1))
dev.off()


cat(latex_table(round(rbind(com_gnar_model_diff$coefficients, c(0.08584, 0.06390, 0.08241, 0.07572, 0.05457, 0.06899,  0.07993, 0.05383, 0.07229)), 3),
                c("$\\hat{\\alpha}_{21}$", "\\hat{\\alpha}_{21}$", "$\\hat{\\alpha}_{31}$", "$\\hat{\\alpha}_{12}$",
                  "$\\hat{\\alpha}_{22}$", "$\\hat{\\alpha}_{32}$",
                  "$\\hat{\\alpha}_{13}$", "$\\hat{\\alpha}_{23}$", "$\\hat{\\alpha}_{33}$")))
summary(model_two[[1]])
summary(model_two[[2]])
summary(model_two[[3]])

model_two <- community_gnar_fit_list(votes_vts_matrix_diff, usa_net_GNAR, c(1, 3, 2), list(c(1), c(0, 0, 0), c(1, 0)), W_all, 
                                     list(red_states, blue_states, swing_states))