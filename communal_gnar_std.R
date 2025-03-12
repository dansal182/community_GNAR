votes_vts_matrix_diff <- votes_vts_matrix[2:12, ] - votes_vts_matrix[1:11, ]
votes_vts_matrix_std <- vapply(seq(1:51), function(x) {return((votes_vts_matrix[, x] - mean(votes_vts_matrix[, x])) / sd(votes_vts_matrix[, x]))}, rep(0, 11))
colnames(votes_vts_matrix_std) <- usa_net_data$Tag
red_vts_std <- t(vapply(seq(1:11), function(x) {return(red_states * votes_vts_matrix_std[x, ])}, rep(0, 51)))
blue_vts_std <- t(vapply(seq(1:11), function(x) {return(blue_states * votes_vts_matrix_std[x, ])}, rep(0, 51)))
swing_vts_std <- t(vapply(seq(1:11), function(x) {return(swing_states * votes_vts_matrix_std[x, ])}, rep(0, 51)))
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
# Wagner Plots
wagner_plot(list(red_vts_std, blue_vts_std, swing_vts_std), list(usa_net_GNAR), 8, 6, list(W_red, W_blue, W_swing), c("Red States", "Blue States", "Swing States"), 
            same_net = "yes")
wagner_plot(list(red_vts_std, blue_vts_std, swing_vts_std), list(usa_net_GNAR), 5, 5, list(W_red, W_blue, W_swing), c("Red States", "Blue States", "Swing States"), 
            same_net = "yes", partial = "yes")

com_gnar_model_std_vts <-  community_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, c(2, 2, 2), list(c(1, 0), c(0, 0), c(0, 1)), W_all, 
                                      list(red_states, blue_states, swing_states))
summary(com_gnar_model_std_vts)
hist(com_gnar_model_std_vts$residuals)
mean(com_gnar_model_std_vts$residuals)
qqnorm(com_gnar_model_std_vts$residuals)
quantile((com_gnar_model_std_vts$residuals + 0.02) / 0.76 , probs = c(0.025, 0.975))

global_gnar_model <- global_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, 2, c(1, 0), W_all)
summary(global_gnar_model)
hist(global_gnar_model$residuals)
mean(global_gnar_model$residuals)
quantile(global_gnar_model$residuals / 0.78, probs = c(0.025, 0.975))
