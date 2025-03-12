# First we build the network

get_state_neighbour_tags <- function(neighbour_string) {
  state_list = strsplit(neighbour_string, split = ", ")
  num_neighbours = length(state_list[[1]])
  out = vapply(seq(1:num_neighbours), function(x) {get_state_neighbour_tag(state_list[[1]][x],  usa_net_data[, 1:2])}, c(""))
  return (out)
}

get_state_neighbour_tag <- function(state_list_i, tag_names) {
  out = ""
  if (!is.na(state_list_i)) {
    for (i in 1:nrow(tag_names)) {
      if (state_list_i == tag_names[i, 2]) {
        out = tag_names[i, 1]
        break
        } 
    }
  }
  return(out)
}

get_state_neighbour_ids <- function(tag_list) {
  out = vapply(seq(1:length(tag_list)), function(x) {get_state_neighbour_tag(tag_list[x], usa_net_data[, c(5, 1)])}, 0)
  return(out)
}

get_state_one_stage_neighbours <- function(neighbour_tags, tag_list) {
  out = rep(0, 51)
  for (i in 1:length(neighbour_tags)) {
    out = out + ifelse(tag_list == neighbour_tags[i], 1, 0)
  }
  return(out)
}


usa_net_data <- read.csv("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/data_sets/usa_net_names/Sheet 3-Neighbours.csv", 
                         header = TRUE)

get_state_neighbour_tag("Alabama", usa_net_data[, 1:2])
get_state_neighbour_tags(usa_net_data[2, 3])
edge_list = lapply(usa_net_data[, 3], get_state_neighbour_tags)

usa_net_data$state_id <- seq(1:51)

usa_adjacency_matrix = vapply(seq(1:51), function(x) { get_state_one_stage_neighbours(edge_list[[x]], usa_net_data[, 1])}, rep(0, 51))
colnames(usa_adjacency_matrix) <- usa_net_data[, 1]
rownames(usa_adjacency_matrix) <- usa_net_data[, 1]
usa_net <- graph_from_adjacency_matrix(usa_adjacency_matrix, mode = "undirected")
usa_net_GNAR <- igraphtoGNAR(usa_net)
plot(usa_net)

# We can start here now

library(corbit)
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/communal_gnar_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/communal_gnar_simulation_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/global_gnar_fitting.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/nacf_level_plots.R")
load("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/data_sets/USA_presidential_GNAR.RData")
save(usa_net_data, usa_net, usa_net_GNAR, rep_data, republican_votes_vts, votes_vts_matrix, red_states, blue_states,
     swing_states, W_all, W_red, W_blue, W_swing, 
     file = "/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/data_sets/USA_presidential_GNAR.RData")
y76_y20_election_data <- read.csv("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/data_sets/1976-2020-president.csv")
rep_data <- y76_y20_election_data[y76_y20_election_data[, 9] == "REPUBLICAN", ]
dem_data <-  y76_y20_election_data[y76_y20_election_data[, 9] == "DEMOCRAT", ]
rep_data$republican_percentage <- rep_data$candidatevotes / rep_data$totalvotes
republican_votes_vts <- lapply(seq(1:51), function(x) {rep_data[rep_data[, 3] == usa_net_data$Tag[x], c(1, 3, 11, 12, 16)]})
# Fix the odd observation (did not add the votes in 2016 for Maryland)
# We add the two 2016 Republican counts in Maryland
republican_votes_vts[[20]][11, 3] = republican_votes_vts[[20]][11, 3] + republican_votes_vts[[20]][12, 3]
republican_votes_vts[[20]][11, 5] = republican_votes_vts[[20]][11, 3] / republican_votes_vts[[20]][11, 4]
republican_votes_vts[[20]] = republican_votes_vts[[20]][c(1:11, 13), ]
# We build the network time series
votes_vts_matrix <- vapply(seq(1:51), function(x) {return(republican_votes_vts[[x]][, 5])}, rep(0, 12))
colnames(votes_vts_matrix) <- usa_net_data$Tag
rownames(votes_vts_matrix) <- republican_votes_vts[[20]]$year
votes_vts_matrix <- 100.0 * votes_vts_matrix
# vapply(seq(1:51), function(x) {length(republican_votes_vts[[x]][, 5])}, 0)
# We cluster the network
state_wins <- vapply(seq(1:51), function(x) {sum(ifelse(republican_votes_vts[[x]][, 5] > 0.5, 1, 0))}, 0)
names(state_wins) <- usa_net_data$Tag
max(state_wins)
min(state_wins)
community <- vapply(seq(1:51), function(x) {ifelse(state_wins[x] >= 9, "Red", ifelse(state_wins[x] > 3, "Swing", "Blue"))}, "")
names(community) <- usa_net_data$Tag

red_states <- vapply(seq(1:51), function(x) {ifelse(community[x] == "Red", 1, 0)}, 0)
blue_states <- vapply(seq(1:51), function(x) {ifelse(community[x] == "Blue", 1, 0)}, 0)
swing_states <- rep(1, 51) - red_states - blue_states


names(red_states) <- usa_net_data[, 1]
names(blue_states) <- usa_net_data[, 1]
names(swing_states) <- usa_net_data[, 1]

num_red_states = sum(red_states)
num_blue_states = sum(blue_states)
num_swing_states = sum(swing_states)

V(usa_net)$community <- community
V(usa_net)$color <- ifelse(V(usa_net)$community == "Red", "orange", ifelse(V(usa_net)$community == "Blue", "lightblue", "gray"))
V(usa_net)$frame.color <- ifelse(V(usa_net)$community == "Red", "orange", ifelse(V(usa_net)$community == "lightblue", "blue", "gray"))
plot(usa_net)
legend(x=1.0, y=1.0, c("Red States","Blue States", "Swing States"), pch=21, pt.bg=c("orange", "lightblue", "gray"), pt.cex=2, 
       cex=.8, bty="n", ncol=1)


# GNAR modelling
# Longest shortest path has a length of 11 and is between ME and NV, CA and WA
W_all = weights_matrix_cpp(usa_net_GNAR, 11)
colnames(W_all) <- usa_net_data$Tag
rownames(W_all) <- usa_net_data$Tag
corbit_plot(votes_vts_matrix, usa_net_GNAR, 5, 11, W_all)
corbit_plot(votes_vts_matrix, usa_net_GNAR, 5, 5, W_all, partial = "yes")
W_red = t(vapply(seq(1:51), function(x) {return(red_states[x] * W_all[x, ])}, rep(0, 51)))
W_blue = t(vapply(seq(1:51), function(x) {return(blue_states[x] * W_all[x, ])}, rep(0, 51)))
W_swing = t(vapply(seq(1:51), function(x) {return(swing_states[x] * W_all[x, ])}, rep(0, 51)))
rownames(W_red) <- usa_net_data$Tag
rownames(W_blue) <- usa_net_data$Tag
rownames(W_swing) <- usa_net_data$Tag
# The communal decomposition of the weight matrix
sum(W_all == W_red + W_blue + W_swing) == 51 * 51
# Communal decomposition for each vector time series
red_vts <- t(vapply(seq(1:12), function(x) {return(red_states * votes_vts_matrix[x, ])}, rep(0, 51)))
blue_vts <- t(vapply(seq(1:12), function(x) {return(blue_states * votes_vts_matrix[x, ])}, rep(0, 51)))
swing_vts <- t(vapply(seq(1:12), function(x) {return(swing_states * votes_vts_matrix[x, ])}, rep(0, 51)))
# Plot the time series by blocks
# Red States
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/red_states_vts.pdf")
par(mfrow = c(3, 5))
for(i in 1:51) {
  if (sum(red_vts[, i]) != 0) {
    plot(seq(1976, 2020, 4), red_vts[, i], xlab = "year", ylab = expression(X["i, t"]), main = usa_net_data$Tag[i], type = "l")
  }
}
par(mfrow = c(1, 1))
dev.off()
# Blue States
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/blue_states_vts.pdf")
par(mfrow = c(4, 5))
for(i in 1:51) {
  if (sum(blue_vts[, i]) != 0) {
    plot(seq(1976, 2020, 4), blue_vts[, i], xlab = "year", ylab = expression(X["i, t"]), main = usa_net_data$Tag[i], type = "l")
  }
}
# Swing States
par(mfrow = c(1, 1))
dev.off()
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/swing_states_vts.pdf")
par(mfrow = c(4, 5))
for(i in 1:51) {
  if (sum(swing_vts[, i]) != 0) {
    plot(seq(1976, 2020, 4), swing_vts[, i], xlab = "year", ylab = expression(X["i, t"]), main = usa_net_data$Tag[i], type = "l")
  }
}
par(mfrow = c(1, 1))
dev.off()
# Corbit plots by block
# Red States
corbit_plot(red_vts, usa_net_GNAR, 5, 5, W_red)
corbit_plot(red_vts, usa_net_GNAR, 5, 5, W_red, partial = "yes")
# Blue States
corbit_plot(blue_vts, usa_net_GNAR, 5, 5, W_blue)
corbit_plot(blue_vts, usa_net_GNAR, 5, 5, W_blue, partial = "yes")
# Swing States
corbit_plot(swing_vts, usa_net_GNAR, 5, 5, W_swing)
corbit_plot(swing_vts, usa_net_GNAR, 5, 5, W_swing, partial = "yes")
# Wagner Plots

wagner_plot(list(red_vts, blue_vts, swing_vts), list(usa_net_GNAR, usa_net_GNAR, usa_net_GNAR), 5, 5, list(W_red, W_blue, W_swing), c("Red States", "Blue States", "Swing States"), 
            same_net = "no")
par(mfrow=c(2, 2))
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/wagner_pnacf_both.pdf", height = 8.5, width = 14.00)
wagner_plot(list(red_vts, blue_vts, swing_vts), list(usa_net_GNAR, usa_net_GNAR, usa_net_GNAR), 5, 11, list(W_red, W_blue, W_swing), 
            c("Red States", "Blue States", "Swing States"), 
            same_net = "no", partial = "yes")
wagner_plot(list(red_vts_diff, blue_vts_diff, swing_vts_diff), list(usa_net_GNAR, usa_net_GNAR, usa_net_GNAR), 5, 11, list(W_red, W_blue, W_swing), c("Red States", "Blue States", "Swing States"), 
            same_net = "no", partial = "yes")
dev.off()
par(mfrow=c(1, 1))

nacf_level_plot(votes_vts_matrix, usa_net_GNAR, 2, 5, W_all)
pnacf_level_plot(red_vts, usa_net_GNAR, 5, 5, W_red)
pnacf_level_plot(blue_vts, usa_net_GNAR, 5, 5, W_blue)
pnacf_level_plot(swing_vts, usa_net_GNAR, 5, 5, W_swing)

nacf_level_plot(red_vts, usa_net_GNAR, 1, 5, W_red)
nacf_level_plot(blue_vts, usa_net_GNAR, 1, 5, W_blue)
nacf_level_plot(swing_vts, usa_net_GNAR, 1, 5, W_swing)
# Standardised data
votes_vts_matrix_std <- vapply(seq(1:51), function(x) {return((votes_vts_matrix[, x] - mean(votes_vts_matrix[, x])) / sd(votes_vts_matrix[, x]))}, rep(0, 12))

colnames(votes_vts_matrix_std) <- usa_net_data$Tag

com_gnar_model <- community_gnar_fit(votes_vts_matrix_std, usa_net_GNAR, c(2, 2, 2), list(c(1, 0), c(1, 0), c(1, 0)), 
                                     W_all, list(red_states, blue_states, swing_states))
summary(com_gnar_model)
cat(latex_table(round(rbind(com_gnar_model$coefficients, c(0.10411, 0.19966, 0.06649, 0.10805, 0.17842, 0.06309,0.08868, 0.16244, 0.06107)), 3),
                c("$\\hat{\\alpha}_{11}$", "$\\hat{\\beta}_{111}$", "$\\hat{\\alpha}_{21}$", "$\\hat{\\alpha}_{12}$",
                  "$\\hat{\\beta}_{112}$","$\\hat{\\alpha}_{22}$", "$\\hat{\\alpha}_{13}$", "$\\hat{\\beta}_{113}$", "$\\hat{\\alpha}_{23}$")))
hist(com_gnar_model$residuals, breaks = 20, xlab = expression(hat(u[t])), 
     main = "Residual histogram of GNAR(2, {[1, 0]}, 3) fit")
mean(com_gnar_model$residuals)
qqnorm(com_gnar_model$residuals)
quantile(com_gnar_model$residuals, probs = c(0.025, 0.975))
global_gnar_model <- global_gnar_fit(votes_vts_matrix, usa_net_GNAR, 2, c(1, 0), W_all)
summary(global_gnar_model)
