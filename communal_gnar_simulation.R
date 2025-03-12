# Simulation for communal alpha GNAR
library(GNAR) # Cobrit and Wagner plots
# community-alpha GNAR fitting code
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_gnar_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/communal_gnar_simulation_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/global_gnar_fitting.R")
# PNACF level plots
source("/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/nacf_level_plots.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_interactions_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/corbit_paper/R_scripts/R_scripts_revised/R-Corbit_Covid.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/latex_output_and_utilities/latex_outputs_utilities.R")
source("/Users/danielsalnikov/Documents/code_scripts/fastCorbit/R/fast_corbit_plot_cpp.r")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit/R/correlation_structure_plots.R")


W = weights_matrix_cpp(fiveNet, 3)
sim0 <- gnar_sim_community(200, fiveNet, W, c(1, 1), c(1, 1), list(c(0.15), c(0.30)), list(list(c(0.30)), list(c(0.12))), list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)))
W1 <- t(vapply(seq(1:ncol(sim0)), function(x) {return(c(0, 1, 1, 1, 0)[x] * W[x, ])}, rep(0, ncol(W))))
W1 = t(t(W1) * c(0, 1, 1, 1, 0))
W2 <- t(vapply(seq(1:ncol(sim0)), function(x) {return(c(1, 0, 0, 0, 1)[x] * W[x, ])}, rep(0, ncol(W))))
W2 = t(t(W2) * c(1, 0, 0, 0, 1))
sum(W1 + W2 == W) == 3 * 3
fiveNet_stages <- get_r_stages_adjacency_list(as.matrix(as.matrix(GNARtoigraph(fiveNet))), 3)
# get_interactions_time_step_dmat(fiveNet_stages, 1, 2, 4, c(3, 3), sim0, W, c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1))
#par(mfrow=c(2, 3))
#for (l in 1:5) {
#  plot(sim0[, l], ylab = expression(X["i, t"]), xlab = "t", main = paste0("Time Series: i = ", as.character(l)), type = "l")
#}
#par(mfrow=c(1, 1))

gnar_model <- community_gnar_fit(sim0, fiveNet, rep(1, 2), list(c(1), c(1)), W, list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)))
summary(gnar_model)     
my_seed = 1983 + c(0:9)
params = matrix(rep(0, 60), nrow = 6, ncol = 10)

for (i in 1:10) {
  set.seed(2024)
  sim1 <- gnar_sim_community(1000, fiveNet, W, c(1, 2), c(1, 1), list(c(0.23), c(0.20, 0.18)), list(list(c(0.47)), list(c(0.30), c(0.27))), 
                           list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)))
  community_one_sim1 <- t(vapply(seq(1:nrow(sim1)), function(x) {return(c(0, 1, 1, 1, 0) * sim1[x, ])}, rep(0, ncol(sim1))))
  community_two_sim1 <- t(vapply(seq(1:nrow(sim1)), function(x) {return(c(1, 0, 0, 0, 1) * sim1[x, ])}, rep(0, ncol(sim1))))
  cross_correlation_plot(1, sim1)
  cross_correlation_plot(2, sim1)
  cross_correlation_plot(4, sim1)
  cross_correlation_plot(8, sim1)
  cross_correlation_plot(16, sim1)
# par(mfrow=c(2, 3))
# for (l in 1:5) {
#  plot(sim1[, l], ylab = expression(X["i, t"]), xlab = "t", main = paste0("Time Series: i = ", as.character(l)), type = "l")
# }
#par(mfrow=c(1, 1))
s0 <- Sys.time()
gnar_com_diff_order <- community_gnar_fit(sim1, fiveNet, c(1, 2), list(c(1), c(1, 1)), W, list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)))
summary(gnar_com_diff_order)
#sum(abs(gnar_com_diff_order$coefficients[1:2]))
#sum(abs(gnar_com_diff_order$coefficients[3:6]))
#hist(gnar_com_diff_order$residuals)
#qqnorm(gnar_com_diff_order$residuals)
s1 <- Sys.time()
print(round(s1 - s0, 2))
params[, i] = gnar_com_diff_order$coefficients

file_name = paste0("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/R-Corbit_sim", as.character(i), ".pdf")
pdf(file_name, width = 14.00, height = 8.50)
r_corbit_plot_cpp(list(community_one_sim1, community_two_sim1), list(fiveNet, fiveNet), 6, 3, list(W, W), frame_names = c("Community One", "Community Two"), 
            same_net = 'no', partial = 'yes', viridis_color_option = 'turbo')
dev.off()
}

params_std <- vapply(seq(1:6), function(x) { sqrt(1 / 6 * sum((params[x, ] - rowMeans(params)[x])^2))}, 0)
cat(latex_table(cbind(round(rowMeans(params), 3), round(params_std, 3), c(0.27, 0.18, 0.25, 0.30, 0.12, 0.20)), c("Estimate", "Std. Error", "True Value")))

# It appears that the non-parallelised works faster.

#community_one_sim0 <- t(vapply(seq(1:nrow(sim0)), function(x) {return(c(0, 1, 1, 1, 0) * sim0[x, ])}, rep(0, ncol(sim0))))
#community_two_sim0 <- t(vapply(seq(1:nrow(sim0)), function(x) {return(c(1, 0, 0, 0, 1) * sim0[x, ])}, rep(0, ncol(sim0))))



#wagner_plot(list(community_one_sim0, community_two_sim0), list(fiveNet), 5, 3, list(W1, W2), c("Community One", "Community Two"), 
#            same_net = 'yes')

#wagner_plot(list(community_one_sim0, community_two_sim0), list(fiveNet), 5, 3, list(W1, W2), c("Community One", "Community Two"), 
#            same_net = 'yes', partial = 'yes')

#corbit_plot(community_one_sim0, fiveNet, 10, 3, W1)
#corbit_plot(community_one_sim0, fiveNet, 10, 3, W1, partial = 'yes')
#corbit_plot(community_two_sim0, fiveNet, 10, 3, W2)
#corbit_plot(community_two_sim0, fiveNet, 10, 3, W2, partial = 'yes')


#community_one_sim1 <- t(vapply(seq(1:nrow(sim1)), function(x) {return(c(0, 1, 1, 1, 0) * sim1[x, ])}, rep(0, ncol(sim1))))
#community_two_sim1 <- t(vapply(seq(1:nrow(sim1)), function(x) {return(c(1, 0, 0, 0, 1) * sim1[x, ])}, rep(0, ncol(sim1))))
#pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/wagner_nacf.pdf")
#wagner_plot(list(community_one_sim1, community_two_sim1), list(fiveNet, fiveNet), 10, 3, list(W1, W2), c("Community One", "Community Two"), 
#            same_net = 'no')
#dev.off()
#pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/wagner_pnacf.pdf", width = 14.00, height = 8.50)
#wagner_plot(list(community_one_sim1, community_two_sim1), list(fiveNet, fiveNet), 8, 3, list(W1, W2), c("Community One", "Community Two"), 
#            same_net = 'no', partial = 'yes')
#pnacf_level_plot(community_one_sim1, fiveNet, 2, 3, W1)
#pnacf_level_plot(community_two_sim1, fiveNet, 2, 3, W2)
#dev.off()
#pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/corbit_com_one_nacf.pdf")
#corbit_plot(community_one_sim1, fiveNet, 10, 3, W1)
#dev.off()
#pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/corbit_com_one_pnacf.pdf")
#corbit_plot(community_one_sim1, fiveNet, 10, 3, W1, partial = 'yes')
#dev.off()
#pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/corbit_com_two_nacf.pdf")
#corbit_plot(community_two_sim1, fiveNet, 10, 3, W2)
#dev.off()
#pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/corbit_com_two_pnacf.pdf")
#corbit_plot(community_two_sim1, fiveNet, 10, 3, W2, partial = 'yes')
#dev.off()
#pnacf_level_plot(community_two_sim1, fiveNet, 1, 3, W)#

#edge_list <- data.frame(from = c("1", "1", "2", "2", "3"), to = c("5", "4", "4", "3", "4"))
#node_list <- data.frame(names = c("1", "2", "3", "4", "5"), community = c("two", "one", "one", "one", "two"))
#com_net <- graph_from_data_frame(d = edge_list, vertices = node_list, directed = FALSE)
#plot(com_net)
#V(com_net)$color <- ifelse(V(com_net)$community == "one", "lightblue", "orange")
#V(com_net)$frame.color <- ifelse(V(com_net)$community == "one", "lightblue", "orange")
# pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/images/fiveNet_communities.pdf")
#plot(com_net)
#legend(x='bottomleft', c("Community One","Community Two"), pch=21, pt.bg=c("lightblue", "orange"), pt.cex=2, 
#       cex=.8, bty="n", ncol=1)
# dev.off()
#plot(fiveNet, mark.groups = list(c(2, 3, 4), c(1, 5)))
