# One-step prediction for community alpha GNAR models

community_GNARpred <- function(network, weight_matrix, model_orders, theta_hat, design_data, covariate_groups, interaction_groups) {
    # get the number of communities
    C = length(covariate_groups)
    # compute the autoregressive and network regression orders 
    alpha_order = model_orders[[1]]
    beta_orders = model_orders[[2]]
    d = ncol(design_data)
    p_max = max(alpha_order)
    r_star = max(vapply(seq(1:C), function(x) {max(beta_orders[[x]])}, 0))
    # compute the r-adjacency matrices
    if (r_star > 0) {
        dummy_net = GNARtoigraph(network)
        adj_mat = as.matrix(dummy_net)
        stages_list <- get_r_stages_adjacency_list(as.matrix(adj_mat), r_star)
    }
    # compute the corresponding community blocks
    community_blocks = get_covariate_data_blocks(design_data, covariate_groups)
    # compute the community-wise weight matrices
    Wt_list = community_tv_weight_matrices(list(weight_matrix), covariate_groups) 
    # community-wise linear model for one-step prediction (i.e., one time step design construction)
    community_designs <- lapply(seq(1:C), function(x) {build_community_design_time_step(p_max + 1, stages_list, alpha_order[x], p_max,
                                                                                                  beta_orders[[x]],
                                                                                                  Wt_list, community_blocks, x, 
                                                                                                  interaction_groups[[x]])})
    # merge the communities into the complete model
    time_step_design = merge_community_time_step_designs(community_designs)
    # print(time_step_design)
    # substract intercept term from columns
    #time_step_design = time_step_design - theta_hat[1]
    # include intercept term
    #time_step_design = cbind(rep(1, d), time_step_design)
    # compute the estimated mean at the next time step
    if (ncol(time_step_design) != length(theta_hat)) {
        warning("Model order incompatible with estimated coefficients.")
        break
    } else {
        forecast_value = c(time_step_design %*% theta_hat)
    }                                                                               
    names(forecast_value) <- vapply(seq(1:d), function(x) {paste0("X", as.character(x))}, "")
    return (forecast_value)
}

community_GNAR_AIC <- function(model_residuals, number_of_parameters) {
    squared_error = sum(model_residuals^2) / length(model_residuals)
    out = log(squared_error) + 2 * number_of_parameters / length(model_residuals)
    return (out)
}

community_GNAR_BIC <- function(model_residuals, number_of_parameters) {
    squared_error = sum(model_residuals^2)
    out = (squared_error + log(length(model_residuals)) * number_of_parameters) / length(model_residuals)
    return (out)
}