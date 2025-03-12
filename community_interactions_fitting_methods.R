# fitting methods for community-alpha with interactions GNAR models

get_interactions_time_step_inner <- function(stages_tensor, current_lag, p_lag, t, r_stages, W_covariate, covariate_data, interaction_data){
  d = ncol(covariate_data)
  # covariate_data = t(vapply(seq(1:nrow(NTS_data)), function(x) {return(covariate_group * NTS_data[x, ])}, rep(0, ncol(NTS_data))))
  # W_covariate =  t(vapply(seq(1:ncol(NTS_data)), function(x) {return(covariate_group[x] * weights_matrix[x, ])}, rep(0, ncol(weights_matrix))))
  # interaction_data = t(vapply(seq(1:nrow(NTS_data)), function(x) {return(interaction_group * NTS_data[x, ])}, rep(0, ncol(NTS_data))))
  group_component = get_time_step_inner(stages_tensor, current_lag, p_lag, t, r_stages, covariate_data, W_covariate)
  if (is.null(interaction_data) || r_stages[current_lag] == 0) {
    time_step_dmat <- group_component 
  } else {
  # print(interaction_columns)
  interaction_columns = lapply(seq(1:length(interaction_data)), function(x) { get_neighbourhood_regression(stages_tensor, current_lag, t, r_stages, 
                                                                                                             interaction_data[[x]], W_covariate)})
  interaction_column_merged = recursive_dmat_inner(interaction_columns)
  time_step_dmat = cbind(group_component, interaction_column_merged)
  }
  # print(time_step_dmat)
  return(time_step_dmat)
}

get_covariate_data_blocks <- function(nts_data, covariate_groups) {
  output = lapply(seq(1:length(covariate_groups)), function(x) {get_covariate_data_block(nts_data, covariate_groups[[x]])})
  return (output)
}

get_covariate_data_block <- function(nts_data, covariate_group) {
  covariate_data = t(vapply(seq(1:nrow(nts_data)), function(x) {return(covariate_group * nts_data[x, ])}, rep(0, ncol(nts_data))))
  return (covariate_data)
}

get_interactions_time_step_dmat <- function(stages_tensor, p_lag, t, r_stages, W_covariate, covariate_data, interaction_data) {
  dmat <- lapply(seq(1:p_lag), function(x) {get_interactions_time_step_inner(stages_tensor, x, p_lag, t, r_stages, W_covariate, covariate_data, interaction_data)})
  dmat_merged <- recursive_dmat_inner(dmat)
  # print(dmat_merged)
  return(dmat_merged)
}

get_covariate_weight_matrices <- function(weights_matrix, covariate_groups) {
  W_covariate_list <- lapply(seq(1:length(covariate_groups)), function(x) {get_covariate_weight_matrix(weights_matrix, covariate_groups[[x]])})
    return (W_covariate_list)
}

get_covariate_weight_matrix <- function(W, covariate_group) {
  output <- t(vapply(seq(1:ncol(W)), function(x) {return(covariate_group[x] * W[x, ])}, rep(0, ncol(W))))
  return (output)
}

build_covariate_design <- function(stages_tensor, p_lag, r_stages, W_list, covariate_data_blocks, current_covariate, interaction_covariates) { # nolint
  num_covariates = length(covariate_data_blocks)
  vts_data = covariate_data_blocks[[current_covariate]]
  n_steps = nrow(vts_data) - p_lag
  response_vec = lapply(seq((p_lag + 1), nrow(vts_data), 1), function(x) {return(cbind(vts_data[x, ]))})
  yvec = recursive_dmat(response_vec)
  Wc = W_list[[current_covariate]]
  if (missing(interaction_covariates)) {
    cat("Warning: no interactions specified returning no interactions.")
    # how to select all other covariates 
    # interaction_covariates = c(1:num_covariates)[c(1:num_covariates) != current_covariate]
    vts_interaction = NULL
  } else if (interaction_covariates[1] == 0) {
    vts_interaction = NULL
    col_names <- c('yvec', get_cov_block_coefficient_names(p_lag, r_stages, current_covariate))
  } else {
    vts_interaction = covariate_data_blocks[interaction_covariates]
    col_names <- c('yvec', get_cov_interaction_block_coefficient_names(p_lag, r_stages, current_covariate, interaction_covariates)) # nolint
  }
  design_mats = lapply(seq(1:n_steps), function(x) { 
    get_interactions_time_step_dmat(stages_tensor, p_lag, x + p_lag, 
    r_stages, Wc, vts_data, vts_interaction) }) # nolint
  dmat <- recursive_dmat(design_mats)
  covariate_model = cbind(yvec, dmat)
  colnames(covariate_model) <- col_names
  # print(col_names)
  # print(c('yvec', get_cov_interaction_block_coefficient_names(p_lag, r_stages, current_covariate, interaction_covariates))) # nolint: line_length_linter.
  return(data.frame(covariate_model, check.names = FALSE))
}

get_lag_interactions_block_names <- function(k, sk, c, interacting_covariates) {
  alpha = paste0("alpha.", as.character(k),'.', as.character(c))
  if (sk > 0) {
    betas = vapply(seq(1:sk), function(x) {paste0("beta.", as.character(k), ".", as.character(x), '.', 
                                                  as.character(c))}, c(""))
    gammas = vapply(seq(1:length(interacting_covariates)), function(x) {get_gamma_interaction_names(k, sk, c, interacting_covariates[[x]])}, 
                    rep("", sk))
    return(c(alpha, betas, gammas))
  } else {
    return(alpha)
  }
}

get_gamma_interaction_names <- function(k, sk, c, interacting_covariate) {
  gamma_interaction = vapply(seq(1:sk), function(x) {paste0("gamma.", as.character(k), ".", as.character(x), '.', 
                                                            as.character(c), ":", as.character(interacting_covariate))}, c(""))
  return (c(gamma_interaction))
}

get_cov_interaction_block_coefficient_names <- function(alpha_order, beta_order, c, interacting_covariates) {
  out = lapply(seq(1:alpha_order), function(x) {get_lag_interactions_block_names(x, beta_order[x], c, interacting_covariates)})
  return(merge_name_blocks(out))
}

get_com_interaction_gnar_names <- function(alpha_orders, beta_orders, num_covariates, interacting_covariates_list) {
  out = lapply(seq(1:num_covariates), function(x) {get_cov_interaction_block_coefficient_names(alpha_orders[x], beta_orders[[x]], x,
                                                                                               interacting_covariates_list[[x]])})
  return(merge_name_blocks(out))
}

community_interaction_gnar_fit <- function(vts, network, alpha_order, beta_order, weight_matrix, covariate_groups, interaction_groups) {
  beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0)
  stages_tensor = get_r_stages_adjacency_list(as.matrix(as.matrix(GNARtoigraph(network))), max(beta_depths))
  cov_data_blocks = get_covariate_data_blocks(vts, covariate_groups)
  W_list = get_covariate_weight_matrices(weight_matrix, covariate_groups)
  max_lag = max(alpha_order)
  d = ncol(vts)
  covariate_models <- lapply(seq(1:length(covariate_groups)), function(x) {build_covariate_design(stages_tensor, alpha_order[x],
                                                                                                  beta_order[[x]],
                                                                                                  W_list, cov_data_blocks, x, 
                                                                                                  interaction_groups[[x]])})
  if (length(unique(alpha_order)) != 1) {
  covariate_balanced_designs <- lapply(seq(1:length(covariate_models)), function(x) {grab_valid_observations(covariate_models[[x]],
                                                                                                             alpha_order[x], 
                                                                                                             max_lag, d)})
  full_model <- merge_covariate_designs(covariate_balanced_designs)
  } else {
    full_model <- merge_covariate_designs(covariate_models)
  }
  return(lm(yvec ~ . + 0, data = full_model))
}


community_interaction_GNAR_RIDGEfit <- function(vts, network, alpha_order, beta_order, weight_matrix, covariate_groups, interaction_groups) {
  beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0)
  stages_tensor = get_r_stages_adjacency_list(as.matrix(as.matrix(GNARtoigraph(network))), max(beta_depths))
  cov_data_blocks = get_covariate_data_blocks(vts, covariate_groups)
  W_list = get_covariate_weight_matrices(weight_matrix, covariate_groups)
  max_lag = max(alpha_order)
  d = ncol(vts)
  covariate_models <- lapply(seq(1:length(covariate_groups)), function(x) {build_covariate_design(stages_tensor, alpha_order[x],
                                                                                                  beta_order[[x]],
                                                                                                  W_list, cov_data_blocks, x,  # nolint: trailing_whitespace_linter.
                                                                                                  interaction_groups[[x]])})
  if (length(unique(alpha_order)) != 1) {
  covariate_balanced_designs <- lapply(seq(1:length(covariate_models)), function(x) {grab_valid_observations(covariate_models[[x]],
                                                                                                             alpha_order[x], 
                                                                                                             max_lag, d)}) # nolint: line_length_linter.
  full_model <- merge_covariate_designs(covariate_balanced_designs)
  } else {
    full_model <- merge_covariate_designs(covariate_models)
  }
  return(glmnet(full_model[, 2:ncol(full_model)], full_model[, 1], alpha = 0, lambda = 0, family = "gaussian", intercept = FALSE ))
}