# Fiiting community-alpha with interactions and time-varying weights
com_GNARfit_tv_weights <- function(vts, network, alpha_order, beta_order, weight_matrix_list, covariate_groups, interaction_groups) {
  beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0)
  if (max(beta_depths) > 0) {
  stages_tensor = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(network)), max(beta_depths))
  }
  cov_data_blocks = get_covariate_data_blocks(vts, covariate_groups)
  Wt_list = community_tv_weight_matrices(weight_matrix_list, covariate_groups)
  max_lag = max(alpha_order)
  d = ncol(vts)
  covariate_models <- lapply(seq(1:length(covariate_groups)), function(x) {build_covariate_design_tv_weights(stages_tensor, alpha_order[x],
                                                                                                  beta_order[[x]],
                                                                                                  Wt_list, cov_data_blocks, x, 
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

community_tv_weight_matrices <- function(weight_matrix_list, covariate_groups) {
  W_list_tv_weights = lapply(seq(1:length(weight_matrix_list)), function(x) {get_covariate_weight_matrices(weight_matrix_list[[x]], covariate_groups)})
  return(W_list_tv_weights)
}

build_linear_model_tv_weights <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix_list, ar = 'yes') {
  n_steps = nrow(NTS_data) - p_lag
  response_vec = lapply(seq((p_lag + 1), nrow(NTS_data), 1), function(x) {return(cbind(NTS_data[x, ]))})
  yvec = recursive_dmat(response_vec)
  if (ar == 'yes') {
    design_mats = lapply(seq(1:n_steps), function(x) { get_time_step_dmat_edt(stages_tensor, p_lag, x + p_lag, r_stages, NTS_data, weights_matrix_list[[x + p_lag]]) })
  } else {
    design_mats = lapply(seq(1:n_steps), function(x) { get_neighbourhood_regression(stages_tensor, p_lag, x + p_lag, 
                                                                                    r_stages, NTS_data, weights_matrix_list[[x + p_lag]]) })
  }
  dmat = recursive_dmat(design_mats)
  return (cbind(yvec, dmat))
}

build_covariate_design_tv_weights <- function(stages_tensor, p_lag, r_stages, Wt_list, covariate_data_blocks, current_covariate, interaction_covariates) {
  num_covariates = length(covariate_data_blocks)
  vts_data = covariate_data_blocks[[current_covariate]]
  n_steps = nrow(vts_data) - p_lag
  response_vec = lapply(seq((p_lag + 1), nrow(vts_data), 1), function(x) {return(cbind(vts_data[x, ]))})
  yvec = recursive_dmat(response_vec)
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
    col_names <- c('yvec', get_cov_interaction_block_coefficient_names(p_lag, r_stages, current_covariate, interaction_covariates))
  }
  design_mats = lapply(seq(1:n_steps), function(x) { get_interactions_time_step_dmat(stages_tensor, p_lag, x + p_lag, r_stages, Wt_list[[x + p_lag]][[current_covariate]], vts_data, vts_interaction) })
  dmat <- recursive_dmat(design_mats)
  covariate_model = cbind(yvec, dmat)
  colnames(covariate_model) <- col_names
  # print(col_names)
  # print(c('yvec', get_cov_interaction_block_coefficient_names(p_lag, r_stages, current_covariate, interaction_covariates)))
  return(data.frame(covariate_model, check.names = FALSE))
}

sim_time_step_community_tv_weights <- function(stages_tensor, weight_matrix_list, max_lag, alpha_vec, beta_vec, gamma_vec, covariate_group, current_cov, lags_data) {
  d = ncol(weight_matrix)
  if(max_lag > 1) {
    cov_lags = t(vapply(seq(1:max_lag), function(x) {return(covariate_group * lags_data[x, ])}, rep(0, d)))
  } else {
    cov_lags = matrix(covariate_group * lags_data[1, ], nrow = 1, ncol = d)
  }
  out = rep(0, d)
  for(k in 1:max_lag) {
    within_community_term = vapply(seq(1:length(beta_vec[[k]])), function(x) {(weight_matrix_list[[time_step]][[current_cov]] * stages_tensor[[r]]) %*% cov_lags[k, ]},
                                       rep(0, d))
    between_community_term = 
    out = out + alpha_vec[k] * cov_lags[k, ] + neighbourhood_term
  }
  return(as.numeric(out))
}
