get_community_dmat <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, covariate_group, ar = "yes") {
  covariate_data = t(vapply(seq(1:nrow(NTS_data)), function(x) {return(covariate_group * NTS_data[x, ])}, rep(0, ncol(NTS_data))))
  # print(covariate_data)
  Wc = t(vapply(seq(1:ncol(NTS_data)), function(x) {return(covariate_group[x] * weights_matrix[x, ])}, rep(0, ncol(weights_matrix))))
  # print(Wc)
  covariate_design = build_linear_model_edt(stages_tensor, p_lag, r_stages, covariate_data, Wc, ar)
  return(covariate_design)
}

merge_covariate_designs <- function(covariate_models) {
  if (length(covariate_models) == 1) {
    return (covariate_models[[1]])
  } else {
    n_dmats = length(covariate_models)
    mid_point = n_dmats %/% 2
    left_split = merge_covariate_designs(covariate_models[1:mid_point])
    right_split = merge_covariate_designs(covariate_models[(mid_point + 1):n_dmats])
    return(merge_covariate_models(left_split, right_split))
  }
}

merge_covariate_models <- function(mod0, mod1) {
  yvec = mod0[, 1] + mod1[, 1]
  dmat = cbind(mod0[, 2:ncol(mod0)], mod1[, 2:ncol(mod1)])
  return(cbind(yvec, dmat))
}

merge_name_blocks <- function(names_list) {
  if (length(names_list) == 1) {
    return(names_list[[1]])
  } else {
    n_lists = length(names_list)
    mid_point = n_lists %/% 2
    left = merge_name_blocks(names_list[1:mid_point])
    right = merge_name_blocks(names_list[(1 + mid_point):n_lists])
    return(c(left, right))
  }
}

grab_valid_observations <- function(covariate_model, alpha_covariate, max_lag, vts_dimension) {
  if (max_lag == alpha_covariate) {
    return(covariate_model)
  } else {
    cutoff = (vts_dimension * (max_lag - alpha_covariate)) + 1
    return(covariate_model[cutoff:nrow(covariate_model), ])
  }
}

community_gnar_fit <- function(vts, network, alpha_order, beta_order, weight_matrix, covariate_groups, ar = 'yes') {
  beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0)
  stages_tensor = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(network)), max(max(beta_depths), 1))
  covariate_models <- lapply(seq(1:length(covariate_groups)), function(x) { grab_valid_observations(get_community_dmat(stages_tensor, alpha_order[x], beta_order[[x]], vts, 
                                                                                               weight_matrix, covariate_groups[[x]], ar), 
                                                                                               alpha_order[x], max(alpha_order), ncol(vts))})
  
  linear_model <- merge_covariate_designs(covariate_models)
  # print(linear_model)
  colnames(linear_model) <- c("yvec", get_com_gnar_names(alpha_order, beta_order, length(covariate_groups)))
  return (lm(yvec ~ . + 0, data = data.frame(linear_model)))
}

build_lm_data_frame <- function(covariate_model, alpha_order, beta_order, c) {
  colnames(covariate_model) <- c("yvec", get_cov_block_coefficient_names(alpha_order, beta_order, c))
  return(data.frame(covariate_model))
}

community_gnar_fit_list <- function(vts, network, alpha_order, beta_order, weight_matrix, covariate_groups, ar = 'yes') {
  beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0)
  stages_tensor = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(network)), max(beta_depths))
  covariate_models <- lapply(seq(1:length(covariate_groups)), function(x) { get_community_dmat(stages_tensor, 
                                                                                               alpha_order[x], beta_order[[x]], vts, 
                                                                                                weight_matrix, covariate_groups[[x]], ar)})
  
  linear_models <- lapply(seq(1:length(covariate_models)), function(x) {build_lm_data_frame(covariate_models[[x]], alpha_order[x],
                                                                                              beta_order[[x]], x)})
  # print(linear_models)
  linear_fits <- lapply(seq(1:length(linear_models)), function(x) {lm(yvec ~ . + 0, data = linear_models[[x]])})
  return(linear_fits)
}

get_cov_block_coefficient_names <- function(alpha_order, beta_order, c) {
  out = lapply(seq(1:alpha_order), function(x) {get_lag_block_names(x, beta_order[x], c)})
  return(merge_name_blocks(out))
}

get_lag_block_names <- function(k, sk, c) {
  alpha = paste0("alpha.", as.character(k),'.', as.character(c))
  if (sk > 0) {
    betas = vapply(seq(1:sk), function(x) {paste0("beta.", as.character(k), ".", as.character(x), '.', 
                                                  as.character(c))}, c(""))
    return(c(alpha, betas))
  } else {
    return(alpha)
  }
}

get_com_gnar_names <- function(alpha_orders, beta_orders, num_covariates) {
  out = lapply(seq(1:num_covariates), function(x) {get_cov_block_coefficient_names(alpha_orders[x], beta_orders[[x]], x)})
  return(merge_name_blocks(out))
}

