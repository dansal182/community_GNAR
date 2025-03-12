get_neighbourhood_regression <- function(stages_tensor, current_lag, t, r_stages, NTS_data, weights_matrix) {
  max_stage = r_stages[current_lag]
  Z = (weights_matrix * stages_tensor[[1]]) %*% NTS_data[t - current_lag, ]
  if (max_stage > 1) {
    for (stage in 2:max_stage) {
        Z_r = (weights_matrix * stages_tensor[[stage]]) %*% NTS_data[t - current_lag, ]
        Z = cbind(Z, Z_r)
    }
  }
  return(Z)
}


get_time_step_dmat <- function(stages_tensor, p_lag, t, r_stages, NTS_data, weights_matrix) {
  for (current_lag in 1:p_lag) {
    xvec = NTS_data[t - current_lag, ]
    Z = get_neighbourhood_regression(stages_tensor, current_lag, t, r_stages, NTS_data, weights_matrix)
    if (current_lag == 1) {
      dmat = cbind(xvec, Z)
    } else {
      aux = cbind(xvec, Z)
      dmat = cbind(dmat, aux)
    }
  }
  return (as.matrix(dmat))
}


recursive_dmat <- function(dmat_list) {
  if (length(dmat_list) == 1) {
    return (dmat_list[[1]])
  } else {
    n_dmats = length(dmat_list)
    mid_point = n_dmats %/% 2
    left_split = recursive_dmat(dmat_list[1:mid_point])
    right_split = recursive_dmat(dmat_list[(mid_point + 1):n_dmats])
    return (rbind(left_split, right_split))
  }
}


build_linear_model <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar = 'yes') {
  n_steps = nrow(NTS_data) - p_lag
  response_vec = lapply(seq(1:n_steps), function(x) {return(cbind(NTS_data[x + 1, ]))})
  yvec = recursive_dmat(response_vec)
  if (ar == 'yes') {
    design_mats = lapply(seq(1:n_steps), function(x) { get_time_step_dmat(stages_tensor, p_lag, x + p_lag, r_stages, NTS_data, weights_matrix) })
  } else {
    design_mats = lapply(seq(1:n_steps), function(x) { get_neighbourhood_regression(stages_tensor, p_lag, x + p_lag, 
                                                                                             r_stages, NTS_data, weights_matrix) })
  }
  dmat = recursive_dmat(design_mats)
  return (cbind(yvec, dmat))
}


fit_gnar <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar = 'yes') {
  model_data = build_linear_model(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar)
  gnar_fit <- lm(model_data[, 1] ~ model_data[, 2:ncol(model_data)] + 0)
  # print(summary(gnar_fit))
  return (gnar_fit)
}


get_community_dmat <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar, covariate_group) {
  covariate_data = t(vapply(seq(1:nrow(NTS_data)), function(x) {return(covariate_group * NTS_data[x, ])}, rep(0, ncol(NTS_data))))
  Wc = t(vapply(seq(1:ncol(NTS_data)), function(x) {return(covariate_group[x] * weights_matrix[x, ])}, rep(0, ncol(weights_matrix))))
  covariate_design = build_linear_model(stages_tensor, p_lag, r_stages, covariate_data, Wc, ar)
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

community_gnar_fit <- function(vts, network, alpha_order, beta_order, weight_matrix, covariate_groups, ar = 'yes') {
  beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0)
  stages_tensor = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(network)), max(beta_depths))
  covariate_models <- lapply(seq(1:length(covariate_groups)), function(x) { get_community_dmat(stages_tensor, alpha_order[[x]], beta_order[[x]], vts, 
                                                                                                      weight_matrix, ar, covariate_groups[[x]])})
  linear_model <- merge_covariate_designs(covariate_models)
  colnames(linear_model) <- c("yvec", get_com_gnar_names(alpha_order, beta_order, length(covariate_groups)))
  return (lm(yvec ~ . + 0, data = data.frame(linear_model)))
}


gnar_sim_community <- function(nsims, network, weight_matrix, max_lag, r_stages, alpha_vec, beta_vec, covariate_groups, sigma=1) {
  d = ncol(weight_matrix)
  out = matrix(rep(0, d * (nsims + max_lag)), nrow = nsims + max_lag, ncol = d)
  stages_tensor = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(network)), max(r_stages))
  lag_data = matrix(rep(0, d * max_lag), nrow = max_lag, ncol = d)
  for (k in 1:max_lag) {
    out[k, ] = rnorm(d, 0, 1)
  }
  for (t in max_lag:(nsims + max_lag - 1)) {
      for (p in 1:max_lag) {
        lag_data[p, ] = out[(t - p + 1), ]
      }
      X = sim_all_communities(stages_tensor, weight_matrix, max_lag, alpha_vec, beta_vec, covariate_groups, lag_data)
      out[t + 1, ] = X + rnorm(d, 0, sigma)
  }
  return(out[(max_lag+1):(nsims + max_lag), ])
  }


sim_time_step_community <- function(stages_tensor, weight_matrix, max_lag, alpha_vec, beta_vec, covariate_group, lags_data) {
  d = ncol(weight_matrix)
  Wc = t(vapply(seq(1:d), function(x) {return(covariate_group[x] * weight_matrix[x, ])}, rep(0, d)))
  if(max_lag > 1) {
    cov_lags = t(vapply(seq(1:max_lag), function(x) {return(covariate_group * lags_data[x, ])}, rep(0, d)))
  } else {
    cov_lags = matrix(covariate_group * lags_data[1, ], nrow = 1, ncol = d)
  }
  out = rep(0, d)
  for(k in 1:max_lag) {
    neighbourhood_term = rep(0, d)
    for (r in 1:length(beta_vec[[k]])) {
      neighbourhood_term = neighbourhood_term + beta_vec[[k]][r] * (Wc * stages_tensor[[r]]) %*% cov_lags[k, ]
    }
    out = out + alpha_vec[k] * cov_lags[k, ] + neighbourhood_term
  }
  return(as.numeric(out))
}


sim_all_communities <- function(stages_tensor, weight_matrix, max_lag, alpha_vec, beta_vec, covariate_groups, lags_data) {
  d = ncol(weight_matrix)
  Xc = vapply(seq(1:length(covariate_groups)), function(x) {return (sim_time_step_community(stages_tensor, weight_matrix, max_lag,
                                                                                            alpha_vec[[x]], beta_vec[[x]], covariate_groups[[x]], 
                                                                                            lags_data))},
         rep(0, d))
  X = rowSums(Xc)
  return(X)
}


get_cov_block_coefficient_names <- function(alpha_order, beta_order, c) {
 out = lapply(seq(1:alpha_order), function(x) {get_lag_block_names(x, beta_order[x], c)})
 return(merge_name_blocks(out))
}

get_lag_block_names <- function(k, sk, c) {
  alpha = paste0("alpha.", as.character(k), as.character(c))
  betas = vapply(seq(1:sk), function(x) {paste0("beta.", as.character(k), as.character(x), as.character(c))}, c(""))
  return(c(alpha, betas))
}

get_com_gnar_names <- function(alpha_orders, beta_orders, num_covariates) {
  out = lapply(seq(1:num_covariates), function(x) {get_cov_block_coefficient_names(alpha_orders[x], beta_orders[[x]], x)})
  return(merge_name_blocks(out))
}

#######################
# Old slow methods
get_time_dmat <- function(stages_tensor, p_lag, t, r_stages, NTS_data, weights_tensor, ar) {
  if (ar == 'yes') {
    current_lag = 1
    xvec = cbind(NTS_data[t - current_lag, ])
    Z = (weights_tensor * stages_tensor[[1]]) %*% NTS_data[t - current_lag, ]
    if (r_stages[current_lag] > 1) {
      for (r in 2:r_stages[current_lag]) {
        Z_r = (weights_tensor * stages_tensor[[r]]) %*% NTS_data[t - current_lag, ]
        Z = cbind(Z, Z_r)
      }
    }
    dmat = cbind(xvec, Z)
    if (p_lag > 1) {
      for (current_lag in 2:p_lag) {
        xvec = cbind(NTS_data[t - current_lag, ])
        Z = (weights_tensor * stages_tensor[[1]]) %*% NTS_data[t - current_lag, ]
        if (r_stages[current_lag] > 1) {
          for (r in 2:r_stages[current_lag]) {
            Z_r = (weights_tensor * stages_tensor[[r]]) %*% NTS_data[t - current_lag, ]
            Z = cbind(Z, Z_r)
          }
        }
        aux = cbind(xvec, Z)
        dmat = cbind(dmat, aux)
      }
    }
  } else {
    current_lag = 1
    Z = (weights_tensor * stages_tensor[[1]]) %*% NTS_data[t - 1, ]
    if (r_stages[current_lag] > 1) {
      for (r in 2:r_stages[current_lag]) {
        Z_r = (weights_tensor * stages_tensor[[r]]) %*% NTS_data[t - 1, ]
        Z = cbind(Z, Z_r)
      }
    }
    dmat = Z
  }
  return (dmat)
}


get_sample_dmat <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar = 'yes') {
  y = cbind(NTS_data[1 + p_lag, ])
  Z = get_time_dmat(stages_tensor, p_lag, 1 + p_lag, r_stages, NTS_data, weights_matrix, ar)
  for (t  in 2:(nrow(NTS_data) - p_lag)) {
    y = rbind(y, cbind(NTS_data[t + p_lag, ]))
    aux = get_time_dmat(stages_tensor, p_lag, t + p_lag, r_stages, NTS_data, weights_matrix, ar)
    Z = rbind(Z, aux)
  }
  dmat = as.matrix(cbind(y, Z))
  return (data.frame(dmat))
}


ls_gnar <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar='yes') {
  dummy_data = get_sample_dmat(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar)
  y = dummy_data[, 1]
  dmat = dummy_data[, 2:ncol(dummy_data)]
  out = lm(y ~ dmat + 0)
  # print(summary(out))
  return(out)
}

