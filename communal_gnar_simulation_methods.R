gnar_sim_community <- function(nsims, network, weight_matrix, max_lags, r_stages, alpha_list, beta_list, covariate_groups, sigma=1) {
  d = ncol(weight_matrix)
  max_lag = max(max_lags)
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
    X = sim_all_communities(stages_tensor, weight_matrix, max_lags, alpha_list, beta_list, covariate_groups, lag_data)
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


sim_all_communities <- function(stages_tensor, weight_matrix, max_lags, alpha_list, beta_list, covariate_groups, lags_data) {
  d = ncol(weight_matrix)
  Xc = vapply(seq(1:length(covariate_groups)), function(x) { sim_time_step_community(stages_tensor, weight_matrix, max_lags[x],
                                                                                            alpha_list[[x]], beta_list[[x]], covariate_groups[[x]], 
                                                                                            lags_data)},
              rep(0, d))
  X = rowSums(Xc)
  return(X)
}

