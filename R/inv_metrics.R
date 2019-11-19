getHessian = function(fit, q) {
  Aqx = function(fit, q, r) {
    dx = 1e-5
    dr = dx * r
    (grad_log_prob(fit, q + dr / 2, adjust_transform = FALSE) -
        grad_log_prob(fit, q - dr / 2, adjust_transform = FALSE)) / dx
  }

  N = length(q)
  A = matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    r = rep(0, N)
    r[i] = 1
    A[, i] = Aqx(fit, q, r)
  }
  0.5 * (A + t(A))
}

diag_inv_metric = function(samples) {
  diag(diag(cov(samples)))
}

dense_inv_metric = function(samples, rank_check = TRUE) {
  c = cov(samples)

  rank = 0
  if(rank_check) {
    for(i in 1:(nrow(samples) - 1)) {
      r2 = sum(samples[i,] - samples[i + 1,])^2
      if(r2 > 1e-16 * ncol(samples)) {
        rank = rank + 1
      }
    }
  } else {
    rank = min(ncol(samples), nrow(samples) - 1)
  }

  nkeep = rank

  if(nkeep < ncol(samples)) {
    e = eigen(c, T)
    mine = e$values[nkeep]
    c = e$vectors[, 1:nkeep] %*% diag(e$values[1:nkeep] - mine) %*% t(e$vectors[, 1:nkeep])
    c = c + mine * diag(ncol(samples))
  }

  return(c)
}

lw_linear_corr_inv_metric = function(samples) {
  sqrt_D = diag(sqrt(diag(cov(samples))))
  sqrt_Dinv = diag(1 / diag(sqrt_D))
  sqrt_D %*% linshrink_cov(samples %*% sqrt_Dinv) %*% sqrt_D
}

compute_inv_metric = function(stan_fit, usamples) {
  Ntest = max(nrow(usamples) / 2, 50)
  Ntrain = nrow(usamples) - Ntest

  Ytrain = head(usamples, Ntrain)
  top_evec = eigen(cov(Ytrain), T)$vectors[, 1]
  Ytest = tail(usamples, Ntest)

  inv_metric_options = list(diag = diag_inv_metric,
                            dense = dense_inv_metric,
                            lw2004 = lw_linear_corr_inv_metric)

  perfdf = lapply(names(inv_metric_options), function(inv_metric_name) {
    inv_metric = inv_metric_options[[inv_metric_name]](Ytrain)

    cov_test = cov(Ytest)
    L = t(chol(inv_metric))
    el = eigen(solve(L, t(solve(L, t(cov_test)))), T)
    H = t(L) %*% getHessian(stan_fit, tail(usamples, 1)) %*% L
    eh = eigen(H, T)

    tibble(name = inv_metric_name,
           c_hybrid = sqrt(max(abs(eh$values)) * max(abs(el$values))))
  }) %>%
    bind_rows %>%
    arrange(-c_hybrid)

  print("Metric calculation info (lower c better, minimum is 1.0):")
  print(perfdf)

  name = perfdf %>% tail(1) %>% pull(name)

  inv_metric = inv_metric_options[[name]](usamples)

  return(inv_metric)
}
