package_get = function(package_name, install_fun = function(c) install.packages(package_name)){
  tryCatch(library(package_name, character.only = TRUE), error = install_fun, finally = library(package_name, character.only = TRUE))
}

lap_diag_from_tmb = function(obj){
  
  d = length(obj$env$random)
  
  last.par.best = obj$env$last.par.best
  
  mode = last.par.best[obj$env$random]
  
  theta_hat = as.numeric(last.par.best[obj$env$lfixed()])
  
  logf_at_mode = -obj$env$f(last.par.best)
  
  logf = function(x) -obj$env$f(c(theta_hat, x))
  
  neg_H = as.matrix(obj$env$spHess(last.par.best))[obj$env$random, obj$env$random]
  
  start = proc.time()
  eig = eigen(neg_H, symmetric = TRUE)
  T_mat = eig$vectors %*% diag(1/sqrt(eig$values))
  log_T_det = -sum(log(eig$values))/2
  eig_time = proc.time() - start
  
  return(list(d = d, mode = mode, logf_at_mode = logf_at_mode, logf = logf, T_mat = T_mat, log_T_det = log_T_det, eig_time = eig_time))
}

get_log_interrs = function(logf, T_mat, mode, s_star){
  s = sweep(T_mat %*% s_star, 1, mode, '+')
  apply(s, 2, logf)
}

imp_sampler_t = function(N, nu, logf, T_mat, mode, num_outliers = 10, thresh = 0.9, sigma = 1){
  scale_mat = sigma^2 * T_mat %*% t(T_mat)
  
  start = proc.time()
  X = rmvt(N, sigma = scale_mat, df = nu, delta = mode, type = 'shifted')
  # mode_diff = logf_at_mode - dmvt(mode, sigma = scale_mat, df = nu, delta = mode, type = 'shifted',
  #                                         log = TRUE, checkSymmetry = FALSE)
  weights = exp(apply(X, 1, logf) - dmvt(X, sigma = scale_mat, df = nu, delta = mode, type = 'shifted',
                                       log = TRUE, checkSymmetry = FALSE))
  time = proc.time() - start
  
  # Get indices of num_outliers largest weights to check tail behaviour
  tops = which(weights > quantile(weights, 1-num_outliers/N))
  
  thresh_quant = quantile(weights, thresh)
  
  # Only return top 100*(1-thresh)% of weights to save space
  return(list(topweights = weights[weights >= thresh_quant], time = time, thresh_quant = thresh_quant,
              outliers = t(X[tops,]), nans = t(X[is.nan(weights),]), mean = mean(weights), var = var(weights)))
}

xi_score = function(xi, beta, Z){
  sum(log(1+xi*Z/beta))/xi^2 - (1+1/xi)*sum(Z/(beta+xi*Z))
}

beta_score = function(xi, beta, Z){
  ((xi+1)*sum(Z/(beta+xi*Z)) - length(Z))/beta
}

loglik = function(xi, beta, Z){
  -length(Z)*log(beta) - (1+1/xi)*sum(log(1+xi*Z/beta))
}

imp_score_test = function(weights, thresh){
  Z = weights[weights >= thresh] - thresh
  # beta_hat = nleqslv(x = 0.5*mean(Z), fn = function(beta) beta_score(0.5, beta, Z),
  #                    control = list(ftol = 1e-18, xtol = 1e-18, btol = 1e-18))$x
  
  beta_hat = optim(par = 0.5*mean(Z), fn = function(beta) loglik(0.5, beta, Z),
                   gr = function(beta) beta_score(0.5, beta, Z), method = 'L-BFGS-B',
                   control = list(factr = 1e-8, maxit = 1000, fnscale = -1, parscale = 0.0001,
                                  pgtol = 0), lower = 1e-18)$par
  score_stat = xi_score(0.5, beta_hat, Z)/sqrt(2*length(Z))
  return(list(score_stat = score_stat, beta_hat = beta_hat))
}

tail_checker = function(outliers, logf, logf_at_mode, T_mat, mode, sigma = 1, deccheck = 20, zerocheck= 10){
  num_outliers = ncol(as.matrix(outliers))
  tails = vector('list', num_outliers)
  is_decreasing = logical(num_outliers)
  is_zero = logical(num_outliers)
  
  for(i in 1:num_outliers){
    out_cent = solve(T_mat, outliers[,i] - mode)
    out_norm = sqrt(sum(out_cent^2))
    out_dir = out_cent/out_norm
    out_vec = T_mat %*% out_dir
    scale_mat = sigma^2 * T_mat %*% t(T_mat)
    
    eval_grid = seq(floor(out_norm), 1000, by = 5)
    diff_logs = rep(dmvnorm(as.numeric(mode), sigma = scale_mat, mean = as.numeric(mode), log = TRUE), length(eval_grid))
    for(j in seq_along(eval_grid)){
      diff_logs[j] = diff_logs[j] + 2*(logf(eval_grid[j]*out_vec + mode) - logf_at_mode) -
        dmvnorm(as.numeric(eval_grid[j]*out_vec + mode), sigma = scale_mat, mean = as.numeric(mode), log = TRUE)
    }
    tails[[i]] = diff_logs[!is.nan(diff_logs)]
    
    is_decreasing[i] = all(diff(tail(tails[[i]], deccheck)) < 0)
    is_zero[i] = all(exp(tail(tails[[i]], zerocheck)) == 0)
  }

  if(all(is_decreasing) & all(is_zero)) print('Tails in directions of outliers look okay')
  else{
    if(!all(is_decreasing)) print('One of the tails does not monotonically decrease')
    if(!all(is_zero)) print('One of the tails does not decay to zero')
  }
  
  return(list(tails = tails, is_decreasing = is_decreasing, is_zero = is_zero))
}
  
  
