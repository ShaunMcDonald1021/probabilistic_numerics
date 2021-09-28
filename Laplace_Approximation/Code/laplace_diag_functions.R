# Dependencies: TMB, nleqslv, parallel, rmatio, mvnfast

package_get = function(package_name, install_fun = function(c) install.packages(package_name)){
  # Function to automatically attach packages from library
  # If they're not already installed, it'll install them, THEN attach
  
  # package_name: name of package as character object
  # install_fun: function to install package if not already installed
  
  tryCatch(library(package_name, character.only = TRUE), error = install_fun,
           finally = library(package_name, character.only = TRUE))
}

# Get all dependencies
package_get('nleqslv')
package_get('parallel')
package_get('TMB')
package_get('rmatio')
package_get('mvnfast')

lap_diag_from_tmb = function(obj){
  # Given a TMB object (obj), extracts all the necessary components for Laplace diagnostic
  
  # Returns:
  # d: dimensionality of x
  # mode: numeric vector of size d. Mode w.r.t. x (i.e. $\hat{x}(\theta)$)
  # theta: Value of fixed effects/parameters
  # logf_at_mode: $\log f(\hat{x})$
  # logf: $\log f$ as a function of x with $\theta$ fixed
  # T_mat: transformation matrix based on eigendecomposition of Hessian
  # log_T_det: log of determinant of T
  # eig_time: time (in seconds) taken for eigendecomposition of Hessian
  
  d = length(obj$env$random)
  
  last.par = obj$env$last.par
  # Using last.par instead of last.par.best makes it easier to fold into model fitting
  # If you want to use the diagnostic on the fitted parameter ($\hat{\theta}$),
  # make sure obj$env$last.par == obj$env$last.par.best before calling this function
  
  theta = as.numeric(last.par[obj$env$lfixed()])
  # theta also required for checkConsistency so we exclude that from timing
  start = proc.time()[3] 
  mode = as.numeric(last.par[obj$env$random]) 
  
  logf_at_mode = -obj$env$f(last.par) 
  
  logf = function(x) -obj$env$f(c(theta, x)) 
  
  neg_H = as.matrix(obj$env$spHess(last.par))[obj$env$random, obj$env$random]
  # Negative of Hessian (w.r.t. x) evaluated at $\hat{x}(\theta)$
  
  # Eigendecomposition of -H
  # -H has the same eigenvectors as $-H^{-1}$, and its eigenvalues are simply the inverses
  # of those of -H^{-1}
  eig = eigen(neg_H, symmetric = TRUE)
  T_mat = eig$vectors %*% diag(1/sqrt(eig$values))
  log_T_det = -sum(log(eig$values))/2
  
  time = proc.time()[3] - start
  return(list(d = d, mode = mode, theta = theta, logf_at_mode = logf_at_mode,
              logf = logf, T_mat = T_mat, log_T_det = log_T_det, time = time))
}


get_log_interrs = function(logf, T_mat, mode, s_star, write_mat = FALSE,
                           logf_at_mode, log_T_det, filename){
  # Returns interrogation values on the log scale and optionally writes
  # everything required for the diagnostic to a .mat file.
  
  # logf: a function of x only
  # T_mat: the matrix T defined in the manuscript
  # mode: vector equal to dimension of x giving the mode $\hat{x}$
  # s_star: the grid of preliminary interrogation points
  # write_mat: If TRUE, writes everything required for the diagnostic to a .mat file
  # logf_at_mode & log_T_det, are as in lap_diag_from_tmb and are only required
  # filename: a character string for the name of the .mat file to write. Only required
  # when write_mat = TRUE
  # if write_mat = TRUE
  s = sweep(T_mat %*% s_star, 1, mode, '+')
  logf_interrs = apply(s, 2, logf)
  
  if(write_mat){
    write.mat(list(logf_interrs = logf_interrs, mode = mode, s_star = s_star,
                   logf_at_mode = logf_at_mode, log_T_det = log_T_det), filename = filename)
  }
  
  return(logf_interrs)
}


tail_checker = function(outliers, logf, logf_at_mode, T_mat, mode, sigma = 1,
                        scale_mat = sigma^2*T_mat%*%t(T_mat), deccheck = 20, zerocheck= 10){
  # Suppose we estimate $\int f(x) \mathrm{d}x$ using importance sampling (say, using the
  # multivariate T as the proposal distribution). Given a vector $x \in \mathbb{R}^d$ for which
  # the importance weight w(x) is large, this function checks if the tail of $f^2$ decays
  # faster than a Gaussian density in the direction of x (from the mode).
  # This provides *some* empircal assurance that the importance sampler has finite variance.
  
  # outliers: d-by-n matrix of n vectors in $\mathbb{R}^d$, corresponding to the samples
  # which gave the largest importance weights
  # logf, logf_at_mode, T_mat, mode are as in lap_diag_from_tmb
  # sigma: optional scaling factor for Gaussian density
  # scale_mat should ALWAYS be sigma^2 * T_mat %*% t(T_mat). Can be precomputed for efficiency
  # deccheck: number of tail values used to check that tails are monotonically decreasing
  # zerocheck: number of tail values used to check that f^2/Gaussian eventually goes to 0
  
  num_outliers = ncol(as.matrix(outliers))
  tails = vector('list', num_outliers)
  is_decreasing = logical(num_outliers)
  is_zero = logical(num_outliers)
  #scale_mat = sigma^2 * T_mat %*% t(T_mat)
  
  for(i in 1:num_outliers){
    # First, "standardize" the outlier
    out_cent = solve(T_mat, outliers[,i] - mode)
    out_norm = sqrt(sum(out_cent^2))
    out_dir = out_cent/out_norm # Unit vector in the direction of out_cent
    # Scale and rotate out_dir by T_mat (see Sec. 4.1 of manuscript)
    out_vec = as.numeric(T_mat %*% out_dir)
    
    # Check starts at the outlier itself and goes WAAAAAAY out into the tail
    eval_grid = seq(floor(out_norm), 1500, by = 5)
    # Divide both $f^2$ and normal density by their values at mode for numerical convenience
    diff_logs = rep(dmvn(mode, sigma = scale_mat, mu = mode, log = TRUE), length(eval_grid))
    for(j in seq_along(eval_grid)){
      diff_logs[j] = diff_logs[j] + 2*(logf(eval_grid[j]*out_vec + mode) - logf_at_mode) -
        dmvn(eval_grid[j]*out_vec + mode, sigma = scale_mat, mu = mode, log = TRUE)
    }
    tails[[i]] = diff_logs[!is.nan(diff_logs)] # logf may give NaN's for values way out in tails
    
    is_decreasing[i] = all(diff(tail(tails[[i]], deccheck)) < 0)
    is_zero[i] = all(exp(tail(tails[[i]], zerocheck)) == 0) # log values should be LARGE negatives
  }
  
  if(all(is_decreasing) & all(is_zero)){
    print('Tails in directions of outliers look okay')
  } else{
    if(!all(is_decreasing)){
      print('Some tails do not monotonically decrease')
      print('Outlier sample(s):')
      print(apply(outliers[,!is_decreasing], 2, function(x) paste(x, collapse = ", ")))
      print('Diffs in tail(s):')
      lapply(tails[!is_decreasing], function(x) print(diff(tail(x, deccheck))))
    }
    if(!all(is_zero)) {
      print('Some tails do not decay to zero')
      print('Outlier sample(s):')
      print(apply(outliers[,!is_zero], 2, function(x) paste(x, collapse = ", ")))
      print('Tail(s):')
      lapply(tails[!is_zero], function(x) print(exp(tail(x, zerocheck))))
    }
  }
  
  #return(list(tails = tails, is_decreasing = is_decreasing, is_zero = is_zero))
}


imp_sampler_t = function(N, nu, logf, logf_at_mode, T_mat, mode, sigma = 1,
                         scale_mat = sigma^2*T_mat%*%t(T_mat), num_outliers = 10){
  # Generates importance weights to estimate $\int f(x) \mathrm{d}x$ using a multivariate
  # T distribution as the proposal. Also checks tail behaviour in directions of samples
  # for which weights are large
  
  # N: number of samples
  # nu: degrees of freedom for T distribution
  # logf, logf_at_mode, T_mat, mode are as in lap_diag_from_tmb
  # scale_mat should ALWAYS be sigma^2 * T_mat %*% t(T_mat). Can be precomputed for efficiency
  # num_outliers: number of outliers to check (i.e. check tail behaviour in directions of
  # samples giving the [num_outliers] largest weights)
  # sigma: optional scaling factor for proposal distribution
  
  X = rmvt(N, sigma = scale_mat, df = nu, mu = mode)
  weights = exp(apply(X, 1, logf) - dmvt(X, sigma = scale_mat, df = nu, mu = mode, log = TRUE))
  
  tops = which(weights > quantile(weights, 1-num_outliers/N, na.rm = TRUE))
  outliers = t(X[tops,])
  # Get rid of X to free up memory
  rm(X)
  gc()
  
  # Check tail behaviour in directions of outliers
  tail_checker(outliers, logf, logf_at_mode, T_mat, mode, sigma, scale_mat)
  
  return(weights)
}


imp_sampler_parallel = function(N, nu, logf, logf_at_mode, T_mat, mode, splitup = 1,
                                num_outliers = 10, sigma = 1,
                                thresh = 0.9, cores = 6){
  # Parallelized importance sampling with a multivariate T proposal distribution.
  # This relies on forking, so it only works on Windows if cores = 1.
  
  # N, nu, logf, logf_at_mode, T_mat, mode, num_outliers, sigma are as in imp_sampler_t
  # splitup: How many serial "rounds" of importance sampling? e.g. if N = 2e+6
  # and splitup = 2, it runs a parallel sampler of size 1e+6, then does it again and
  # appends the results. Useful for large N when there might not be enough memory to do
  # it all at once.
  # thresh: Only keep the top 100*(1-thresh)% of weights for later diagnostics. Saves space.
  # cores: number of cores to use for parallel sampling. MUST BE 1 ON WINDOWS MACHINES
  
  # Returns:
  # topweights: the largest 100*(1-thresh)% of importance weights
  # time: How long it took in seconds
  # thresh_quant: (thresh)-th quantile of importance weights
  # mean: the actual estimate of the integral
  # var: the variance of the weights. var/N is a sensible estimate of the integral variance
  
  start = proc.time()[3]
  scale_mat = sigma^2 * T_mat %*% t(T_mat)
  
  weights = numeric(0) 
  for(i in 1:splitup){
    weights = c(weights,
                pvec(1:(N/splitup),
                     function(x) imp_sampler_t(length(x), nu, logf, logf_at_mode, T_mat,
                                               mode, sigma, scale_mat, num_outliers),
                     mc.cores = cores))
    assign('.Random.seed', nextRNGSubStream(.Random.seed), pos = .GlobalEnv) 
    # Otherwise we get repeated weights from multiple pvec calls
  }
  time = proc.time()[3] - start # Pretty rough estimate of timing b/c of tail checking
  
  thresh_quant = quantile(weights, thresh, na.rm = TRUE)
  
  return(list(topweights = weights[weights >= thresh_quant], time = time,
              thresh_quant = thresh_quant, mean = mean(weights, na.rm = TRUE),
              var = var(weights, na.rm = TRUE)))
}


loglik = function(xi, beta, Z){
  # log-likelihood of a sample Z from a generalized Pareto distribution with
  # shape parameter xi and scale parameter beta
  -length(Z)*log(beta) - (1+1/xi)*sum(log(1+xi*Z/beta))
}


xi_score = function(xi, beta, Z){
  # Derivative of loglik w.r.t. xi
  sum(log(1+xi*Z/beta))/xi^2 - (1+1/xi)*sum(Z/(beta+xi*Z))
}


beta_score = function(xi, beta, Z){
  # Derivative of loglik w.r.t. beta
  ((xi+1)*sum(Z/(beta+xi*Z)) - length(Z))/beta
}


imp_score_test = function(weights, thresh){
  # The score test from Jan Koopman et al. (2009) to test if an IS has finite variance
  # Assumes the excesses of the weights follow a GPD, then tests $\xi \leq 0.5$
  # weights: vector of (largest) weights from an importance sampler
  # thresh: The threshold to use in modelling the tail of the weight distribution
  # Must have thresh >= min(weights)
  
  # Returns:
  # score_stat: the test statistic
  # beta_hat: the estimated scale parameter to check that it's a true optimum
  # (i.e. can manually check that beta_score(0.5, beta_hat, Z) is close to 0)
  Z = weights[weights >= thresh] - thresh
  
  beta_hat = optim(par = 0.5*mean(Z), fn = function(beta) loglik(0.5, beta, Z),
                   gr = function(beta) beta_score(0.5, beta, Z), method = 'L-BFGS-B',
                   control = list(factr = 1e-8, maxit = 1000, fnscale = -1,
                                  parscale = 0.0001, pgtol = 0), lower = 1e-18)$par
  score_stat = xi_score(0.5, beta_hat, Z)/sqrt(2*length(Z))
  return(list(score_stat = score_stat, beta_hat = beta_hat))
}
