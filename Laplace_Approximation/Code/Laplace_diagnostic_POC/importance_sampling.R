library(mvtnorm)
library(nleqslv)
library(ismev)
library(SpatialExtremes)
library(TMB)
library(stockassessment)

# Fit SSM on last 4 years (2012-2015) of North Sea cod data
check_model_good = function(n){
  # Due to small amounts of data (I'm assuming), it's hard to get convergence and 
  poss_start_years = numeric(0)
  
  dat = reduce(nscodData, year = (1963+n):2015, conf = nscodConf)
  conf = attr(dat, "conf")
  param = defpar(dat, conf)
  tryCatch({
    fit = tryCatchsam.fit(dat, conf, param, silent = TRUE, newtonsteps = 0)
    if(all(!is.nan(unlist(fit$plsd))) & fit$opt$convergence == 0 & fit$sdrep$pdHess){
      poss_start_years = c(poss_start_years, 1963)
    }
  }, error = function(e){})
  
  for(i in 1963:(2015-n-1)){
    dat = reduce(nscodData, year = c(1963:i, (i+n+1):2015), conf = nscodConf)
    conf = attr(dat, "conf")
    tryCatch({
      param = defpar(dat, conf)
      fit = sam.fit(dat, conf, param, silent = TRUE, newtonsteps = 0)
      if(all(!is.nan(unlist(fit$plsd))) & fit$opt$convergence == 0 & fit$sdrep$pdHess){
        poss_start_years = c(poss_start_years, i+1)
      }
    }, error = function(e){})
  }
  dat = reduce(nscodData, year = 1963:(2015-n), conf = nscodConf)
  conf = attr(dat, "conf")
  param = defpar(dat, conf)
  tryCatch({
    fit = sam.fit(dat, conf, param, silent = TRUE, newtonsteps = 0)
    if(all(!is.nan(unlist(fit$plsd))) & fit$opt$convergence == 0 & fit$sdrep$pdHess){
      poss_start_years = c(poss_start_years, 2015 - n + 1)
    }
  }, error = function(e){})
  if(length(poss_start_years) == 0) stop(paste("No window of size", n, "converged properly. Try larger window."))
  else return(poss_start_years)
}

# Eventually we find that we can use the 6-year window from 2005-2010 (so 72-dimensional integral for LA)
# Can also try 6-year windows starting from 1970 or 1972
dat = reduce(nscodData, year = c(1963:1969, 1976:2015), conf = nscodConf)
# dat = nscodData
# conf = nscodConf
conf = attr(dat, "conf")
param = defpar(dat, conf)

fit = sam.fit(dat, conf, param, newtonsteps = 0, sim.condRE = FALSE,
              inner.control = list(tol = 1e-10, tol10 = 1e-6, maxit = 5000))

d = 12*dat$noYears

mode = fit$obj$env$last.par.best[length(fit$opt$par)+(1:d)]
neglog_f_mode = fit$obj$env$f(fit$obj$env$last.par.best)
logf = function(x) -fit$obj$env$f(c(fit$opt$par, x))

H = as.matrix(fit$obj$env$spHess(fit$obj$env$last.par.best))[length(fit$opt$par)+(1:d), length(fit$opt$par)+(1:d)]
U = chol(H)
lap_app = exp(-sum(log(diag(U))) - neglog_f_mode + d*log(2*pi)/2)

eig = eigen(-H)
G = eig$vectors %*% diag(1/sqrt(-eig$values))
# s_template = -3:3
# l = length(s_template)
# s = matrix(0, nrow = d, ncol = d*(l-1) + 1)
# s[1,1:l] = s_template
# for(i in 2:d) s[i,(i-2)*(l-1)+l+(1:(l-1))] = s_template[c(1:((l-1)/2), ((l+3)/2):l)]
s = as.matrix(read.csv('X.csv', header = FALSE))
s_trans = sweep(G %*% s, 1, mode, '+')
gauss_check = sapply(1:ncol(s_trans), FUN = function(i) logf(s_trans[,i]) + neglog_f_mode)

gauss_points = -colSums(s^2)/2
summary(exp(gauss_check) - exp(gauss_points))
chk = checkConsistency(fit$obj)
summary(chk)$marginal

imp_sampler = function(n, sigma){
  X = rmvnorm(n, numeric(d), diag(nrow = d))
  samps = sweep(backsolve(U/sigma, t(X)), 1, mode, "+")
  return(exp(apply(samps, 2, logf) - sum(log(diag(U/sigma)))  - dmvnorm(X, log = TRUE)))
}

imp_sampler_t = function(n, sigma, nu){
  X = rmvt(n, sigma = sigma^2*solve(H), df = nu, delta = mode, type = 'shifted')
 # samps = sweep(backsolve(U/sigma, t(X)), 1, mode, "+")
  return(exp(apply(X, 1, logf)  - dmvt(X, sigma = solve(H), df = nu, 
                                                       delta = mode, type = 'shifted', log = TRUE)))
}

# imp_sampler = function(n){
#   X = rmvnorm(n, numeric(d), diag(nrow = d))
#   samps = sweep(backsolve(U, t(X)), 1, mode, "+")
#   return(exp(apply(samps, 2, function(x) logf(x)) + neglog_f_mode - dmvnorm(X, log = TRUE)))
# }

# Wald test for infinite variance of importance sampler (Koopman et al. 2009)
xi_tau = function(n, tau, Z){
  return(mean(log(1+tau*Z)))
  print(tau)
  print(1+tau*Z)
}

tau_fun = function(tau, n, Z){
  return((1 + 1/xi_tau(n, tau, Z))*sum(Z/(1+tau*Z)) - n/tau)
}

wald_test = function(imp_samps, thresh){
  Z = imp_samps[imp_samps > thresh] - thresh
  n = length(Z)

  tau_hat = nleqslv(x = 1/mean(Z), fn = function(tau) tau_fun(tau, n, Z))$x

  xi_hat = xi_tau(n, tau_hat, Z)

  beta_hat = xi_hat/tau_hat

  t_stat = sqrt(n/(3*beta_hat^2))*(xi_hat - 0.5)

  return(list(xi_hat = xi_hat, beta_hat = beta_hat, t_stat = t_stat))
}

wald_test2 = function(imp_samps, thresh){
  mle = gpd.fit(imp_samps, thresh)$mle
  
  xi_hat = mle[2]
  
  beta_hat = mle[1]
  
  t_stat = sqrt(length(which(imp_samps > thresh))/(3*beta_hat^2))*(xi_hat - 0.5)
  
  return(list(xi_hat = xi_hat, beta_hat = beta_hat, t_stat = t_stat))
}

imp_samp_test = function(n, eps){
  X = rmvnorm(n, numeric(d), diag(1/(1+eps), nrow = d))
  
  samps = exp(dmvnorm(X, log = TRUE) - dmvnorm(X, sigma = diag(1/(1+eps), nrow = d), log = TRUE))
  
  wald_test2(samps, quantile(samps, 0.55))
}

eval_grid = seq(-3, 3, by = 0.01)

gauss_dir_checker = function(weights, sigma = 1){
  weights = weights/sqrt(sum(weights^2))
  direction = as.numeric(G %*% weights)
  # direction = direction/sqrt(sum(direction^2))
  fpoints = numeric(length(eval_grid))
  quadpoints = numeric(length(eval_grid))
  for(i in 1:length(eval_grid)){
    fpoints[i] = logf(mode + eval_grid[i]*direction) + neglog_f_mode
    quadpoints[i] = -(eval_grid[i]*direction)%*%(H/sigma^2)%*%(eval_grid[i]*direction)/2
  }
  # plot(eval_grid, fpoints, type = 'l')
  # lines(eval_grid, quadpoints, col = 'red')
  
  plot(eval_grid, exp(fpoints), type = 'l')
  lines(eval_grid, exp(quadpoints), col = 'red')
}

