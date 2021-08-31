#source('probabilistic_numerics/Laplace_Approximation/Code/Laplace_diagnostic_POC/laplace_diag_functions.R')
source('laplace_diag_functions.R')

#package_get('mvtnorm')
package_get('TMB')
package_get('devtools')
package_get('nleqslv')
package_get('parallel')
package_get('stockassessment',
            install_fun = function(c) devtools::install_github("fishfollower/SAM/stockassessment",
                                                               INSTALL_opts=c("--no-multiarch")))
package_get('rmatio')
package_get('mvnfast')

RNGkind("L'Ecuyer-CMRG")
set.seed(24601)
# 1970 model
dat1970 = reduce(nscodData, year = c(1963:1969, 1976:2015), conf = nscodConf)
conf1970 = attr(dat1970, 'conf')
param1970 = defpar(dat1970, conf1970)
fit1970 = sam.fit(dat1970, conf1970, param1970, newtonsteps = 0, sim.condRE = FALSE,
                  inner.control = list(tol = 1e-10, tol10 = 1e-6, maxit = 5000))$obj

diag_1970 = lap_diag_from_tmb(fit1970)
d = diag_1970$d
rm(dat1970, conf1970, param1970)
gc()
#s_star = cbind(numeric(d), diag(sqrt(d), d), diag(-sqrt(d), d))
#s_star[,2:(d+1)] = s_star[,(d+1):2]
#logf_interrs_1970 = get_log_interrs(diag_1970$logf, diag_1970$T_mat, diag_1970$mode, s_star)

samp_sizes = 1e6*2^(-2:6)
imp_list_1970 = vector('list', length(samp_sizes))
for(i in rev(seq_along(samp_sizes))){
  imp_list_1970[[i]] = imp_sampler_parallel(samp_sizes[i], 5, diag_1970$logf, diag_1970$logf_at_mode, diag_1970$T_mat, diag_1970$mode,
splitup = 1 + (i >= 9))
  gc()
  print(paste(samp_sizes[i], 'done'))
 .Random.seed = nextRNGStream(.Random.seed)
}

small_checkcon_times1970 = numeric(100)
small_checkcon_pvals1970 = numeric(100)
large_checkcon_times1970 = numeric(100)
large_checkcon_pvals1970 = numeric(100)
diag_times1970 = numeric(100)

# Repetition legitimizes
for(i in 1:100){
  small_start = proc.time()[3]
  small_checkcon = checkConsistency(fit_1970, par = diag_1970$theta_hat, n = 100)
  small_checkcon_times1970[i] = proc.time()[3] - small_start
  small_checkcon_pvals1970[i] = summary(small_checkcon)$marginal$p.value
  
  large_start = proc.time()[3]
  large_checkcon = checkConsistency(fit_1970, par = diag_1970$theta_hat, n = 1000)
  large_checkcon_times1970[i] = proc.time()[3] - large_start
  large_checkcon_pvals1970[i] = summary(large_checkcon)$marginal$p.value
  
  diag_start = proc.time()[3]
  s_star = cbind(numeric(d), diag(sqrt(d), d), diag(-sqrt(d), d))
  s_star[,2:(d+1)] = s_star[,(d+1):2]
  logf_interrs_1970 = get_log_interrs(diag_1970$logf, diag_1970$T_mat, diag_1970$mode, s_star)
  write.mat(append(list(logf_interrs = logf_interrs_1970, s_star = s_star), diag_1970[c(2,3,5,6)]),
            filename = 'rmatiotest.mat')
  diag_times1970[i] = proc.time()[3] - diag_start
}

save.image('1970_stuff')
