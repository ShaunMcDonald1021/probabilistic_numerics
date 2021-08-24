source('probabilistic_numerics/Laplace_Approximation/Code/Laplace_diagnostic_POC/laplace_diag_functions.R')
#source('laplace_diag_functions.R')

package_get('mvtnorm')
package_get('TMB')
package_get('devtools')
package_get('nleqslv')
package_get('stockassessment',
            install_fun = function(c) devtools::install_github("fishfollower/SAM/stockassessment",
                                                               INSTALL_opts=c("--no-multiarch")))
package_get('rmatio')

# 1970 model
dat1970 = reduce(nscodData, year = c(1963:1969, 1976:2015), conf = nscodConf)
conf1970 = attr(dat1970, 'conf')
param1970 = defpar(dat1970, conf1970)
fit1970 = sam.fit(dat1970, conf1970, param1970, newtonsteps = 0, sim.condRE = FALSE,
                  inner.control = list(tol = 1e-10, tol10 = 1e-6, maxit = 5000))

diag_1970 = lap_diag_from_tmb(fit1970$obj)
d = diag_1970$d
#s_star = cbind(numeric(d), diag(sqrt(d), d), diag(-sqrt(d), d))
#s_star[,2:(d+1)] = s_star[,(d+1):2]
#logf_interrs_1970 = get_log_interrs(diag_1970$logf, diag_1970$T_mat, diag_1970$mode, s_star)

samp_sizes = 1e6*2^(-2:5)
imp_list = vector('list', length(samp_sizes))
tail_check_list = vector('list', length(samp_sizes))
for(i in seq_along(samp_sizes)){
  imp_list[[i]] = imp_sampler_t(samp_sizes[i], 25921, diag_1970$logf, diag_1970$T_mat, diag_1970$mode)
  tail_check_list[[i]] = tail_checker(imp_list[[i]]$outliers, diag_1970$logf, diag_1970$logf_at_mode,
                                      diag_1970$T_mat, diag_1970$mode)
}








