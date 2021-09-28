# Code for stock assessment experiments as shown in Laplace diagnostic manuscript
# NOTE: this demo is designed for Unix systems. The importance sampler uses forking
# and WILL NOT WORK IN WINDOWS unless you manually change the number of cores to 1

# Assuming you run this demo from the same location where laplace_diag_functions.R is.
# If not, change the below line accordingly
source('laplace_diag_functions.R')

# Get necessary packages
package_get('devtools')
package_get('stockassessment',
            install_fun = function(c) devtools::install_github("fishfollower/SAM/stockassessment",
                                                               INSTALL_opts=c("--no-multiarch")))
# Set seed for consistency in parallelized PRNG
RNGkind("L'Ecuyer-CMRG")
set.seed(24601)

#### 1970 model ####
dat1970 = reduce(nscodData, year = c(1963:1969, 1976:2015), conf = nscodConf)
conf1970 = attr(dat1970, 'conf')
param1970 = defpar(dat1970, conf1970)
fit_1970 = sam.fit(dat1970, conf1970, param1970, newtonsteps = 0, sim.condRE = FALSE, silent = TRUE,
                   inner.control = list(tol = 1e-10, tol10 = 1e-6, maxit = 5000))$obj
# Ensure last.par == last.par.best
fit_1970$fn(fit_1970$env$last.par.best[fit_1970$env$lfixed()])

diag_1970 = lap_diag_from_tmb(fit_1970)
d = diag_1970$d

# Clear up some memory
rm(dat1970, conf1970, param1970)
gc()

# Importance samplers of increasing sizes
# IMPORTANT: manually set cores = 1 in imp_sampler_parallel if running on Windows
samp_sizes = 1e6*2^(-2:7)
splitup_factor = c(rep(1, 7), 2, 4, 8)
threshes = c(rep(0.9, 7), rep(0.95, 2), 0.975)
imp_list_1970 = vector('list', length(samp_sizes))
for(i in rev(seq_along(samp_sizes))){
  imp_list_1970[[i]] = imp_sampler_parallel(samp_sizes[i], 5, diag_1970$logf,
                                            diag_1970$logf_at_mode, diag_1970$T_mat,
                                            diag_1970$mode, splitup = splitup_factor[i],
                                            thresh = threshes[i])
  gc() # Can never have too much memory
  print(paste('Importance sampler of size', samp_sizes[i], 'done for 1970 model'))
}

# Repeatedly run and time various Laplace diagnostics:
# checkConsistency with a small number of samples (100) and a large number (1000),
# as well as our BQ diagnostic with two different interrogation grids
small_checkcon_times1970 = numeric(100)
small_checkcon_pvals1970 = numeric(100)
large_checkcon_times1970 = numeric(100)
large_checkcon_pvals1970 = numeric(100)
diag_times1970 = numeric(100) # Time for diagnostic w/CKF grid
diag_times1970_GH = numeric(100) # Time for diagnostic w/Gauss-Hermite grid

# Repetition legitimizes
for(i in 1:100){
  small_start = proc.time()[3]
  small_checkcon = checkConsistency(fit_1970, par = diag_1970$theta, n = 100)
  small_checkcon_pvals1970[i] = summary(small_checkcon)$marginal$p.value
  small_checkcon_times1970[i] = proc.time()[3] - small_start
  
  large_start = proc.time()[3]
  large_checkcon = checkConsistency(fit_1970, par = diag_1970$theta, n = 1000)
  large_checkcon_pvals1970[i] = summary(large_checkcon)$marginal$p.value
  large_checkcon_times1970[i] = proc.time()[3] - large_start
  
  # Recalculate diagnostic components for more accurate estimates of
  # variability in computation time
  component_time = lap_diag_from_tmb(fit_1970)$time
  
  diag_start = proc.time()[3]
  s_star = cbind(numeric(d), diag(sqrt(d), d), diag(-sqrt(d), d))
  s_star[,2:(d+1)] = s_star[,(d+1):2]
  logf_interrs_1970 = get_log_interrs(diag_1970$logf, diag_1970$T_mat, diag_1970$mode, s_star,
                                      TRUE, diag_1970$logf_at_mode, diag_1970$log_T_det, '1970_diag.mat')
  # This is somewhat inefficient, but saving the interrogation values is certainly part of
  # the computational cost so it should be timed as well.
  diag_times1970[i] = proc.time()[3] - diag_start + component_time
  
  # Repeat for the GH grid
  diag_start = proc.time()[3]
  GH_stuff = read.mat('GH_stuff.mat')
  logf_interrs_1970 = get_log_interrs(diag_1970$logf, diag_1970$T_mat, diag_1970$mode, GH_stuff$s_star,
                                      TRUE, diag_1970$logf_at_mode, diag_1970$log_T_det, '1970_diag_GH.mat')
  diag_times_1970_GH[i] = proc.time()[3] - diag_start + component_time + GH_stuff$time[i]
  
  print(paste('Diagnostic replication', i, 'done for 1970 model'))
}

save.image('1970_stuff.RData')
write.mat(list(imp_means = unlist(lapply(imp_list_1970, '[[', 'mean')),
               imp_vars = unlist(lapply(imp_list_1970, '[[', 'var')), diag_times = diag_times1970,
               diag_times_GH = diag_times_1970_GH,small_pvals = small_checkcon_pvals1970,
               large_pvals = large_checkcon_pvals1970))

# Get rid of all the 1970 stuff to save space
rm(list = grep('1970', ls(), value = TRUE))
gc()

#### 2005 model ####
# Same things as for the 1970 model
dat2005 = reduce(nscodData, year = c(1963:2004, 2011:2015), conf = nscodConf)
conf2005 = attr(dat2005, 'conf')
param2005 = defpar(dat2005, conf2005)
# This one needs lower tolerance for good numerical behaviour
fit_2005 = sam.fit(dat2005, conf2005, param2005, newtonsteps = 0, sim.condRE = FALSE, silent = TRUE,
                   inner.control = list(tol = 1e-10, tol10 = 1e-6, maxit = 5000))$obj
# Ensure last.par == last.par.best
fit_2005$fn(fit_2005$env$last.par.best[fit_2005$env$lfixed()])

diag_2005 = lap_diag_from_tmb(fit_2005)

# Clear up some memory
rm(dat2005, conf2005, param2005)
gc()

# Importance samplers of increasing sizes
# IMPORTANT: manually set cores = 1 in imp_sampler_parallel if running on Windows
imp_list_2005 = vector('list', length(samp_sizes))
for(i in rev(seq_along(samp_sizes))){
  imp_list_2005[[i]] = imp_sampler_parallel(samp_sizes[i], 5, diag_2005$logf,
                                            diag_2005$logf_at_mode, diag_2005$T_mat,
                                            diag_2005$mode, splitup = splitup_factor[i],
                                            thresh = threshes[i])
  gc()
  print(paste('Importance sampler of size', samp_sizes[i], 'done for 2005 model'))
}

# Repeatedly run and time various Laplace diagnostics:
small_checkcon_times2005 = numeric(100)
small_checkcon_pvals2005 = numeric(100)
large_checkcon_times2005 = numeric(100)
large_checkcon_pvals2005 = numeric(100)
diag_times2005 = numeric(100)
diag_times2005_GH = numeric(100)

# Repetition legitimizes
for(i in 1:100){
  small_start = proc.time()[3]
  small_checkcon = checkConsistency(fit_2005, par = diag_2005$theta, n = 100)
  small_checkcon_pvals2005[i] = summary(small_checkcon)$marginal$p.value
  small_checkcon_times2005[i] = proc.time()[3] - small_start
  
  large_start = proc.time()[3]
  large_checkcon = checkConsistency(fit_2005, par = diag_2005$theta, n = 1000)
  large_checkcon_pvals2005[i] = summary(large_checkcon)$marginal$p.value
  large_checkcon_times2005[i] = proc.time()[3] - large_start
  
  component_time = lap_diag_from_tmb(fit_2005)$time
  diag_start = proc.time()[3]
  s_star = cbind(numeric(d), diag(sqrt(d), d), diag(-sqrt(d), d))
  s_star[,2:(d+1)] = s_star[,(d+1):2]
  logf_interrs_2005 = get_log_interrs(diag_2005$logf, diag_2005$T_mat, diag_2005$mode, s_star,
                                      TRUE, diag_2005$logf_at_mode, diag_2005$log_T_det, '2005_diag.mat')
  # This is somewhat inefficient, but saving the interrogation values is certainly part of
  # the computational cost so it should be timed as well.
  diag_times2005[i] = proc.time()[3] - diag_start + component_time
  
  # Repeat for the GH grid
  diag_start = proc.time()[3]
  GH_stuff = read.mat('GH_stuff.mat')
  logf_interrs_2005 = get_log_interrs(diag_2005$logf, diag_2005$T_mat, diag_2005$mode, GH_stuff$s_star,
                                      TRUE, diag_2005$logf_at_mode, diag_2005$log_T_det, '2005_diag_GH.mat')
  diag_times_2005_GH[i] = proc.time()[3] - diag_start + component_time + GH_stuff$time[i]
  
  print(paste('Diagnostic replication', i, 'done for 2005 model'))
}

save.image('2005_stuff.RData')
write.mat(list(imp_means = unlist(lapply(imp_list_2005, '[[', 'mean')),
               imp_vars = unlist(lapply(imp_list_2005, '[[', 'var')), diag_times = diag_times2005,
               diag_times_GH = diag_times_2005_GH,small_pvals = small_checkcon_pvals2005,
               large_pvals = large_checkcon_pvals2005))
