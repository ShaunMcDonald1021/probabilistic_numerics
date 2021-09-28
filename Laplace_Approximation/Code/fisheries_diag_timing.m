% This script calculates computation time for the diagnostic applied to the
% fisheries data. The repetitions are to give us an accurate idea of
% variability in computation time. Make sure FSKQ is downloaded prior to
% running this (see README and Karvonen et al. 2018) and that
% GH_stuff_timing.m and fisheries_demo.R have been ran (in that order).
load('1970_demo_stuff.mat')
load('2005_demo_stuff.mat')

times_1970 = zeros([1 100]);
times_1970_GH = zeros([1 100]);
times_2005 = zeros([1 100]);
times_2005_GH = zeros([1 100]);

d = 72;
v = 25921;
gam = sqrt(1.5*(v+d)/(v+d-3));
lambda = 3.7;
lambda_GH = 1.73;

% No need to time this part, since it was already timed in
% GH_stuff_timing.m
XS = gh_seq(2);
us = sparse_gens(XS, d);
Us_GH = fss_gen(us);
s_star_GH = cell2mat(Us)';

for i = 1:100
    tic
    load('1970_diag.mat')
    Us = fss_gen(s_star(:,[1 73]));
    s_star = s_star';
    logf_interrs = logf_interrs';
    [~, w, wce, alph] = diag_calib(gam, lambda, d, v, Us, s_star);
    [post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, false, Us, wce, s_star);
    times_1970(i) = toc + diag_times1970(i);
    
    tic
    load('2005_diag.mat')
    Us = fss_gen(s_star(:,[1 73]));
    s_star = s_star';
    logf_interrs = logf_interrs';
    [~, w, wce, alph] = diag_calib(gam, lambda, d, v, Us, s_star);
    [post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, false, Us, wce, s_star);
    times_2005(i) = toc + diag_times2005(i);
    
    tic
    load('1970_diag_GH.mat')
    logf_interrs = logf_interrs';
    [~, w, wce, alph] = diag_calib(gam, lambda_GH, d, v, Us_GH, s_star_GH);
    [post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, false, Us_GH, wce, s_star_GH);
    times_1970_GH(i) = toc + diag_times1970_GH(i);
    
    tic
    load('2005_diag_GH.mat')
    logf_interrs = logf_interrs';
    [~, w, wce, alph] = diag_calib(gam, lambda_GH, d, v, Us_GH, s_star_GH);
    [post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, false, Us_GH, wce, s_star_GH);
    times_2005_GH(i) = toc + diag_times2005_GH(i);
end

save('diag_times.mat', 'times_1970', 'times_1970_GH', 'times_2005',...
    'times_2005_GH')
