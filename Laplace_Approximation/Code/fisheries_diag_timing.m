% This script calculates computation time for the diagnostic applied to the
% fisheries data. The repetitions are to give us an accurate idea of
% variability in computation time. Make sure FSKQ is downloaded prior to
% running this (see README and Karvonen et al. 2018) and that
% GH_stuff_timing.m and fisheries_demo.R have been ran (in that order).

% From fisheries_demo.R
load('1970_demo_stuff.mat')
load('2005_demo_stuff.mat')
% From GH_stuff_timing.m
load('GH_times.mat')

times_1970 = zeros([1 100]);
times_1970_GH = zeros([1 100]);
times_2005 = zeros([1 100]);
times_2005_GH = zeros([1 100]);

% Diagnostic ingredients
d = 72;
v = 25921;
gam = sqrt(1.5*(v+d)/(v+d-3));
% Lambdas determined empirically
lambda = 3.7;
lambda_GH = 2.06;

% No need to time this part, since it was already timed in
% GH_stuff_timing.m
XS = gh_seq(2);
us = 3.6*sparse_gens(XS, d);
Us_GH = fss_gen(us);
Us_GH = Us_GH(2:end);
s_star_GH = cell2mat(Us_GH)';

% 100 repetitions
for i = 1:100
    % 1970 data w/simple grid
    tic
    load('1970_diag.mat')
    Us = fss_gen(s_star(:,[1 73]));
    s_star = s_star';
    logf_interrs = logf_interrs';
    [~, w, wce, alph] = diag_calib(gam, lambda, d, v, Us, s_star);
    [post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, false, Us, wce, s_star);
    times_1970(i) = toc + diag_times1970(i);
    
    % 2005 data w/simple grid
    tic
    load('2005_diag.mat')
    Us = fss_gen(s_star(:,[1 73]));
    s_star = s_star';
    logf_interrs = logf_interrs';
    [~, w, wce, alph] = diag_calib(gam, lambda, d, v, Us, s_star);
    [post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, false, Us, wce, s_star);
    times_2005(i) = toc + diag_times2005(i);
    
    % 1970 data w/higher-order Gauss-Hermite grid
    tic
    load('1970_diag_GH.mat')
    logf_interrs = logf_interrs';
    [~, w, wce, alph] = diag_calib(gam, lambda_GH, d, v, Us_GH, s_star_GH);
    [post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, false, Us_GH, wce, s_star_GH);
    times_1970_GH(i) = toc + diag_times1970_GH(i) + GH_times(i);
    
    % 2005 data w/higher-order Gauss-Hermite grid
    tic
    load('2005_diag_GH.mat')
    logf_interrs = logf_interrs';
    [~, w, wce, alph] = diag_calib(gam, lambda_GH, d, v, Us_GH, s_star_GH);
    [post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, false, Us_GH, wce, s_star_GH);
    times_2005_GH(i) = toc + diag_times2005_GH(i) + GH_times(i);
end

save('diag_times.mat', 'times_1970', 'times_1970_GH', 'times_2005',...
    'times_2005_GH')
