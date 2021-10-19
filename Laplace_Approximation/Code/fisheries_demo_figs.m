
load('1970_demo_stuff.mat')
load('2005_demo_stuff.mat')
load('diag_times.mat')

d = 72;
v = 25921;
gam = sqrt(1.5*(v+d)/(v+d-3));
lambda = 3.7;
lambda_GH = 2.06;

XS = gh_seq(2);
us = 3.6*sparse_gens(XS, d);
%us = zeros([72 3]);
%us(1:2,2) = 6;
%us(1, 3) = sqrt(d);
Us_GH = fss_gen(us);
Us_GH = Us_GH(2:end);
s_star_GH = cell2mat(Us_GH)';

load('1970_diag.mat')
Us = fss_gen(s_star(:,[1 73]));
s_star = s_star';
logf_interrs = logf_interrs';
[~, w, wce, alph] = diag_calib(gam, lambda, d, v, Us, s_star);
[post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
    log_T_det, d, gam, alph, w, true, true, Us, wce, s_star);

load('2005_diag.mat')
logf_interrs = logf_interrs';
s_star = s_star';
[post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, true, Us, wce, s_star);
    
load('1970_diag_GH.mat')
logf_interrs = logf_interrs';
[~, w, wce, alph] = diag_calib(gam, lambda_GH, d, v, Us_GH, s_star_GH);
[post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, true, Us_GH, wce, s_star_GH);
    
load('2005_diag_GH.mat')
logf_interrs = logf_interrs';
[post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, w, true, true, Us_GH, wce, s_star_GH);
    
    