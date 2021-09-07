d = 72;
v = 25921;
gam = sqrt(1.5*(v+d)/(v+d-3));
lambda = 3.7;

load('1970_diag.mat')
s_star = s_star';
logf_interrs = logf_interrs';
Us = fss_gen(s_star([1 73],:)');
[~, w, wce, alph] = diag_calib(gam, lambda, d, v, Us, s_star);

[post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode, log_T_det,...
    d, gam, alph, w, true, true, Us, wce, s_star);