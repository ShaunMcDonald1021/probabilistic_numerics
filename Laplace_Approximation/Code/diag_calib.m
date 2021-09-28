function [post_mean, w, wce, alph] =...
    diag_calib(gam, lambda, d, v, Us, s_star)

% Calibrate the LA diagnostic in d dimensions, using a multivariate T
% density as the test function. Make sure you download FSKQ prior to
% running this (see README).

% INPUTS:
% gam, lambda: values of hyperparameters $\gamma$ and $\lambda$, resp.
% d: dimensionality
% v: degrees of freedom for multivariate T density
% Us: The preliminary grid as a cellular array. Will usually be obtained by
% running fss_gen (from FSKQ) on a matrix of generator vectors (see demo.m
% in FSKQ).
% s_star: the preliminary grid in the form of an n*d matrix. Not necessary,
% but it'll save time if you precomputed it

% OUTPUTS:
% post_mean: posterior integral mean for diagnostic applied to test fn
% w: BQ weights. Since these only depend on hyperparameters and s_star,
% they can be used for any function to save time
% wce: worst-case error, or unscaled posterior integral standard deviation.
% wce^2*[function-specific factors]/alph^d = posterior integral variance.
% alph: the precision hyperparameter, calculated to put the test fn on the
% boundary of rejection

% First, make sure fskq is on path
if exist('fskq', 'dir') == 0
    addpath('fskq')
end

% Calculate s_star if it wasn't supplied already
if nargin == 5
    s_star = cell2mat(Us)';
end

logf_inter = log(mvtpdf(s_star.*sqrt(v/(v+d)), eye(d), v));
logf_at_mode = (gammaln((v+d)/2) - gammaln(v/2) - d*log(v*pi)/2);
logdet_T = d*log(v/(v+d))/2;
lap_app = exp(logdet_T + logf_at_mode + d*log(2*pi)/2);

% Re-weight interrogations and subtract prior mean interrogation values
Y = exp(d*log(gam) + logdet_T + logf_inter - log(mvnpdf(s_star/gam)))-...
    exp(d*log(2*pi*gam) + logdet_T + logf_at_mode +...
    log(mvnpdf(sqrt(gam^2-1)*s_star/gam)));

% Kernel stuff
k = @(r)exp(-r.^2/(2*lambda^2));
kmean = @(x)(lambda^2/(gam^2+lambda^2))^(d/2)*...
    exp(-norm(x)^2/(2*(gam^2+lambda^2)));
Ikmean = (lambda^2/(2*gam^2 + lambda^2))^(d/2);

[Q, wce, w] = kq_fss(Y, Us, k, kmean, Ikmean);
post_mean = Q + lap_app;

alph = (Q/(1.96*wce*exp(logdet_T + logf_at_mode)))^(-2/d);

