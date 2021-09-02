
function lambda = L2_minimization_2d(gam, v)
% Minimize (w.r.t. lambda) the "posterior L2 error" of the LA diagnostic
% when applied to a bivariate T density with v degrees of freedom, given
% the parameter $\gamma$ (gam). gam = sqrt(1.5*(v+d)/(v+d-3)) is a sensible
% default choice. See Section 5.1 of the manuscript for details.
% Outputs the approximately "optimal" lambda value.

d = 2;

tau = @(x) mvtpdf(x, eye(d), v);
tau_at_mode = tau(zeros([1 d]));
hess_at_mode = -diag(repelem((v+d)/v, d));

% Gaussian approximation to tau
gauss_approx = @(x) tau_at_mode*exp(sum((x*hess_at_mode).*x,2)/2);
% Integrating measure
g = @(x) mvnpdf(x, zeros([1 d]), -gam^2*hess_at_mode);
% Prior mean function for GP
mx_0 = @(x) gauss_approx(x)./g(x);

% Preliminary grid used in Section 5.1
s_star = [[-3:3 zeros([1 6])]; [zeros([1 7]) -3:-1 1:3]]';
s = s_star/sqrt(-hess_at_mode(1,1));
interrs = (tau(s)./g(s)) - mx_0(s);

% Set up grid for Riemann sum approximation to L2 integral
[xgrid, ygrid] = meshgrid(-10:0.01:10, -10:0.01:10);
x = [xgrid(:) ygrid(:)];
clear xgrid ygrid % Saves space

% Mahalanobis distance between interrogation points
mah_dist_s = pdist2(s_star, s_star);
% Mahalanobis distance between grid points and interrogation points
mah_dist_x_s = pdist2(x*sqrt(-hess_at_mode(1,1)), s_star);

% These lines more or less borrowed from kq_kernel.m in FSKQ code
% (Karvonen et al. 2018)
Cx_0 = @(r, l) exp(-r.^2/(2*l^2));
dCx_0 = @(r, l) r.^2 .* exp(-r.^2/(2*l.^2))./l.^3;

% Precompute these to save time
gvals = g(x);
prior_diffs =  gauss_approx(x) - tau(x);

% L2 error and its derivative
L2 = @(l) sum((prior_diffs + gvals.*...
    Cx_0(mah_dist_x_s, l)*(Cx_0(mah_dist_s, l)\interrs)).^2);
dL2 = @(l) 2*sum((prior_diffs + gvals.*...
    Cx_0(mah_dist_x_s, l)*(Cx_0(mah_dist_s, l)\interrs)).*gvals.*...
    ((dCx_0(mah_dist_x_s, l) - (Cx_0(mah_dist_x_s, l)/...
    Cx_0(mah_dist_s, l))*dCx_0(mah_dist_s, l))*(Cx_0(mah_dist_s, l)\...
    interrs)));
L2_with_grad = @(y) deal(L2(y),dL2(y));

[lambda,~, ~, ~, ~, ~] = fminunc(L2_with_grad, .1,...
    optimset('GradObj','on', 'TolX', 1e-15, 'Display', 'iter'));

