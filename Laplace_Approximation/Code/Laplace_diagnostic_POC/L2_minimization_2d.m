function lambda = L2_minimization_2d(gam, v)
% Minimize (w.r.t. lambda) the "posterior L2 error" of the LA diagnostic
% when applied to a bivariate T density with v degrees of freedom, given
% the parameter $\gamma$ (gam). gam = sqrt(1.5*(v+d)/(v+d-3)) is a sensible
% default choice. Preliminary interrogation grid consists of points at 0,
% +/-1, +/-2, +/-3 along each axis. See Section 5.1 of the manuscript for
% more details.
% Outputs the approximately "optimal" lambda value.

d = 2;

% Test function/T density
tau = @(x) mvtpdf(x, eye(d), v);
% [mode, ~, tau_at_mode, hess_at_mode] = laplace_ingredients(tau, d,
% false);
mode = zeros([1 d]);
tau_at_mode = tau(mode);
hess_at_mode = -diag(repelem((v+d)/v, d));

% Gaussian approximation to tau
gauss_approx = @(x) tau_at_mode*exp(sum((x*hess_at_mode).*x,2)/2);

% Integrating measure
inv_hess_at_mode = -diag(repelem(v/(v+d), d));
% Could just do inv(hess_at_mode) but I hate looking at that warning MATLAB
% gives you
g = @(x) mvnpdf(x, zeros([1 d]), -gam^2*inv_hess_at_mode);

% Prior mean function for GP
mx_0 = @(x) gauss_approx(x)./g(x);

% Preliminary grid used in Section 5.1 of manuscript (eqn. (19))
s_star = [[-3:3 zeros([1 6])]; [zeros([1 7]) -3:-1 1:3]]';
s = s_star/sqrt(-hess_at_mode(1,1));

interrs = (tau(s)./g(s)) - mx_0(s);

% Set up grid for Riemann sum approximation to L2 integral
delta = 0.005;
[xgrid, ygrid] = meshgrid(-8:delta:8, -8:delta:8);
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
clear x % Saves space

% L2 error and its derivative
function [L2, dL2] = L2_with_grad(l)
    weighted_interrs = Cx_0(mah_dist_s, l)\interrs;
    L2_points = prior_diffs + gvals.*...
        (Cx_0(mah_dist_x_s, l)*weighted_interrs);
    L2 = delta^2*sum(L2_points.^2);
    if nargout == 2
        dL2 = delta^2*2*sum(L2_points.*gvals.*(dCx_0(mah_dist_x_s, l) -...
            (Cx_0(mah_dist_x_s, l)*...
            (Cx_0(mah_dist_s,l)\dCx_0(mah_dist_s, l))))*...
            weighted_interrs);
    end
end

% Find optimal lambda
% For gam = sqrt(1.5*(v+d)/(v+d-3) it ends up being 4.2240. The optimizer
% doesn't seem to work with initial values higher than this, but L2 and dL2
% were both larger for the values of l > 4.2240 that I tried
[lambda,~, ~, ~, ~, ~] = fminunc(@(l) L2_with_grad(l), 1,...
    optimoptions('fminunc', 'CheckGradients', true,...
    'FiniteDifferenceType', 'central', 'OptimalityTolerance', 1e-10,...
    'StepTolerance', 1e-9, 'Display', 'iter-detailed',...
    'SpecifyObjectiveGradient', true));
end
