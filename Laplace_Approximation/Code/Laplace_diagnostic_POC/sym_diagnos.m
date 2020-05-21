function [lambda, gam, C0, m0, bias, var_corr] = sym_diagnos(d, v, is_indep)

% This function allows us to assess the effects of the covariance kernel
% hyperparameters on the Laplace diagnostic, showing how they affect bias
% and variance when the full function is a multivariate T distribution.
% Everything here uses the symbolic toolkit (i.e. slow).

% INPUTS:
% d: the dimension. MUST be 1 or 2, since higher dimensions would likely
% have too many interrogation points for symbolic inversion to be viable.
% TODO: add error message for d > 2
% v: degrees of freedom for the d-dimensional T-distribution.
% is_indep: logical only needed for d == 2. Determines if we use a product
% of 1-dimensional pdf's (independent), or the actual bivariate
% T-distribution (uncorrelated but NOT independent)

% OUTPUTS
% lambda, gam: symbolic variables corresponding to covariance length-scale
% and decay, respectively.
% C0, m0: symbolic functions of lambda and gam corresponding to prior
% variance and mean of integral (over whole domain), respectively. 
% The latter is, of course, the Laplace approximation.
% bias: symfun for the difference between the prior and posterior means
% (i.e. weighted sum of f(s) - m0_t(s)).
% var_corr: symfun for the variance correction (i.e. C0 - var_corr =
% posterior variance).

syms lambda gam real
assumeAlso(lambda > 0 & gam > 0);

if d == 1
    % Interrogation grid is +-2 standard deviations
    s = sym(-3:1.5:3);
    f = @(x) gamma(sym((v+1)/2))*(1+x.^2/v).^(-(v+1)/2)/...
        (gamma(sym(v/2))*sqrt(v*sym(pi)));
    
elseif d == 2
    % Cross-shaped grid: +-2 s.d.'s along each axis
    s = sym(zeros([2 9]));
    s(1,1:5) = [-3 -1.5];
    s(2, 8:9) = [1.5 3];
    
    if is_indep
        % Product of Cauchy distributions...
        f = @(x) gamma(sym((v+1)/2))^2*((1+x(:,1).^2/v).*...
            (1+x(:,2).^2/v)).^(-(v+1)/2)/(gamma(sym(v/2))^2*v*sym(pi));
    
    else
        % ...or actual bivariate T
        f = @(x) gamma(sym((v+2)/2))*(1+(x(:,1).^2 + x(:,2).^2)/v).^(-(v+2)/2)/...
            (gamma(sym(v/2))*v*sym(pi));
    end
end

[mode, ~, f_at_mode, hess_at_mode, ~] = laplace_ingredients(f, d,...
    false);
[V, D] = eig(hess_at_mode);
rot_scale_mat = sqrt(-D)\V';
rot_det = abs(det(rot_scale_mat));

s_trans = s'*rot_scale_mat + mode;
m0_t = @(x) f_at_mode*exp(sum(((x-mode)*hess_at_mode).*(x-mode), 2)/2);
m0 = rot_det*(f_at_mode*...
    (sym(2)^(1/2)*sym(pi)^(1/2))^d);
gauss_interr = m0_t(s_trans);
f_interr = f(s_trans);

% TODO: peel out variance stuff, since this doesn't really depend on the
% function values. It would be faster to save this stuff separately and
% then scale as needed to more quickly run diagnostics for many different
% values of v
C0_t = se_decay_gram(s, lambda, gam);
gram_inv = inv(C0_t);
C0_cross = se_decay_cross_covar(s, lambda, gam, -Inf, Inf);
C0 = rot_det^2*((4*gam^2*lambda^2*sym(pi)^(3/2))/...
    (2*gam^2 + lambda^2)^(1/2))^d;

weights = C0_cross*gram_inv;
var_corr = rot_det^2*weights*C0_cross';
bias = rot_det*weights*(f_interr - gauss_interr);

var_corr = simplify(var_corr);
bias = simplify(bias);
end
