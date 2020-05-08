function [mean, covar, true_int, lap_approx] ...
    = laplace_poc_sym(f, a, b, d, s, lambda, alpha, do_sym_calcs, true_int)
% THIS FUNCTION IS CURRENTLY UNFINISHED. IT WILL EVENTUALLY BE AN
% ADAPTATION OF LAPLACE_POC.M WITH FULLY SYMBOLIC COMPUTATIONS, BUT THIS IS
% NOT CURRENTLY A PRIORITY.

%Covariance stuff
C0_t = se_decay_gram(s, lambda);
C0_cross = se_decay_cross_covar(s, lambda, a, b);
if isa(s, 'sym') && ~do_sym_calcs %do_sym_calcs will make things VERY slow
    C0_t = double(C0_t);
    C0_cross = double(C0_cross);
end
C0_multiple = (4*pi^(3/2)*lambda^2/(lambda^2 + 2)^(1/2))^d;
if a == -Inf && b == Inf
    C0 = C0_multiple;
else
    C_sigma = 2*[lambda^2+1 1; 1 lambda^2+1]/(lambda^2+2);
    if ~do_sym_calcs
        C0 = C0_multiple*mvncdf([a a], [b b], [0 0], C_sigma)^d;
    else
        syms x real;
        C0 = C0_multiple*...
            vpaintegral(se_decay_cross_covar(x, lambda, a, b), x, a, b)
    end
end

weights = C0_cross/C0_t;
post_var = C0 - weights*C0_cross';

[mode, f_sym, f_at_mode, hess_at_mode] = laplace_ingredients(f, d);

if isa(s, 'gpuArray')
    hess_at_mode = gpuArray(double(hess_at_mode));
else
    hess_at_mode = double(hess_at_mode);
end %TODO: add routine for symbolic computations (which DO NOT play nice w/GPU)

[V, D] = eig(hess_at_mode);
rot_scale_mat = sqrt(-D)\V';

s_trans = s'*rot_scale_mat + double(mode);

diff_vals = f(s_trans) - double(f_at_mode)*exp(-vecnorm(s).^2/2)';

m0 = double(f_at_mode)*(2^(1/2)*pi^(1/2)*(erf(2^(1/2)*b/2) - erf(2^(1/2)*a/2))/2)^d;

post_mean = m0 + weights*diff_vals;

