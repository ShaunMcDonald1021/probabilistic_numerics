function cross_vals = se_decay_cross_covar(x, lambda, gamma, a, b)

% Cross-covariance (i.e. between full and integrated/target function)
% based on sq. exp. kernel modified to decay with magnitude (equivalent to
% MVN pdf with 45-degree cotours)

% INPUTS:
% x: matrix of evaluation points in the non-integrated argument. Should be
% a d*n matrix: n points of dimension d.
% lambda: the length-scale parameter
% a, b: integration limits (possibly infinite) for integrated argument.
% Currently, all d dimensions of integrated argument must have same limits.
% Should be fine WLOG by scaling the axes of the target function.

% OUTPUT:
% Cross covariance C([b...b], x), where first arg is a d-dimensional copy
% of b and we evaluate vectorwise for all rows of x

% Want to use sym(pi) for symbolic vectors to maximize accuracy
if isa(x, 'sym')
    cross_val_mat = -(gamma*lambda^2*sym(pi)*exp(-(x.^2*(2*gamma^2 + lambda^2))./...
        (4*gamma^2*(gamma^2 + lambda^2))).*(erf((a*gamma^2 + a*lambda^2 ...
    - gamma^2*x)./(2*gamma*lambda*(gamma^2 + lambda^2)^(1/2))) - ...
    erf((b*gamma^2 + b*lambda^2 - gamma^2*x)./(2*gamma*lambda*...
    (gamma^2 + lambda^2)^(1/2)))))./(gamma^2 + lambda^2)^(1/2);

    % Make it nice and pretty
    cross_vals = simplify(prod(cross_val_mat, 1));
else
    if gamma == Inf
        cross_val_mat = -lambda^2*pi*(erf((a - x)./(2*lambda)) - erf((b - x)./(2*lambda)));
    else
        cross_val_mat = -(gamma*lambda^2*pi*exp(-(x.^2*(2*gamma^2 + lambda^2))./...
            (4*gamma^2*(gamma^2 + lambda^2))).*(erf((a*gamma^2 + a*lambda^2 ...
        - gamma^2*x)./(2*gamma*lambda*(gamma^2 + lambda^2)^(1/2))) - ...
        erf((b*gamma^2 + b*lambda^2 - gamma^2*x)./(2*gamma*lambda*...
        (gamma^2 + lambda^2)^(1/2)))))./(gamma^2 + lambda^2)^(1/2);
    end
    
    cross_vals = prod(cross_val_mat, 1);
end

end
