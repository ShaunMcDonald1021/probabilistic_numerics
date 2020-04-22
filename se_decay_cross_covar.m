function cross_vals = se_decay_cross_covar(y, lambda, a, b)

% Cross-covariance (i.e. between full and integrated/target function)
% based on sq. exp. kernel modified to decay with magnitude (equivalent to
% MVN pdf with 45-degree cotours)

% INPUTS:
% y: matrix of evaluation points in the non-integrated argument. Should be
% a d*n matrix: n points of dimension d.
% lambda: the length-scale parameter
% a, b: integration limits (possibly infinite) for integrated argument.
% Currently, all d dimensions of integrated argument must have same limits.
% Should be fine WLOG by scaling the axes of the target function.

% OUTPUT:
% Cross covariance C([b...b], y), where first arg is a d-dimensional copy
% of b and we evaluate vectorwise for all rows of y

% Want to use sym(pi) for symbolic vectors to maximize accuracy
if isa(y, 'sym')
    cross_val_mat = -(lambda^2*sym(pi)*exp(-(y.^2.*(lambda^2 + 2))./...
        (4*(lambda^2 + 1))).*(erf((a*lambda^2 + a - y)./...
        (2*lambda*(lambda^2 + 1)^(1/2))) - erf((b*lambda^2 + b - y)./...
        (2*lambda*(lambda^2 + 1)^(1/2)))))./(lambda^2 + 1)^(1/2);

    % Make it nice and pretty
    cross_vals = simplify(prod(cross_val_mat));
else
    cross_val_mat = -(lambda^2*pi*exp(-(y.^2.*(lambda^2 + 2))./...
        (4*(lambda^2 + 1))).*(erf((a*lambda^2 + a - y)./...
        (2*lambda*(lambda^2 + 1)^(1/2))) - erf((b*lambda^2 + b - y)./...
        (2*lambda*(lambda^2 + 1)^(1/2)))))./(lambda^2 + 1)^(1/2);

    cross_vals = prod(cross_val_mat, 1);
end

end
