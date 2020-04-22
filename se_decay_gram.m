function gram_mat_simp = se_decay_gram(s,lambda)

% Adapted from RR1d_se.m (Chkrebtii et al.) in 'uqdes' folder
% Used to calculate Gram matrix (covariance matrix at interrogation points)
% based on sq. exp. kernel modified to decay with magnitude (equivalent to
% MVN pdf with 45-degree cotours)

% INPUTS:
% s: a set of n interrogation points of dimension d, stored in d*n matrix
% lambda: the length-scale parameter

% OUTPUT:
% Gram matrix of covariances between interrogation points

[d, n] = size(s);
gram_mat = zeros(n, class(s));

% Matrix for quadratic form (see documentation for details)
quad_form = [1/(4*lambda^2)+1/4 -1/(4*lambda^2);...
    -1/(4*lambda^2) 1/(4*lambda^2)+1/4];

% Get lower triangular part of Gram matrix on log scale
for i = 1:n
    for j = 1:i
        point_pair = [s(:,i) s(:,j)]';
        % Probably most direct way to sum up the quadratic forms we need
        gram_mat(i,j) = -sum(point_pair.*(quad_form*point_pair), 'all');
    end
end

% Copy everything from lower triangular part to above diagonal
gram_mat = tril(gram_mat, -1)' + gram_mat;

gram_mat = exp(gram_mat);

% Maintain prettiness and accuracy for as long as possible with sym
if isa(s, 'sym')
    gram_mat_simp = simplify((sqrt(sym(pi))*lambda)^d*gram_mat);
else
    gram_mat_simp = (sqrt(pi)*lambda)^d*gram_mat;
end

end
