function XS = gh_seq(q)
% GH_SEQ - Clenshaw-Curtis point sequence 
%   Given a level q, XS = CC_SEQ(q) generates nested Gauss-Hermite
%   points in the cell array XS. These are the 0 and the q positive
%   roots of (2q+1)th Hermite polynomial.
%
% INPUT
%   - q     sparse grid level
%
% OUTPUT
%   - XS    cell array of Gauss-Hermite points

% Toni Karvonen, 2017

  % Nodes are the roots (2q+1)th Hermite polynomial. They
  % can be obtained as the eigenvalues of the tridiagonal
  % Jacobi matrix J constructed from the three-term recurrence
  % relation of Hermite polynomials.
  
  n = 2*q + 1;
  J = zeros(n, n);
  b = sqrt(1:n-1);
  J(n+1:n+1:end-1) = b;
  J(2:n+1:end-n) = b;
  XS = eig(J)';
  XS = num2cell(XS(q+1:end));
  % The eigenvalue 0 is not computed exactly zero
  XS{1} = 0;

end
