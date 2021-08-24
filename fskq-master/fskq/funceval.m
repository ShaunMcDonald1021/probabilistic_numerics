function Y = funceval(f, X)
% FUNCEVAL - evaluate a function at multiple vectors
%   Given a function f:R^n -> R^m and a (n x N) matrix
%   X of input vectors, Y = funceval(f, X) returns a
%   (m x N) matrix Y of function evaluations.
%
% INPUT
%   - f   the function f:R^n -> R^m
%   - X   input points, each column is one point
%
% OUTPUT
%   - Y   function evaluation at X, each column is one evaluation

% Toni Karvonen, 2017

  [n, N] = size(X);
  m = size(f(ones(n,1)), 1);
  Y = zeros(m, N);
  
  for i = 1:N
    Y(:,i) = f(X(:,i));
  end
  
end
