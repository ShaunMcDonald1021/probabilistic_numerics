function [w, kmv] = kqw(X, k, kmean, isotropic)
% KQW - weights of kernel quadrature
%   [w, kmv] = KQW(X, k, kmean) computes the weights of
%   a kernel quadrature rule at points X. Also returned is the vector
%   kmv of kernel mean evaluations.
%
%   The argument isotropic determines if the kernel is to be
%   assumed isotropic or not: 'true' (default) if yes, 'false'
%   if no.
%
% INPUT
%   - X           nodes, one columnes being the nodes
%   - k           kernel, given as k(r) if isotropic or k(x,y) if not
%   - kmean       kernel mean, given as kmean(x)
%   - isotropic   'true' for isotropic (default), 'false' for not isotropic
%
% OUTPUT
%   - w           column vector of the weights
%   - kmv         column vector of the kernel mean evaluations

% Toni Karvonen, 2017

  if ~exist('isotropic', 'var')
    isotropic = 'true';
  else
    if ~strcmp(isotropic, 'false')
      isotropic = 'true';
    end
  end
  
  K = kmat(X, k, isotropic);
  N = size(X, 2);
  kmv = zeros(N, 1);
  for n = 1:N
    kmv(n) = kmean(X(:,n));
  end
  w = K\kmv;

end
