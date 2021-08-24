function [Q, wce, w] = kq(Y, X, k, kmean, Ikmean, isotropic)
% KQ - compute an integral with kernel quadrature
%   Given function evaluations Y at points X (each column being one point),
%   a kernel k, its kernel mean kmean and the integral Ikmean of the
%   kernel mean, [Q, V, w] = kquad_fss(f, Us, k, kmean, Ikmean)
%   computes the kernel quadrature integral estimate Q, 
%   the associated RKHS worst-case error wce and also returns the weights w.
%
%   The argument isotropic determines if the kernel is to be
%   assumed isotropic or not: 'true' (default) if yes, 'false'
%   if no.
%
% INPUT
%   - Y           a row vector of function evaluation at points of Us
%   - X           the nodes
%   - k           kernel, given as k(r) if isotropic or k(x,y) if not
%   - kmean       kernel mean, given as kmean(x)
%   - Ikmean      integrated kernel mean (i.e. initial WCE)
%   - isotropic   'true' for isotropic (default), 'false' for not isotropic
%
% OUTPUT
%   - Q           kernel quadrature integral approximation
%   - wce         RKHS worst-case error (i.e. standard deviation)
%   - w           column vector of the weights

% Toni Karvonen, 2017

  if ~exist('isotropic', 'var')
    isotropic = 'true';
  else
    if ~strcmp(isotropic, 'false')
      isotropic = 'true';
    end
  end
  
  [w, kmv] = kqw(X, k, kmean, isotropic);
  Q = w' * Y(:);
  wce = sqrt(Ikmean - w' * kmv);
  
end
