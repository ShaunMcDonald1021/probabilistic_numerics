function l = eq_fit(X, Y, k, isotropic, l0, kdl)
% EQ_FIT - isotropic length-scale optimization with empirical Bayes
%   Given an isotropic kernel k, its derivative kdl (as a function 
%   kdl(r,l) of distance r of points and the length-scale l), data 
%   locations X (each column represents one point) and real data 
%   Y, l = EQ_FIT(X, Y, k, isotropic, l0, kdl) uses empirical Bayes 
%   to select the  optimal length-scale l.
%
%   With the optional paramter l0, the root finding is started at 
%   l0. The default value is l0 = 1.
%
% INPUT
%   - X           input points
%   - Y           function evaluations
%   - k           an isotropic kernel (without length-scale)
%   - isotropic   'true' or 'false' depending on the kernel being isotropic
%   - l0          initial length-scale (optional), default l0 = 1
%   - kdl         derivative of the kernel w.r.t. length-scale
%
% OUTPUT
%   - l           length-scale fitted via empirical Bayes

% Toni Karvonen, 2017

  if ~exist('l0', 'var')
    l0 = 1;
  end

  if ~exist('isotropic', 'var')
    isotropic = 'true';
  else
    if ~strcmp(isotropic, 'false')
      isotropic = 'true';
    end
  end
  
  Y = Y(:);
  lmld = @(l) logmld(l, k, isotropic, kdl, X, Y);
  l = fzero(lmld, l0, optimset('TolX', 1e-10));
  
end

function lmld = logmld(l, k, isotropic, kdl, X, Y)
  K = kmat(X/l, k, isotropic);
  Kdl = kmat(X, @(r) kdl(r,l), isotropic);
  %a = K\*Y;
  lmld = sum(sum(((K\Y*Y' - eye(size(K,1)))/K).* Kdl', 2));
end
