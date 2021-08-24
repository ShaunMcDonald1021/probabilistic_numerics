function [Q, wce, wr] = kq_fss(Y, Us, k, kmean, Ikmean, isotropic)
% KQUAD_FSS - fully symmetric kernel quadrature
%   Given function evaluations Y at a cell array Us of fully symmetric sets,
%   an isotropic kernel k, its kernel mean kmean and the integral Ikmean of the
%   kernel mean, [Q, V, w] = kquad_fss(f, Us, k, kmean, Ikmean)
%   computes the fully symmetric kernel quadrature integral estimate Q, 
%   the associated RKHS worst-case error wce and also returns the FSSKQ
%   unique weight vector wr.
%
%   The argument isotropic determines if the kernel is to be
%   assumed isotropic or not: 'true' (default) if yes, 'false'
%   if no.
%
% INPUT
%   - Y           a row vector of function evaluation at points of Us
%   - Us          a cell array of FSSs with columns being the nodes
%   - k           kernel, given as k(r) if isotropic or k(x,y) if not
%   - kmean       kernel mean, given as kmean(x)
%   - Ikmean      integrated kernel mean (i.e. initial WCE)
%   - isotropic   'true' for isotropic (default), 'false' for not isotropc
%
% OUTPUT
%   - Q           FSSKQ integral approximation
%   - wce         FSSKQ RKHS worst-case error (i.e. standard deviation)
%   - wr          column vector of the unique weights

% Toni Karvonen, 2017

  if ~exist('isotropic', 'var')
    isotropic = 'true';
  else
    if ~strcmp(isotropic, 'false')
      isotropic = 'true';
    end
  end

  % Generate the unique weights and kernel mean
  % evaluations
  [wr, kmvr] = kqw_fss(Us, k, kmean, isotropic);
  
  % Compute the FSSKQ approximation
  Q = 0;
  J = length(wr);
  Ls = zeros(J,1);
  ind = 0;
  for i = 1:J
    Ls(i) = size(Us{i}, 2);
    Q = Q + wr(i) * sum(Y(ind+1:ind+Ls(i)));
    ind = ind + Ls(i);
  end
  
  % Compute the worst-case error
  wce = sqrt(Ikmean - wr' * (Ls .* kmvr));
  
end
