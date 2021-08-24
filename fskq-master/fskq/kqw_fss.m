function [wr, kmvr] = kqw_fss(Us, k, kmean, isotropic)
% KQW_FSS - unique weights of fully symmetric kernel quadrature
%   [wr, kmvr] = KQW_FSS(Us, k, kmean) computes the unique weights,
%   one for each fully symmetric set in the cell array Us of
%   fully symmetric sets constituting the whole node set,
%   of fully symmetric kernel quadrature with the kernel
%   that has the kernel mean kmean. Also returned is the vector
%   kmvr of kernel mean evaluations at the norm of each FSS.
%
%   The argument isotropic determines if the kernel is to be
%   assumed isotropic or not: 'true' (default) if yes, 'false'
%   if no.
%
% INPUT
%   - Us          a cell array of FSSs with columns being the nodes
%   - k           kernel, given as k(r) if isotropic or k(x,y) if not
%   - kmean       kernel mean, given as kmean(x)
%   - isotropic   'true' for isotropic (default), 'false' for not isotropic
%
% OUTPUT
%   - wr          column vector of the unique weights
%   - kmvr        column vector of the unique kernel mean evaluations

% Toni Karvonen, 2017
  
  if ~exist('isotropic', 'var')
    isotropic = 'true';
  else
    if ~strcmp(isotropic, 'false')
      isotropic = 'true';
    end
  end
  
  J = length(Us);
  
  % Form the matrix S containing row sums between
  % different fully symmetric sets and the vector
  % kmvr of kernel mean evaluations at norms of 
  % each generator.
  S    = fss_rsmat(Us, k, isotropic);
  kmvr = zeros(J, 1);
  for j = 1:J
    kmvr(j) = kmean(Us{j}(:,1));
  end
  
  % Solve for the set of distinct weights
  wr = S\kmvr;

end
