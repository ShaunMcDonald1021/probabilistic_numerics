function K = kmat(X, k, isotropic)
% KMAT - kernel matrix
%   Constructs the kernel matrix given a kernel
%   k and a set of points X. Each column of X
%   represents one point.
%
%   The argument isotropic determines if the kernel is to be
%   assumed isotropic or not: 'true' (default) if yes, 'false'
%   if no.
%
% INPUT
%   - X           the set of points, each column being one point
%   - k           kernel, given as k(r) if isotropic or k(x,y) if not
%   - isotropic   'true' for isotropic (default), 'false' for not isotropc
%
% OUTPUT
%   - K           the kernel matrix

% Toni Karvonen, 2017

  if ~exist('isotropic', 'var')
    isotropic = 'true';
  else
    if ~strcmp(isotropic, 'false')
      isotropic = 'true';
    end
  end
  
  if strcmp(isotropic, 'true')
    D = squareform(pdist(X'));
    K = k(D);  
  else
    N = size(X, 2);
    K = zeros(N, N);
    for i = 1:N
      for j = 1:i
        foo = k(X(:,i), X(:,j));
        K(i,j) = foo;
        K(j,i) = foo;
      end
    end
  end

end
