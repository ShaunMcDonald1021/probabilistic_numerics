function S = fss_rsmat(Us, k, isotropic)
% FSS_RSMAT - compute the FSS row sum matrix
%   Given an isotropic kernel k and a cell array
%   Us of fully symmetric sets (each column in the
%   cell array elements being one point),
%   S = REDUCE_KMAT(Us, k) computes the row sum
%   matrix Kr.
%
%   The argument isotropic determines if the kernel is to be
%   assumed isotropic or not: 'true' (default) if yes, 'false'
%   if no.
%
% INPUT
%   - Us          a cell array of FSSs, columns being the points
%   - k           kernel, given as k(r) if isotropic or k(x,y) if not
%   - isotropic   'true' for isotropic (default), 'false' for not isotropic
%
% OUTPUT
%   - S           the FSS reduced row sum matrix

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
  % different fully symmetric sets
  S = zeros(J, J);
  
  % Do kernel evaluations differently based on the kernel
  % being isotropic or not.
  if strcmp(isotropic, 'true')
  
    for i = 1:J
      x = Us{i}(:,1);
      for j = 1:J
        U = Us{j};
        S(i,j) = sum( k( sqrt(sum( (repmat(x,1,size(U,2))-U).^2)) ) );
      end
    end
    
  else
  
    for i = 1:J
      x = Us{i}(:,1);
      for j = 1:J
        U = Us{j};
        % This is way slower than what is done with isotropic kernels
        foo = 0;
        for p = 1:size(U,2)
          foo = foo + k(x, U(:,p));
        end
        S(i,j) = foo;
      end
    end 
  
  end

end
  
