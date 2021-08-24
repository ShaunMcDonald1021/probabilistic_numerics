function CX = cartesian(XS)
% CARTESIAN - a Cartesian product
%   Given a cell array XS of vectors this function
%   produces the Cartesian product CX of these vectors.
%
% INPUT
%   - XS    a cell array of vectors that may be of
%           different lengths
%
% OUTPUT
%   - CX    the Cartesian product, columns are the points

% Simo Särkkä & Toni Karvonen, 2017
  
  % Dimension of the resulting Cartesian product set
  d = length(XS);
  
  % Lengths of each 1D point set
  Ls = [];
  for l = 1:d
    Ls(l) = length(XS{l});
  end
  L = prod(Ls);
  
  % Form the index set
  ind = zeros(d,L);
  num = 0:(L-1);
  for n = 1:d
    ind(n,:) = rem(num,Ls(n))+1;
    num = floor(num / Ls(n));
  end
  
  % Form the actual Cartesian product
  CX = zeros(size(ind));
  for n = 1:d
    CX(n,:) = XS{n}(ind(n,:));
  end
  
end
