function [Us Ls] = fss_gen(us, d)
% FSS_GEN - generate a number of FSSs in a cell array
%   Given a matrix us of generator vectors, each column
%   representing one such vector [Us Ls] = FSS_GEN(us) 
%   generates a collection of fully symmetric sets based on these generators
%   in dimension Each set is saved in the cell array Us. Also saved in Ls are
%   numbers of elements in each of the fully symmetric sets.
%
%   If also the dimension d is given, [Us Ls] = FSS_GEN(us, d) 
%   generates the sets in that dimension.
%
% INPUT
%   - us  the generator vectors, each column being one such vector
%   - d   dimension (optional), if not set, assumed equal to size(us,1)      
%
% OUTPUT
%   - Us  the fully symmetric sets in a cell array, each element
%         corresponding to one generator
%   - Ls  cardinalities of the FSSs in a column vector

% Toni Karvonen, 2017

  % Some argument checks
  if nargin < 2
    d = size(us, 1);
  end
  
  if size(us, 1) > d
    error('Length of the generator vector u cannot exceed the given dimension d.');
  end
  
  % Generate the FSSs
  u_num = size(us, 2);
  Us = {};
  Ls = zeros(u_num, 1);
  for i = 1:u_num
    U = fss(us(:,i), d);
    Us{i} = U;
    Ls(i) = size(U, 2);
  end

end
