function us = sparse_gens(XS, d)
% SPARSE_GENS - produce generators generating a sparse grid
%   Given a cell array XS of distinct 1D point sets (with XS{1} = 0),
%   us = SPARSE_GENS(XS, d) produces the generator vectors 
%   us that generate the corresponding sparse grid in dimension d. 
%   The generators are given in a matrix where each column
%   corresponds to one generator vector.
%
% INPUT
%   - XS    a cell array of distinct 1D point sets (XS{1} = 0)    
%   - d     the dimension
%
% OUTPUT
%   - us    vectors generating the sparse grid, columns being the vectors

% Toni Karvonen, 2017
  
  % q is one less than the number of 1D point sets
  q = length(XS) - 1;
  
  % Generate the FSS level sequence
  seq = levelseq_fss(q, d);
  
  % Generate a set of index vectors. Each index i
  % represents elements in XS{i}.
  inds = [];
  for i = 1:size(seq,1)
    ZS = {};
    for j = 1:d
      ZS{j} = 1:seq(i,j);
    end
    inds = [inds cartesian(ZS)];
  end
  
  % We only want unique index vectors
  % THIS IS PROBABLY NOT THE SMARTEST WAY OF DOING THIS
  inds = unique(sort(inds,'descend')','rows')';
  
  % Generate the actual generators based on the 1D points
  us = [];
  for j = 1:size(inds,2)
    ZS = {};
    for i = 1:d
      Z = XS{inds(i,j)};
      Z = Z(Z >= 0);
      ZS{i} = Z;
    end
    us = [us cartesian(ZS)];
  end
  
  % We only want unique generator vectors
  % THIS IS PROBABLY NOT THE SMARTEST WAY OF DOING THIS
  us = unique(sort(us,'descend')','rows')';
  
end
