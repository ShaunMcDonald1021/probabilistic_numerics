function XS = cc_seq(q)
% CC_SEQ - Clenshaw-Curtis point sequence 
%   Given a level q, XS = CC_SEQ(q) generates nested Clenshaw-Curtis
%   points in the cell array XS. Only additional non-negative points
%   are saved. That is, each X{i} is distinct and contains only
%   non-negative part of the Clenshaw-Curtis point set. With level q,
%   the cell array XS has q+1 elements, the first being just 0.
%
% INPUT
%   - q     sparse grid level
%
% OUTPUT
%   - XS    cell array of distinct Clenshaw-Curtis point sets

% Toni Karvonen, 2017
  
  n = q + 1;
  XS = {};
  XS{1} = 0;
  XS{2} = 1;
  for i = 3:n
    m = 2^(i-1) + 1;
    XS{i} = cos( ( pi*((2:2:(m-1)/2)-1)  ) / (m-1) );
  end
  
end
