function seq = levelseq_fss(q, d)
% LEVELSEQ_FSS - level sequence for fully symmetric methods
%   Given a sparse grid level q and dimension d,
%   seq = LEVELSEQ_FSS(q, d) generates all possible
%   sequences of levels needed to form generators
%   generating the sparse grid.
%
% INPUT
%   - q     level of the sparse grid
%   - d     dimension
%
% OUTPUT
%   - seq   sequence of levels

% Toni Karvonen, 2017
  
  % Let spinterp generate the level sequence
  seq = spgetseq(q, d);
  m = size(seq, 1);
  seq = seq + ones(m, d, 'uint8'); % uint8 is required due to spinterp
  
  % Remove those not in ascending order
  seq_desc = [];
  for i = 1:m
    if issorted(seq(i,:))
      seq_desc = [seq_desc; seq(i,:)];
    end
  end
  
  % Flip to get in ascending order
  seq = fliplr(seq_desc);
  
end
