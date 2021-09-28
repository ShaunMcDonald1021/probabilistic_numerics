% This script repeatedly creates a 2nd-order sparse Gauss-Hermite grid (as
% in Karvonen et al. 2018) and saves it to a file. The file is ultimately
% used in fisheries_demo.R. The 100 repetitions allow us to more accurately
% assess the variability in computation time for the diagnostic with this
% grid. Make sure you download FSKQ prior to running this (see README).
% Note also that this is heavily based on part of demo.m from that repo.

for i = 1:100
    tic
    XS = gh_seq(2);
    us = sparse_gens(XS, d);
    Us = fss_gen(us);
    s_star = cell2mat(Us);
    time = toc;
    save('GH_stuff.mat', 's_star', 'time')
end
