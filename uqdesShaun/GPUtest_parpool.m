%% Analysis of GPU speed using a parallel pool
% This code does a similar thing to benchmark_test.m, but this is an old
% version designed to split the iterations across a parallel pool of GPU's.
% Intended for something like Dave's remote machine with two GPU's. Rather
% than saving a figure, it saves the computation times themselves for later
% use. Then you can modify benchmark_test.m to load them and rerun the test
% on the CPU for comparison.

clear all; close all; clc;

num_iter = 30; %number of iterations

% Load num_iter*2 matrix of initial conditions, assuming this has been
% premade. Each row is a set of initial conditions for one of the
% iterations.
load('u0mat.mat')

sspan = [0 10];  theta = 2; 
odefn = @simpleode; odesoln = @simpleode_solution;
%nsolvesvec = [50 100 200 500]; 
% nsolves doesn't seem to matter too much anyway,and we're going to VERY
% HIGH N so I'd rather not do more repetitions than I have to
nsolvesvec = 100; kernel = 'sqexp';  % sqexp or uniform
Nvec = [50 100 200 500 1000 2000 5000 10000];

logGPUtimepar = zeros(15, length(Nvec), length(nsolvesvec)); 
% 15 iterations on each GPU, then we'll concatenate the results after

state = 2;
parpool('local', 2) %Set up parpool

spmd
    gpuDevice(labindex);
    % Running it backwards should avoid memory fragmentation issues
    for nind = length(Nvec):-1:1 
        N = Nvec(nind); 
        ds = range(sspan)/(N-1);
        lambda = 1*ds; alpha = N/100;
        for iterind = 1:15
            u0 = u0mat(15*(labindex-1) + iterind,:);

            for solveind = 1:length(nsolvesvec)
                nsolves = nsolvesvec(solveind)
                % Progress message may be weird and not display unless we allow
                % nsolves to print here as well
                tic

                [~,~,logGPUtimepar(iterind, nind, solveind),~] = ...
                    uqdes_mod(sspan, nsolves, N, kernel, lambda, alpha, odefn,...
                    u0, theta);

                % Progress message
                disp(strcat("N = ", num2str(N), ", nsolves = ",...
                    num2str(nsolves), ", iteration ", ...
                    num2str(15*(labindex-1)+iterind), " took ", ...
                    num2str(toc), " seconds."))

            end %end solveind loop
        end %end iterind loop
        logGPUtimegat = gcat(logGPUtimepar, 1);
    end %end nind loop
end %end spmd

logGPUtime = logGPUtimegat{1};
save('logGPUtime.mat', 'logGPUtime');
