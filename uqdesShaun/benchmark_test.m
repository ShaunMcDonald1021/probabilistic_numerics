%% Benchmarking of GPU and/or sparsity using the simple ODE from Figure 1 of Oksana's paper
% We'll run 30 iterations of both the CPU and GPU solvers at various N-values
% with randomly-generated initial conditions.
% We will also loop over different values of nsolves.
% It doesn't seem to make much of a difference here (and theoretically
% there isn't much reason to suspect that it would), but it doesn't hurt
% to check anyways.

clear all; close all; clc;

num_iter = 3; %number of iterations
sspan = [0 10];
theta = 2;
odefn = @simpleode; odesoln = @simpleode_solution;
nsolvesvec = [50 100 200 500]; kernel = 'uniform';  % sqexp or uniform
% (Although only uniform is relevant if we want to see how it performs with
% sparse matrices)
Nvec = [50 100 200 500 1000 2000 5000 10000 20000 50000];
logCPUtime = zeros(num_iter, length(Nvec), length(nsolvesvec));
logGPUtime = zeros(num_iter, length(Nvec), length(nsolvesvec)); 
logCPUSparsetime = zeros(num_iter, length(Nvec), length(nsolvesvec));
logGPUSparsetime = zeros(num_iter, length(Nvec), length(nsolvesvec)); 

state = 2;
for iterind = 1:num_iter
u0 = [unifrnd(-5, 5), unifrnd(-5, 5)];

    for nind = 1:length(Nvec)
        N = Nvec(nind);
        ds = range(sspan)/(N-1);
        lambda = 1*ds; alpha = N/100;

        for solveind = 1:length(nsolvesvec)
            nsolves = nsolvesvec(nind);

            % We're not really interested in the actual results here. Prelminary
            % experiments showed that the GPU should be consistent with the CPU.
            % Only time matters here.
            [~,~,logCPUtime(iterind, nind, solveind),~] = uqde_mod(sspan,...
                nsolves,N,kernel,lambda,alpha,odefn,u0,theta, false, false);

            [~,~,logGPUtime(iterind, nind, solveind),~] = uqdes_mod(sspan,...
                nsolves, N, kernel, lambda, alpha, odefn, u0, theta,...
                true, false);
            
            [~,~,logCPUSparsetime(iterind, nind, solveind),~] = uqdes_mod(sspan,...
                nsolves, N, kernel, lambda, alpha, odefn, u0, theta,...
                false, true);
            
            [~,~,logGPUSparsetime(iterind, nind, solveind),~] = uqdes_mod(sspan,...
                nsolves, N, kernel, lambda, alpha, odefn, u0, theta,...
                true, true);

            % The very first run is always slower.
            % This can't be the initalization overhead, because we start measuring time inside the functions.
            % What is the system doing differently on the first iteration?
        end
    end
end

% Plot results
figure
for solveind = 1:length(nsolvesvec)
    subplot(1,length(nsolvesvec), solveind)
    loglog(Nvec, median(logCPUtime(:,:,solveind),1), 'r--o',...
        Nvec, median(logGPUtime(:,:,solveind),1), 'b--o', ...
        Nvec, median(logCPUSparsetime(:,:,solveind),1), 'g--o', ...
        Nvec, median(logGPUSparsetime(:,:,solveind),1), 'k--o')
    
    title(strcat("numsolves = ", num2str(nsolvesvec(solveind))));
    xlabel('N')
    ylabel('Computation time (seconds)')
    legend('CPU', 'GPU', 'CPU Sparse', 'GPU Sparse', 'Location', 'southeast')
end

% Save the figure. You may need to change this path, depending on where you
% run this code from.
saveas(gcf, 'probabilistic_numerics/uqdesShaun/benchmark_figs/benchmark_test.png')
