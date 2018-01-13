%% Comparison of GPU, GPU using sparse arrays, and CPU using the simple ODE from Figure 1 of Oksana's paper
%We'll run 30 iterations of both the CPU and GPU solvers at various N-values with randomly-generated initial conditions
%We will also loop over different values of nsolves. It doesn't seem to make much of a difference here (and theoretically there isn't much reason to suspect that it would), but it doesn't hurt to check anyways.

clear all; close all; clc;

num_iter = 3; %number of iterations
u0 = [-1 0]; sspan = [0 10]; sspanGPU = gpuArray(sspan); theta = 2; thetaGPU = gpuArray(theta);
odefn = @simpleode; odesoln = @simpleode_solution;
nsolvesvec = [50 100 200 500]; kernel = 'uniform';  % sqexp or uniform
%Nvec = [50 100 200 500 1000 2000 5000 10000];
Nvec = [25 50 100];
logCPUtime = zeros(num_iter, length(Nvec), length(nsolvesvec));
logGPUtime = zeros(num_iter, length(Nvec), length(nsolvesvec)); 
logGPUsparsetime = zeros(num_iter, length(Nvec), length(nsolvesvec));

%figure
state = 2;
for iterind = 1:num_iter
u0 = [unifrnd(-5, 5), unifrnd(-5, 5)]; u0GPU = gpuArray(u0);

for nind = 1:length(Nvec)
 %   subaxis(1,length(Nvec),nind) % you may also use subplot here
   N = Nvec(nind);
   ds = range(sspan)/(N-1);
  lambda = 1*ds; alpha = N/100;

for solveind = 1:length(nsolvesvec)
nsolves = nsolvesvec(nind);

%We're not really interested in the actual results here. Prelminary experiements
%showed that the GPU should be consistent with the CPU. Only time matters here.
    [~,~,logCPUtime(iterind, nind, solveind),~]  = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta);

[~,~,logGPUtime(iterind, nind, solveind),~] = uqdesGPU(sspanGPU, nsolves, N, kernel, lambda, alpha, odefn, u0GPU, thetaGPU);

[~,~, logGPUsparsetime(iterind, nind, solveind), ~] = uqdesGPUsparse(sspanGPU, nsolves, N, kernel, lambda, alpha, odefn, u0GPU, thetaGPU);

%The very first run of both uqdes and uqdesGPU is always slower.
%This can't be the initalization overhead, because we start measuring time inside the functions.
%What is the system doing differently on the first iteration?


end
end
end

%Plot results
figure
for solveind = 1:length(nsolvesvec)
subplot(1,length(nsolvesvec), solveind)
loglog(Nvec, median(logCPUtime(:,:,solveind),1), 'r--o', Nvec, median(logGPUtime(:,:,solveind),1), 'b--o', Nvec, median(logGPUsparsetime(:,:,solveind),1), 'g--o')
title(strcat("numsolves = ", num2str(nsolvesvec(solveind))));
xlabel('N')
ylabel('Computation time (seconds)')
legend('CPU', 'GPU', 'GPU w/sparsity', 'Location', 'southeast')
end

saveas(gcf, '~/uqdesShaun/All-other-forward-problems-and-figures/sparseGPUtest1.png')

%A bunch of commented-out stuff from Oksana's original Figure1 code
 %   [ueuler,teuler] = euler(sspan,N,odefn,u0,theta);
 %   truth = odesoln(t,theta);
 %   truth = truth(1:2,:)';
 %   eu = plot(teuler,ueuler(:,state),'g--');
 %   hold on
 %   for pp = 1:nsolves
 %       tempt = sspan(1):0.01:sspan(2);
 %       tempensemble = spline(t,uensemble(:,state,pp)',tempt);

 %       xflip       = [tempt fliplr(tempt)];
 %       yflip       = [tempensemble fliplr(tempensemble)];
 %       p           = patch(xflip,yflip,'r','EdgeAlpha',0.025,'FaceColor','none','linewidth',1.5);
%   end
  %  eu = plot(teuler,ueuler(:,state),'g--');
  %  tr = plot(t,truth(:,state),'r-');
  %  xlabel('t')
   % if nind == 1; ylabel('u'); end
  %  axis([sspan(1),sspan(2),-4,4])
  %  box off
    %legend([eu,p,tr],'euler','uqdes','true')

