%% Comparison of GPU and CPU using the simple ODE from Figure 1 of Oksana's paper
%We'll run 30 iterations of both the CPU and GPU solvers at various N-values with randomly-generated initial conditions (generated earlier for the sake of consistency)
%This file runs it on the CPU. We'll calculate the times, save them, and plot them against the GPU times later.

clear all; close all; clc;

%Load initial condition matrix
load('u0mat.mat')

num_iter = 30; %number of iterations
sspan = [0 10]; theta = 2; 
odefn = @simpleode; odesoln = @simpleode_solution;
%nsolvesvec = [50 100 200 500]; This one didn't seem to matter too much anyway,
%and we're going to VERY HIGH N so I'd rather not do more repetitions than I have to
nsolvesvec = [100]; kernel = 'sqexp';  % sqexp or uniform
Nvec = [50 100 200 500 1000 2000 5000 10000];
%Nvec = [25];
%logCPUtime = zeros(num_iter, length(Nvec), length(nsolvesvec));
%logGPUtime = zeros(num_iter, length(Nvec), length(nsolvesvec)); 
load('logCPUtime.mat')
%figure
state = 2;

for nind = length(Nvec)
 %   subaxis(1,length(Nvec),nind) % you may also use subplot here
   N = Nvec(nind);
   ds = range(sspan)/(N-1);
  lambda = 1*ds; alpha = N/100;


%Not as efficient/logical with the loops ordered this way, because it means redundancy
%from repeatedly redefining u0. However, it means we get all of the small-N iterations done
%before we start on the really big N's, which will be useful if shit goes sideways there

for iterind = 1:num_iter
u0 = u0mat(iterind,:);

for solveind = 1:length(nsolvesvec)
nsolves = nsolvesvec(solveind);

tic %We're not really interested in the actual results here. Prelminary experiements
%showed that the GPU should be consistent with the CPU. Only time matters here.
[~,~,logCPUtime(iterind, nind, solveind),~]  = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta);

%[~,~,logGPUtime(iterind, nind, solveind),~] = uqdesGPU(sspanGPU, nsolves, N, kernel, lambda, alpha, odefn, u0GPU, thetaGPU);

disp(strcat("N = ", num2str(N), ", nsolves = ", num2str(nsolves), ", iteration ", num2str(iterind), " took ", num2str(toc), " seconds."))


%The very first run of both uqdes and uqdesGPU is always slower.
%This can't be the initalization overhead, because we start measuring time inside the functions.
%What is the system doing differently on the first iteration?

%Save on each iteration to be safe
save('logCPUtime.mat', 'logCPUtime')

end %End of solveind loop
end %End of iterind loop
end %End of nind loop

%Plot results
%figure
%for solveind = 1:length(nsolvesvec)
%subplot(1,length(nsolvesvec), solveind)
%loglog(Nvec, median(logCPUtime(:,:,solveind),1), 'r--o', Nvec, median(logGPUtime(:,:,solveind),1), 'b--o')
%title(strcat("numsolves = ", num2str(nsolvesvec(solveind))));
%xlabel('N')
%ylabel('Computation time (seconds)')
%legend('CPU', 'GPU', 'Location', 'southeast')
%end

%saveas(gcf, '~/uqdesShaun/All-other-forward-problems-and-figures/GPUtest1.png')

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
quit
