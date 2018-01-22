%% Comparison of GPU and CPU using the simple ODE from Figure 1 of Oksana's paper
%%We'll run 30 iterations of both the CPU and GPU solvers at various N-values with randomly-generated initial conditions (generated earlier for the sake of consistency/reproducibility)
%%This file runs 15 iterations on each GPU in parallel to save time.
%We'll record the times, save them, and plot them (with the CPU times) later

clear all; close all; clc;

%Load matrix of initial conditions
load('u0mat.mat')

num_iter = 30; %number of iterations (although I suppose we don't really need to make this a variable with the current setup)
sspan = [0 10];  theta = 2; 
odefn = @simpleode; odesoln = @simpleode_solution;
%nsolvesvec = [50 100 200 500]; This one didn't seem to matter too much anyway,
%and we're going to VERY HIGH N so I'd rather not do more repetitions than I have to
nsolvesvec = [100]; kernel = 'sqexp';  % sqexp or uniform
Nvec = [50 100 200 500 1000 2000 5000 10000];

%Nvec = [25];
%logCPUtime = zeros(num_iter, length(Nvec), length(nsolvesvec));
logGPUtimepar = zeros(15, length(Nvec), length(nsolvesvec)); 
%15 iterations on each GPU, then we'll concatenate after

parpool('local', 2) %Set up parpool
spmd
gpuDevice(labindex);
thetaGPU = gpuArray(theta);
sspanGPU = gpuArray(sspan);
u0matGPU = gpuArray(u0mat);
end %Allocate GPU's and gpuArray's now so we don't have to do it in loops later

state = 2;

for nind = length(Nvec):-1:1 %I think running it backwards will avoid memory fragmentation

   N = Nvec(nind); %For some reason, the progress message doesn't display unless I show N as well
%(presumably any numeric variable would work here, but it's weird that it suppresses output unless you let a numeric variable display)
   ds = range(sspan)/(N-1);
  lambda = 1*ds; alpha = N/100;

%This loop structure is less efficient, as it means we are repeatedly going in and out of 
%spmd structures and redefining the initial conditions more times than necessary.
%However, this way allows us to concatenate, gather, and save the times after every N
%(so unless something goes wrong within a really big N, we can recover some of the results if we crash)

spmd
for iterind = 1:15
u0GPU = u0matGPU(15*(labindex-1) + iterind,:);

for solveind = 1:length(nsolvesvec)
nsolves = nsolvesvec(solveind) %Guess you need the number to be displayed in the same loop in order for the "disp" to work?
tic
%We're not really interested in the actual results here. Prelminary experiements
%showed that the GPU should be consistent with the CPU. Only time matters here.
%[~,~,logCPUtime(iterind, nind, solveind),~]  = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta);

[~,~,logGPUtimepar(iterind, nind, solveind),~] = uqdesGPU(sspanGPU, nsolves, N, kernel, lambda, alpha, odefn, u0GPU, thetaGPU);

disp(strcat("N = ", num2str(N), ", nsolves = ", num2str(nsolves), ", iteration ", num2str(15*(labindex-1)+iterind), " took ", num2str(toc), " seconds."))

%The very first run of both uqdes and uqdesGPU is always slower.
%This can't be the initalization overhead, because we start measuring time inside the functions.
%What is the system doing differently on the first iteration?
%Whatever it is, this behaviour manifests itself on the GPU AND the CPU.

end %end solveind loop
end %end iterind loop
logGPUtimegat = gcat(logGPUtimepar, 1);
end %end spmd

logGPUtime = logGPUtimegat{1};
save('logGPUtime.mat', 'logGPUtime');

end %end nind loop

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
%quit
