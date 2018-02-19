function [uensemble, t, logfntime, m_deriv_svec] ...
    = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta,tgrid)

%% uqdesGPU.m 
%
% DESCRIPTION
% Modification of the uqdes algorithm to allow for the use of gpuArray's. The original comments are below, along with some of my own. -Shaun

% Bayesian updating algorithm by Oksana A. Chkrebtii last updated
% May 20, 2016. The probabilistic solver is described int he paper Bayesian
% Uncertainty Quantification for Differential Equations, by O.A. Chkrebtii,
% D.A. Campbell, B. Calderhead, M.A. Girolami. This implementation allows
% the option to evaluate the forward model outside of the discretization
% grid.
%
% INPUTS
% sspan = 1x2 vector containing the upper and lower bound of the domain of
% integration
% N = an integer for the discretizatio mesh size
% odefn = a function handle for the system of ODEs
% u0 = 1xM vector of initial conditions
% theta = parameters used in the system of ODEs, dimension of theta depends
% on the inputs of the function handle odefn
% tgrid = 1xq optional parameter passed if we wish to evaluate the forward
% model at time points lying outside the discretization grid
% 
% OUTPUTS
% This function returns 
% t: the discretization grid at which the function was evaluated (for the
% moment these are equally spaced)
% logftime: gives the log computation time in seconds
% 
% EXAMPLE
% [umat,tvec,lft] = uqdes([0,10],5,50,'sqexp',1,50,@simpleode,[-1,0],2);


tic % start timing  

% preallocate variables
M = length(u0); B = nsolves; 
s = linspace(sspan(1),sspan(2),N); ds = s(2)-s(1);
if nargin == 10 
    t = sort(unique([s,tgrid]));
    [~,tinds] = intersect(t,s);
    [~,sinds] = intersect(t,tgrid);
else t = s; tinds = (1:N)'; sinds = []; end; 
if size(u0,1)>size(u0,2); u0 = u0'; end;

switch lower(kernel) % sort out which kernel to use
    case {'sqexp'};  kern = 'se'; trim = N;
    case {'uniform'}; kern = 'un'; trim = 2*ceil(lambda/ds);
    otherwise; disp('Unknown kernel -- try agaian'); return;
end

% define appropriate kernel convolutions
QQ = str2func(strcat('QQ1d_',kern));
RR = str2func(strcat('RR1d_',kern));
QR = str2func(strcat('QR1d_',kern));

uensemble = repmat(u0,[length(t),1,B]);
f = odefn(s(1),uensemble(1,:,:),theta);
% there are two possibilities for the prior mean: zero derivative, or
% constant exact derivative f(s(1),u0): either one works
m_deriv_svec = repmat(0*u0,[length(t),1,B]);
%m_deriv_svec = repmat(f,[length(t),1,1]);
m_state_svec =  uensemble + bsxfun(@times,m_deriv_svec,repmat(t',1,M));
C_deriv_ssmat = RR(t,t,lambda,sspan(1),sspan(2))/alpha;
C_state_ssmat = QQ(t,t,lambda,sspan(1),sspan(2))/alpha;
C_cross1_ssmat = QR(t,t,lambda,sspan(1),sspan(2))/alpha;
kinv = 1/(C_deriv_ssmat(1,1));
f_diff = kinv*(f-m_deriv_svec(1,:,:));      

%Sadly, MATLAB doesn't support the use of gpuArray's in if statements. This makes convergence checking at each time iteration problematic. The options are:
%1) pull f_diff off the GPU at every iteration of the loop. The inefficiency of this would defeat the whole purpose of using GPU's in the first place.
%2) Only gather and check convergence at the last time step (after the loop). This the most efficient option. It should be fine if the convergence errors are NaN's or Inf's, as these would presumably propagate through iterations. The problem is the f_diff >= 1e10 condition: theoretically, nothing is stopping this from happening in an earlier iteration but NOT happening in the final one. In this case, just looking at the last iteration may lead us to falsely conclude convergence. I'm not sure how likely this would be, so we'll go with this option for now.
%3) Turn f_diff into a 4-D array (where the fourth dimension is N), then gather it at the end of the loop and check convergence at each time step in another loop. This is arguably a good balance between safety and efficiency, but if N is large it may still be too inefficient. It would also entail quite a bit more code modification. If we decide to come back to this, uncomment the line below.

%f_diff = repmat(kinv*(f-m_deriv_svec(1,:,:)),1,1,1,N);
randnNums  = randn(length(t),M,B);  % generate random numbers outside of loop
counter = 0;

% run one-step uqdes algorithm
for n = tinds(1:end-1)'
    ind = tinds(max(counter+1,max(1,counter-trim)):min(end,min(end,counter+trim)));
	endind = vertcat(tinds(counter+1:N),sinds);
    counter = counter + 1;
    nextind = tinds(counter+1);
    m_state_svec(endind,:,:) = m_state_svec(endind,:,:) + bsxfun(@times,C_cross1_ssmat(endind,n),f_diff(1,:,:));
    m_deriv_svec(ind,:,:) = m_deriv_svec(ind,:,:) + bsxfun(@times,C_deriv_ssmat(ind,n),f_diff(1,:,:));
    C_state_ssmat(endind,endind) = C_state_ssmat(endind,endind) - kinv*C_cross1_ssmat(endind,n)*(C_cross1_ssmat(endind,n))';
    C_cross1_ssmat(endind,ind) = C_cross1_ssmat(endind,ind) - kinv*C_cross1_ssmat(endind,n)*C_deriv_ssmat(n,ind);
    C_deriv_ssmat(ind,ind) = C_deriv_ssmat(ind,ind) - kinv*C_deriv_ssmat(ind,n)*C_deriv_ssmat(n,ind);
    uensemble(nextind,:,:) = m_state_svec(nextind,:,:) + sqrt(C_state_ssmat(nextind,nextind))*randnNums(n,:,:);   
    kinv = 1/(C_deriv_ssmat(nextind,nextind)+C_deriv_ssmat(tinds(counter),tinds(counter)));
    f_diff = kinv*(odefn(t(nextind),uensemble(nextind,:,:),theta) - m_deriv_svec(nextind,:,:));
%    if sum(f_diff >= 1e10 | isnan(f_diff) | isinf(f_diff))>0
 %       disp('Algorithm failed to converge: try increasing mesh size or changing assumptions')
  %      uensemble = [];
   %     logfntime = [];
    %    m_deriv_svec = [];
    %    return
   % end
end

%Gather f_diff to check convergence at the last step
f_diff_gather = gather(f_diff);

%Should there be a condition for f_diff <= -1e10 as well? Large magnitudes in either direction would indicate failure to converge, would they not?
if sum(f_diff_gather >= 1e10 | isnan(f_diff_gather) | isinf(f_diff_gather))>0
 disp('Algorithm failed to converge: try increasing mesh size or changing assumptions')
        uensemble = [];
        logfntime = [];
        m_deriv_svec = [];
        return
end

    
if nargin == 10
    randnNums  = randn(length(sinds),M,B);  % generate random numbers outside of loop
    uensemble(sinds,:,:) = m_state_svec(sinds,:,:) + bsxfun(@times,randnNums,sqrt(diag(C_state_ssmat(sinds,sinds))));   
end

%Gather everything we need off the GPU
%A more effieicnet approach may be to have this function output gpuArray's, then only gather the stuff we need in each application. However, doing it this way allows us to assess the "worst-case" speed of GPU usage
uensemble = gather(uensemble);
t = gather(t);
m_deriv_svec = gather(m_deriv_svec);

wait(gpuDevice)
logfntime = log(toc); % end timer


