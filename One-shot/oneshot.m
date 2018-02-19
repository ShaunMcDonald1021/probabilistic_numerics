function [uensemble,t,logfntime,m_deriv_svec] ...
     = oneshot(tspan, nsolves, N, kernel, lambda, alpha, fdefn, theta, s)


%% MATLAB implementation of the one-shot probabilistic updater to assess the quality of the Laplace approximation.
%% Building on Charlie's thesis, and (sloppily and shamelessly) retooled from Oksana's code.
%% I'll clean up the comments and whatnot once I have a better handle on all this.
%% Just 1-dimensional for now; we'll figure out higher dimensions later.

%% uqdes.m 
%
% DESCRIPTION
% Bayesian updating algorithm by Oksana A. Chkrebtii last updated
% May 20, 2016. The probabilistic solver is described int he paper Bayesian
% Uncertainty Quantification for Differential Equations, by O.A. Chkrebtii,
% D.A. Campbell, B. Calderhead, M.A. Girolami. This implementation allows
% the option to evaluate the forward model outside of the discretization
% grid.
%
% INPUTS
% tspan = 1x2 vector containing the upper and lower bound of the evaluation domain
% nsolves = number of posterior samples to take
% N = an integer for the size of the evaluation grid
% fdefn = a function handle for our likelihood (and will also give the MLE)
% theta = various parameter(s) for the underlying function. Perhaps literally theta
% as meant in the paper? Again, we will come back to this
% kernel = covariance kernel, either squared exponential or uniform.
% Once we generalize to higher dimensions, we will need to consider this more carefully.
% alpha/lambda = kernel hyperparameters.
% s = 1xq optional grid of interrogation points. If not supplied, it will be set
% equal to t.
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
B = nsolves; 
t  = linspace(tspan(1),tspan(2),N); ds = t(2)-t(1);

switch lower(kernel) % sort out which kernel to use
    case {'sqexp'};  kern = 'se'; trim = N;
    case {'uniform'}; kern = 'un'; trim = 2*ceil(lambda/ds);
    otherwise; disp('Unknown kernel -- try agaian'); return;
end

% define appropriate kernel convolutions
QQ = str2func(strcat('QQ1d_',kern));
RR = str2func(strcat('RR1d_',kern));
QR = str2func(strcat('QR1d_',kern));

%uensemble = repmat(u0,[length(t),1,B]);
[f, mle] = fdefn(theta);
laplace_mean = mle;
syms t
laplace_covar = -1/vpa(subs(diff(log(f(t)),t,2),mle));

m_deriv = normpdf(t, laplace_mean, laplace_covar);
m_state = 

m_deriv_svec = repmat(0*u0,[length(t),1,B]);
m_state_svec =  uensemble + bsxfun(@times,m_deriv_svec,repmat(t',1,M));
C_deriv_ssmat = RR(t,t,lambda,sspan(1),sspan(2))/alpha;
C_state_ssmat = QQ(t,t,lambda,sspan(1),sspan(2))/alpha;
C_cross1_ssmat = QR(t,t,lambda,sspan(1),sspan(2))/alpha;
kinv = 1/(C_deriv_ssmat(1,1));
f_diff = kinv*(f-m_deriv_svec(1,:,:));      
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
    if sum(f_diff >= 1e10 | isnan(f_diff) | isinf(f_diff))>0
        disp('Algorithm failed to converge: try increasing mesh size or changing assumptions')
        uensemble = [];
        logfntime = [];
        m_deriv_svec = [];
        return
    end
end
    
if nargin == 10
    randnNums  = randn(length(sinds),M,B);  % generate random numbers outside of loop
    uensemble(sinds,:,:) = m_state_svec(sinds,:,:) + bsxfun(@times,randnNums,sqrt(diag(C_state_ssmat(sinds,sinds))));   
end

logfntime = log(toc); % end timer


