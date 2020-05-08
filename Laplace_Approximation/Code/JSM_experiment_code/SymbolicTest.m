%% SymbolicTest.m
% Just something I used in the experiments for my JSM 2018 poster. Using
% symbolic math to try and optimize the diagnostic in 1-D using 3
% interrogation points. The assumption is you'd want high KL divergence
% from the "null" (i.e. low variance, but bring it as far away as possible
% from the prior mean). For Gaussians of equal variance, KL is simply
% bias^2/(2*variance)

syms t epsilon b a real
% t = argument; a, b = gamma pdf parameters; epsilon = spacing of
% interrogation points in units of standard deviation
assumeAlso(epsilon > 0 & a > 1 & b > 0);
mle = (a-1)/b; mu = a/b; sd = sqrt(a)/b;

% Taking covariance hyperparameters to be functions of epsilon as in
% Charlie's thesis
lambda = 1.5*epsilon*sd;
alpha = 1/(epsilon*sd);

% Gamma pdf full function
f(t) = piecewise(t>0,(b^a)*t^(a-1)*exp(-b*t)/gamma(a),0);

% Hyperparameters and endpoints for GP solver
s1 = mle-epsilon*sd;
% Can force s1 not to go lower than 0 by uncommenting this line 
% s1 = piecewise(epsilon<sqrt(a)-1/sqrt(a),mle-epsilon*sd,0);
s3 = mle+epsilon*sd;

%Prior mean function
m_t(t) = f(mle)*exp(0.5*subs(diff(log(f),2),mle)*(t-mle)^2);

%Interrogation points
s = [s1 mle s3]';

% Here we go
correction(epsilon, a) = RQ1d_se(s3,s,lambda,s1,s3)*...
    inv(RR1d_se(s,s,lambda,s1,s3))*RQ1d_se(s,s3,lambda,s1,s3);
variance(epsilon, a) = simplify(QQ1d_se(s3,s3,lambda,s1,s3) -...
    correction(epsilon, a));
bias(epsilon, a) = simplify(RQ1d_se(s3,s,lambda,s1,s3)*...
    inv(RR1d_se(s,s,lambda,s1,s3))*(f(s)-m_t(s)));

% KL is ultimately a function of epsilon and gamma shape. The scale
% parameter more or less controls scaling, and therefore has no effect on
% the optimal epsilon
KL(epsilon, a) = 0.5*alpha*bias^2/variance
