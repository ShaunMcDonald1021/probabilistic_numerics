%% PowerCalc.m
% Another script I made to play around with symbolic math for the
% experiments for my JSM 2018 presentation. Using the diagnostic in 1-D
% with 3 interrogation points, we want to optimize interrogation spacing
% w.r.t. KL divergence and/or test power. These two should be equivalent,
% and we can investigate that with this script.

% If you want to record values for power in an array at varying parameter
% values. Presumably you'd fill this in with a loop, but for some reason
% 2018 Shaun didn't put that in this script
skew_vec = 0.25:0.05:1.75;
eps_vec = 0.1:0.05:3;
b_vec = 10.^(-1:0.5:2); lambda = 2;
power_arr = zeros(length(skew_vec), length(eps_vec), length(b_vec));

syms t epsilon skew b real
% Define shape parameter in terms of skew for a different interpretation
% that may be somewhat more useful
a = 4/skew^2;
assumeAlso(epsilon > 0 & skew < 2 & b> 0 & skew > 0);

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

%Here we go
correction(epsilon, skew, b) = RQ1d_se(s3,s,lambda,s1,s3)*inv(RR1d_se(s,s,lambda,s1,s3))*RQ1d_se(s,s3,lambda,s1,s3);
variance(epsilon, skew, b) = simplify(QQ1d_se(s3,s3,lambda,s1,s3) - correction(epsilon, skew, b));
bias(epsilon, skew, b) = simplify((RQ1d_se(s3,s,lambda,s1,s3)*inv(RR1d_se(s,s ...
    ,lambda,s1,s3))*(f(s)-m_t(s)))^2);

KL(epsilon, skew, b) = 0.5*alpha*bias(epsilon,skew, b)/variance(epsilon,skew, b);

% Again, power and KL should be 1-1
power = simplify(1+(erf((-1.96-bias/variance)/sqrt(2)) - ...
    erf((1.96-bias/variance)/sqrt(2)))/2);
