skew_vec = 0.25:0.05:1.75;
eps_vec = 0.1:0.05:3;
b_vec = 10.^(-1:0.5:2); lambda = 2;

power_arr = zeros(length(skew_vec), length(eps_vec), length(b_vec));

syms t epsilon skew b real
assumeAlso(epsilon > 0)
%assume(delta > epsilon)
assumeAlso(skew<2); assumeAlso(b>0); assumeAlso(skew > 0);
%assumeAlso(lambda > 0);
%a = 2; b = 4;
a = 4/skew^2;
mle = (a-1)/b; mu = a/b; sd = sqrt(a)/b;

%Gamma pdf
f(t) = (b^a)*t^(a-1)*exp(-b*t)/gamma(a);
%Hyperparameters and endpoints for GP solver
%Make sure that the left onef doesn't go beyond zero
alpha = 1; %lambda = 1;
%delta = 2*epsilon;
s1 = piecewise(epsilon < sqrt(a)-1/sqrt(a),mle-epsilon*sd, 0); 
%s1 = 0;
%s1 = mle-epsilon*sd;
s3 = mle+epsilon*sd;
%Prior mean function:
m_t(t) = f(mle)*exp(0.5*subs(diff(log(f),2),mle)*(t-mle)^2);

%Interrogation points
s = [s1 mle s3]';

%Here we go
correction(epsilon, skew, b) = RQ1d_se(s3,s,lambda,s1,s3)*inv(RR1d_se(s,s,lambda,s1,s3))*RQ1d_se(s,s3,lambda,s1,s3);
variance(epsilon, skew, b) = simplify(QQ1d_se(s3,s3,lambda,s1,s3) - correction(epsilon, skew, b));
bias(epsilon, skew, b) = simplify((RQ1d_se(s3,s,lambda,s1,s3)*inv(RR1d_se(s,s ...
    ,lambda,s1,s3))*(f(s)-m_t(s)))^2);

KL(epsilon, skew, b) = 0.5*alpha*bias(epsilon,skew, b)/variance(epsilon,skew, b);

power = simplify(1+(erf((-1.96-bias/variance)/sqrt(2)) - erf((1.96-bias/variance)/sqrt(2)))/2);


%dKL_eps = diff(KL, epsilon);
%dKL_lam = diff(KL, lambda);

%fixeps = 3;
%B(epsilon) = subs(bias, lambda, fixeps);
%V(epsilon) = subs(variance, lambda, fixeps);
%K(epsilon) = subs(KL, lambda, fixeps);

