syms t epsilon real
%Gamma parameters
a = 2; b = 4;
mle = (a-1)/b; mu = a/b; sd = sqrt(a)/b;

%Gamma pdf
f(t) = (b^a)*t^(a-1)*exp(-b*t)/gamma(a);

%Hyperparameters and endpoints for GP solver
%Make sure that the left one doesn't go beyond zero
alpha = 1; lambda = 1;
%delta = 2*epsilon;
s1 = max(mle-2*epsilon*sd, 0); s2 = max(mle-epsilon*sd,0);
s4 = mle+epsilon*sd; s5 = mle+2*epsilon*sd;

%Prior mean function:
m_t(t) = f(mle)*exp(0.5*subs(diff(log(f),2),mle)*(t-mle)^2);

%Interrogation points
s = [s1 s2 mle s4 s5]';

%Here we go
correction(epsilon) = RQ1d_se(s5,s,lambda,s1,s5)*inv(RR1d_se(s,s,lambda,s1,s5))*RQ1d_se(s,s5,lambda,s1,s5);
variance(epsilon) = simplify(QQ1d_se(s5,s5,lambda,s1,s5) - correction(epsilon));

bias(epsilon) = simplify((RQ1d_se(s5,s,lambda,s1,s5)*inv(RR1d_se(s,s ...
    ,lambda,s1,s5))*(f(s)-m_t(s))))^2;

KL(epsilon) = 0.5*alpha*bias(epsilon)/variance(epsilon);