syms t epsilon b a real
assumeAlso(epsilon > 0)
%assume(delta > epsilon)
assumeAlso(a>1); assumeAlso(b>0);
%assumeAlso(lambda > 0);
%b = 3;
mle = (a-1)/b; mu = a/b; sd = sqrt(a)/b;

lambda = 1.5*epsilon*sd;

%Gamma pdf
f(t) = piecewise(t>0,(b^a)*t^(a-1)*exp(-b*t)/gamma(a),0);
%Hyperparameters and endpoints for GP solver
%Make sure that the left onef doesn't go beyond zero
alpha = 1; %lambda = 1;
%delta = 2*epsilon;
%s1 = piecewise(epsilon<sqrt(a)-1/sqrt(a),mle-epsilon*sd,0);
%s1 = 0;
s1 = mle-epsilon*sd;
s3 = mle+epsilon*sd;
%Prior mean function:
m_t(t) = f(mle)*exp(0.5*subs(diff(log(f),2),mle)*(t-mle)^2);

%Interrogation points
s = [s1 mle s3]';

%Here we go
correction(epsilon, a) = RQ1d_se(s3,s,lambda,s1,s3)*inv(RR1d_se(s,s,lambda,s1,s3))*RQ1d_se(s,s3,lambda,s1,s3);
variance(epsilon, a) = simplify(QQ1d_se(s3,s3,lambda,s1,s3) - correction(epsilon, a));
bias(epsilon, a) = simplify(RQ1d_se(s3,s,lambda,s1,s3)*inv(RR1d_se(s,s,lambda,s1,s3))*(f(s)-m_t(s)));

KL(epsilon, a) = 0.5*alpha*bias^2/variance


