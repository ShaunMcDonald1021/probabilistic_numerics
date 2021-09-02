function [post_mean, w, wce, alph] = diag_calib(gam, lambda, s_star, Us, d, v)


logf_inter = log(mvtpdf(s_star.*sqrt(v/(v+d)), eye(d), v));
logf_at_mode = (gammaln((v+d)/2) - gammaln(v/2) - d*log(v*pi)/2);
logdet_T = d*log(v/(v+d))/2;
lap_app = exp(logdet_T + logf_at_mode + d*log(2*pi)/2);
%disp(lap_app)

Y = exp(d*log(gam)+logdet_T)*...
    (exp(logf_inter -log(mvnpdf(s_star/gam)))-...
    (2*pi)^d*exp(logf_at_mode + log(mvnpdf(sqrt(gam^2-1)*s_star/gam))));

k = @(r)exp(-r.^2/(2*lambda^2));
kmean = @(x)(lambda^2/(gam^2+lambda^2))^(d/2)*...
    exp(-norm(x)^2/(2*(gam^2+lambda^2)));
Ikmean = (lambda^2/(2*gam^2 + lambda^2))^(d/2);

[Q, wce, w] = kq_fss(Y, Us, k, kmean, Ikmean);

post_mean = Q + lap_app;

alph = (Q/(1.96*wce*exp(logdet_T + logf_at_mode)))^(-2/d);

