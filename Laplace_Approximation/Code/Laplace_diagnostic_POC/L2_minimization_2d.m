v = 38;
d = 2;

tau = @(x) mvtpdf(x, eye(d), v);
tau_at_mode = tau(zeros([1 d]));

hess_at_mode = -diag(repelem((v+d)/v, d));
gauss_approx = @(x) tau_at_mode*exp(sum((x*hess_at_mode).*x,2)/2);
%gam = sqrt(1.5*(v+d)/(v+d-3));
gam = 3;

g = @(x) mvnpdf(x, zeros([1 d]), -gam^2*hess_at_mode);
mx_0 = @(x) gauss_approx(x)./g(x);

s_star = [[-3:3 zeros([1 6])]; [zeros([1 7]) -3:-1 1:3]]';

[xgrid, ygrid] = meshgrid(-10:0.01:10, -10:0.01:10);
x = [xgrid(:) ygrid(:)];
clear xgrid ygrid

s = s_star/sqrt(-hess_at_mode(1,1));

interrs = (tau(s)./g(s)) - mx_0(s);

mah_dist_s = pdist2(s_star, s_star);

mah_dist_x_s = pdist2(x*sqrt(-hess_at_mode(1,1)), s_star);

%These lines more or less borrowed from FSSKQ code (kq_kernel.m)
Cx_0 = @(r, l) exp(-r.^2/(2*l^2));
dCx_0 = @(r,l) r.^2 .* exp(-r.^2/(2*l.^2))./l.^3;

gvals = g(x);
prior_diffs =  gauss_approx(x) - tau(x);

MSE = @(l) sum((prior_diffs + gvals.*...
    Cx_0(mah_dist_x_s, l)*(Cx_0(mah_dist_s, l)\interrs)).^2);

dMSE = @(l) 2*sum((prior_diffs + gvals.*...
    Cx_0(mah_dist_x_s, l)*(Cx_0(mah_dist_s, l)\interrs)).*gvals.*...
    ((dCx_0(mah_dist_x_s, l) - (Cx_0(mah_dist_x_s, l)/Cx_0(mah_dist_s, l))*...
    dCx_0(mah_dist_s, l))*(Cx_0(mah_dist_s, l)\interrs)));

MSE_with_grad = @(y) deal(MSE(y),dMSE(y));

