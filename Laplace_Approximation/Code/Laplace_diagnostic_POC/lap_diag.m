function [post_mean, post_var]...
    = lap_diag(logf_inter, logf_at_mode, Us, s_star, logdet_T, d,...
    gam, alph, lap_app, wce, w)

Y = exp(d*log(gam)+logdet_T)*...
    (exp(logf_inter -log(mvnpdf(s_star/gam)))-...
    (2*pi)^d*exp(logf_at_mode + log(mvnpdf(sqrt(gam^2-1)*s_star/gam))));

% This part stolen from the kq_fss source code
% Compute the FSSKQ approximation
  Q = 0;
  J = length(w);
  Ls = zeros(J,1);
  ind = 0;
  for i = 1:J
    Ls(i) = size(Us{i}, 2);
    Q = Q + w(i) * sum(Y(ind+1:ind+Ls(i)));
    ind = ind + Ls(i);
  end
  
  post_mean = Q + lap_app;

post_var = wce^2*exp(2*(logf_at_mode + logdet_T) - d*log(alph));

lower_bound = min([norminv(0.025, post_mean, sqrt(post_var)) lap_app]);
upper_bound = max([norminv(0.975, post_mean, sqrt(post_var)) lap_app]);
lap_dist = fplot(@(x) normpdf(x, post_mean, sqrt(post_var)),...
    [lower_bound - 0.35*(upper_bound - lower_bound) upper_bound+...
    0.35*(upper_bound - lower_bound)], '-k','LineWidth', 1.5);
hold on
lap_line = xline(lap_app, '--b', 'LineWidth', 1.25);
quant_line = xline(norminv(0.025, post_mean, sqrt(post_var)), ':k',...
    'LineWidth', 1.35);
quant_line = xline(norminv(0.975, post_mean, sqrt(post_var)), ':k',...
    'LineWidth', 1.35);
legend([lap_dist, quant_line, lap_line],...
    {'Integral posterior', '95% interval',...
    'Laplace approximation'}, 'Location', 'southoutside');
hold off

