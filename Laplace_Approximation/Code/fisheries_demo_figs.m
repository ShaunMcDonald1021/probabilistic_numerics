load('1970_demo_stuff.mat')
samp_sizes = 1e6*2.^(-2:7);
figure(1)
fig = tiledlayout(1,2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
hSub1 = nexttile(1);
errorbar(samp_sizes, imp_means, 1.96*sqrt(imp_vars./samp_sizes),...
    ':.', 'Color', 'k', 'MarkerSize', 18, 'LineWidth', 1.5)
set(hSub1, 'XScale', 'log', 'XLim',...
    [min(samp_sizes) - 1e5 max(samp_sizes) + 4e7], 'FontSize', 18,...
    'TickLabelInterpreter', 'latex', 'XMinorTick', 'off', 'box', 'off')
xticks(samp_sizes)
xticklabels(samp_sizes)
xlabel('Sample size', 'Interpreter', 'latex', 'FontSize', 20)
title('IS likelihood estimates', ...
        'Interpreter', 'latex', 'FontSize', 22);

d = 72;
v = 25921;
gam = sqrt(1.5*(v+d)/(v+d-3));
lambda = 3.7;
lambda_GH = 2.06;
% lambda_GH = 2.12;

XS = gh_seq(2);
us = 3.6*sparse_gens(XS, d);
%us = zeros([72 3]);
%us(1:2,2) = 6;
%us(1, 3) = sqrt(d);
Us_GH = fss_gen(us);
Us_GH = Us_GH(2:end);
s_star_GH = cell2mat(Us_GH)';

load('1970_diag.mat')
Us = fss_gen(s_star(:,[1 73]));
s_star = s_star';
logf_interrs = logf_interrs';
[~, w, wce, alph] = diag_calib(gam, lambda, d, v, Us, s_star);
nexttile(2);
[post_mean, post_var, ~] = lap_diag(logf_interrs, logf_at_mode,...
    log_T_det, d, gam, alph, Us, w, wce, s_star);

nexttile(1);
lap_line = yline(exp(log_T_det + logf_at_mode + d*log(2*pi)/2),...
    'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.75);
post_line = yline(post_mean, 'Color', 'r', 'LineStyle', '--',...
    'LineWidth', 1.75);
legend([lap_line, post_line],...
        {'Laplace approximation', '$m_1$',},...
        'interpreter', 'latex', 'FontSize', 18);
    
figure(2)
fig = tiledlayout(2,1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
hSub1 = nexttile(1);
histogram(small_pvals, 12, 'FaceColor', [0.75 0.75 0.75]');
set(hSub1, 'FontSize', 18,...
    'TickLabelInterpreter', 'latex', 'xtick', [], 'box', 'off')
title('$n = 100$', ...
        'Interpreter', 'latex', 'FontSize', 22)
pval =  2*(1-normcdf(abs(post_mean - exp(log_T_det + logf_at_mode +...
    d*log(2*pi)/2))/sqrt(post_var)));
pval_line = xline(pval, 'Color', 'r', 'LineStyle', '--',...
    'LineWidth', 1.75);
legend([pval_line],...
        {'Diagnostic p-value',},...
        'interpreter', 'latex', 'FontSize', 18);
hSub2 = nexttile(2);
histogram(large_pvals, 12, 'FaceColor', [0.75 0.75 0.75]');
set(hSub2, 'FontSize', 18,...
    'TickLabelInterpreter', 'latex', 'box', 'off')
title('$n = 1000$', ...
        'Interpreter', 'latex', 'FontSize', 22)
xline(pval, 'Color', 'r', 'LineStyle', '--',...
    'LineWidth', 1.75)

load('2005_demo_stuff.mat')
figure(3)
fig = tiledlayout(1,2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
hSub1 = nexttile(1);
errorbar(samp_sizes, imp_means, 1.96*sqrt(imp_vars./samp_sizes),...
    ':.', 'Color', 'k', 'MarkerSize', 18, 'LineWidth', 1.5)
set(hSub1, 'XScale', 'log', 'XLim',...
    [min(samp_sizes) - 1e5 max(samp_sizes) + 4e7], 'FontSize', 18,...
    'TickLabelInterpreter', 'latex', 'XMinorTick', 'off', 'box', 'off')
xticks(samp_sizes)
xticklabels(samp_sizes)
xlabel('Sample size', 'Interpreter', 'latex', 'FontSize', 20)
title('IS likelihood estimates', ...
        'Interpreter', 'latex', 'FontSize', 22);

load('2005_diag.mat')
s_star = s_star';
logf_interrs = logf_interrs';
nexttile(2);
[post_mean, post_var, ~] = lap_diag(logf_interrs, logf_at_mode,...
    log_T_det, d, gam, alph, Us, w, wce, s_star);

nexttile(1);
lap_line = yline(exp(log_T_det + logf_at_mode + d*log(2*pi)/2),...
    'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.75);
post_line = yline(post_mean, 'Color', 'r', 'LineStyle', '--',...
    'LineWidth', 1.75);
legend([lap_line, post_line],...
        {'Laplace approximation', '$m_1$',},...
        'interpreter', 'latex', 'FontSize', 18);
    
figure(4)
fig = tiledlayout(2,1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
hSub1 = nexttile(1);
histogram(small_pvals, 12, 'FaceColor', [0.75 0.75 0.75]');
set(hSub1, 'FontSize', 18,...
    'TickLabelInterpreter', 'latex', 'xtick', [], 'box', 'off')
title('$n = 100$', ...
        'Interpreter', 'latex', 'FontSize', 22)
pval =  2*(1-normcdf(abs(post_mean - exp(log_T_det + logf_at_mode +...
    d*log(2*pi)/2))/sqrt(post_var)));
pval_line = xline(pval, 'Color', 'r', 'LineStyle', '--',...
    'LineWidth', 1.75);
legend([pval_line],...
        {'Diagnostic p-value',},...
        'interpreter', 'latex', 'FontSize', 18);
hSub2 = nexttile(2);
histogram(large_pvals, 12, 'FaceColor', [0.75 0.75 0.75]');
set(hSub2, 'FontSize', 18,...
    'TickLabelInterpreter', 'latex', 'box', 'off')
title('$n = 1000$', ...
        'Interpreter', 'latex', 'FontSize', 22)
xline(pval, 'Color', 'r', 'LineStyle', '--',...
    'LineWidth', 1.75)

    
load('1970_diag_GH.mat')
logf_interrs = logf_interrs';
figure(5)
fig = tiledlayout(1,2, 'TileSpacing', 'Loose', 'Padding', 'Compact');
hSub1 = nexttile(1);
[~, w, wce, alph] = diag_calib(gam, lambda_GH, d, v, Us_GH, s_star_GH);
[post_mean, post_var, point_contribs] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, Us_GH, w, wce, s_star_GH);
hSub2 = nexttile(2);
radii = unique(sqrt(sum(s_star_GH.^2, 2)));
plot(radii, point_contribs, '.', 'Color', 'k', 'MarkerSize', 18)
set(hSub2, 'FontSize', 18,...
    'TickLabelInterpreter', 'latex', 'box', 'off')
  %  'Ylim', [-2e-4 1e-5])
for i = 1:length(point_contribs)
    line([radii(i) radii(i)], [0 point_contribs(i)], 'LineWidth', 1.75,...
        'Color', 'k')
end
xlabel('$|\!|s^*|\!|$', 'Interpreter',...
    'latex', 'FontSize', 20)
ylabel('Tot. contribution to $m_1 - L\left(p_{xy}\right)$', 'Interpreter',...
    'latex', 'FontSize', 20)
title('Contributions to integral', ...
        'Interpreter', 'latex', 'FontSize', 22)
    
load('2005_diag_GH.mat')
logf_interrs = logf_interrs';
figure(6)
fig = tiledlayout(1,2, 'TileSpacing', 'Loose', 'Padding', 'Compact');
hSub1 = nexttile(1);
[post_mean, post_var] = lap_diag(logf_interrs, logf_at_mode,...
        log_T_det, d, gam, alph, Us_GH, w, wce, s_star_GH);
hSub2 = nexttile(2);
[post_mean, post_var] = lap_diag(logf_interrs([1:2997 2999:end]),...
    logf_at_mode, log_T_det, d, gam, alph, Us_GH, lambda_GH, [],...
    s_star_GH([1:2997 2999:end],:), 'is_w', false, 'is_symm', false);
    
    
    