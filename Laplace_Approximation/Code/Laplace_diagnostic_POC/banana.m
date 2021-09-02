d = 2;
lambda = 2.9224;
gam = 1.2734;
alph = 0.0793;

f = @(x) mvnpdf([x(:,1) x(:,2) - (x(:,1).^2 - 3)/2], [], diag([3 1]));
mode = [0 -1.5];
hess_at_mode = diag([-1/3 -1]);
f_at_mode = f(mode);
gauss_approx = @(x) f_at_mode*exp(sum(((x-mode)*hess_at_mode).*(x-mode),2)/2);
g = @(x) mvnpdf(x, mode, -gam^2*inv(hess_at_mode));
s_star = [[-3:3 zeros([1 6])]; [zeros([1 7]) -3:-1 1:3]]';
mah_dist_s = pdist2(s_star, s_star);
T = diag([sqrt(3) 1]);
s = s_star*T + mode;
f_interrs = f(s);
gauss_interrs = gauss_approx(s);
mx_0 = @(x) gauss_approx(x)./g(x);

hold off
figure(1)
[xgrid_star, ygrid_star] = meshgrid(-4.5:0.05:4.5, -5.5:0.05:10);
xgrid = T(1,1)*xgrid_star + mode(1);
ygrid = T(2,2)*ygrid_star + mode(2);
x = [xgrid(:) ygrid(:)];
mah_dist_x_s = pdist2([xgrid_star(:) ygrid_star(:)], s_star);

grid_bounds = [min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))];
grid_ratio = (grid_bounds(4) - grid_bounds(3))/...
            (grid_bounds(2) - grid_bounds(1));
pbaspect([1 grid_ratio geomean([grid_ratio 1/grid_ratio])])
%axis vis3d
hold on

true_surf = fsurf(@(x, y) f([x y]), grid_bounds, 'FaceAlpha', 0.65);
colormap winter
set(gca, 'FontSize', 15)
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 20);
% This function handle returns NaN if gauss_approx < f_at_mode*10^-4, which
% makes the plot look cleaner
gauss_plot = @(x,y) gauss_approx([x y]).*...
    (gauss_approx([x y]) > f_at_mode*1e-4)./...
    (gauss_approx([x y]) > f_at_mode*1e-4);
gauss_surf = fsurf(gauss_plot, grid_bounds,...
    'FaceAlpha', 0.4, 'FaceColor', 'r', 'EdgeColor', 'y');
legend([true_surf, gauss_surf], {'$\beta$', 'Gaussian approximation'},...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 18);
hold off

figure(2)
fig = tiledlayout(1,3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
hSub1 = nexttile([1 2]);
set(gca, 'FontSize', 18)
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 20);
pbaspect([1 grid_ratio geomean([grid_ratio 1/grid_ratio])])
hold on

Cx_0 = @(r) exp(-r.^2/(2*lambda^2));
gram_mat = Cx_0(mah_dist_s);
pointz = gram_mat\((f_interrs - gauss_interrs)./g(s));
% mx_1 = @(x,y) gauss_approx([x y]) + ...
%     g([x y]).*(Cx_0(pdist2(([x y] - mode)*inv(T), s_star))*pointz);
% mx_1_plot= @(x,y) mx_1(x,y).*...
%     (abs(mx_1(x,y)) > f_at_mode*1e-4)./...
%     (abs(mx_1(x,y)) > f_at_mode*1e-4);
post_points = reshape(gauss_approx(x) + ...
    g(x).*(Cx_0(mah_dist_x_s)*pointz), size(xgrid,1), size(xgrid, 2));
% post_plot_points = post_points.*...
%     (abs(post_points) > f_at_mode*1e-4)./...
%     (abs(post_points) > f_at_mode*1e-4);
post_surf = surf(xgrid, ygrid, post_points ,...
    'FaceAlpha', 0.65, 'EdgeColor', 'none');
colormap winter
for i = 1:10:size(xgrid, 1)
    plot3(xgrid(i,:), ygrid(i,:), post_points(i,:), '-k',...
        'Color', [0.75 0.75 0.75]);
end
for i = 1:5:size(xgrid, 2)
    plot3(xgrid(:,i), ygrid(:,i), post_points(:,i), '-k',...
        'Color', [0.75 0.75 0.75]);
end

f_plot = @(x,y) f([x y]).*...
(f([x y]) > f_at_mode*1e-4)./...
(f([x y]) > f_at_mode*1e-4);
true_surf2 = fsurf(f_plot, grid_bounds,...
    'FaceAlpha', 0.35, 'FaceColor', 'red', 'EdgeColor', 'y');

interrs = scatter3(s(:,1), s(:,2), f_interrs, 64,...
    'filled', 'MarkerFaceColor', 'k');
title('True function and un-normalized GP posterior mean',...
    'Interpreter',"latex", 'FontSize', 22);
legend([post_surf, true_surf2, interrs], {'$m^x_1 \cdot g$', '$\beta$',...
    'Interrogation points'},'Location', 'north',...
    'Interpreter', 'latex', 'FontSize', 18);
hold off

hSub2 = nexttile;
int_weightz = (lambda^2/(gam^2+lambda^2))*...
    exp(-sum(s_star.^2, 2)/(2*(gam^2+lambda^2)));
post_mean = 1 + int_weightz'*pointz;
C0 = lambda^2/(2*gam^2 + lambda^2);
post_var = (f_at_mode*det(T)/alph)^2*...
    (C0 - (int_weightz'/gram_mat)*int_weightz);
lower_bound = 0.25;
upper_bound = 1.05;
hold on
lap_dist = fplot(hSub2, @(x) normpdf(x, post_mean, sqrt(post_var)),...
    [lower_bound upper_bound], '-k','LineWidth', 1.5);
set(gca, 'FontSize', 18)
set(gca, 'TickLabelInterpreter', 'latex')
true_line = xline(1, '--r','LineWidth', 1.25);
post_line = xline(post_mean, '--b','LineWidth', 1.25);
quant_line = xline(norminv(0.025, post_mean, sqrt(post_var)),...
    ':k', 'LineWidth', 1.5);
quant_line = xline(norminv(0.975, post_mean, sqrt(post_var)),...
    ':k', 'LineWidth', 1.5);
legend([lap_dist, quant_line, true_line, post_line],...
    {'Integral posterior', '95\% interval',...
    'True integral value/LA', '$m_1$'},...
    'interpreter', 'latex', 'FontSize', 18);
title(gca, 'Posterior of $F$, $\mathcal{N}\left(m_1, C_1\right)$',...
    'Interpreter', 'latex', 'FontSize', 22);
hold off

