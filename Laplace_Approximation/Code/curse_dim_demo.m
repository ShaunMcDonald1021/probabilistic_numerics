% This script generates Figure 6 from manuscript, showing the discrepancy
% between shape differences and mass differences in high dimensions
d = 72;
v = 25921;

% (Radial) CDF of T density
pT = @(r) fcdf(r.^2/d, d, v);
% (Radial) cumulative mass of its Gaussian approximation
pGauss = @(r) exp(gammaln((v+d)/2) - gammaln(v/2)+...
    d*log(2/(v+d))/2)*chi2cdf((v+d)*r.^2./v, d);

% T density as function of distance from origin
tau = @(r) mvtpdf([r(:) zeros([max(size(r)) d-1])], eye(d), v);
% Its Gaussian approximation, also as fn. of distance/radius
gauss_approx = @(r) exp(gammaln((v+d)/2) - gammaln(v/2) - d*log(v*pi)/2-...
   (v+d)*r.^2./(2*v));

% Difference in (radial) cumulative mass between tau and its Gaussian
% approx.
fig = tiledlayout(2, 1, 'TileSpacing', 'Tight');
hsub1 = nexttile;
hold on
tline = fplot(hsub1, pT, [0 12], 'LineWidth', 1.75);
set(gca, 'FontSize', 14)
set(gca, 'TickLabelInterpreter', 'latex')
gaussline = fplot(pGauss, [0 12], 'Color', 'r', 'Linewidth', 1.25);
ylim([0 1.05]);
h_leg = legend([tline, gaussline],...
    {'$\tau$',...
    '$m^x_0 \cdot g$'},...
    'Location', 'northwest', 'interpreter', 'latex', 'FontSize', 16);
HeightScaleFactor = 2.5;
NewHeight = h_leg.Position(4) * HeightScaleFactor;
h_leg.Position(2) = h_leg.Position(2) - (NewHeight - h_leg.Position(4));
h_leg.Position(4) = NewHeight;
set(gca, 'xtick', [])
ylabel('$$\int_{|\!|x|\!|< r} f(x)\mathrm{d}x$$ ', 'interpreter', 'latex')
hold off

% Difference between shape of tau and its Gaussian approximation
hsub2 = nexttile;
hold on
set(gca, 'FontSize', 14)
set(gca, 'TickLabelInterpreter', 'latex')
diffline = fplot(hsub2, @(r) (tau(r) - gauss_approx(r))./tau(0),...
    [0 12], 'LineWidth', 1.5);
ylim([0 2.1e-05]);
ylabel('$$\frac{\tau(r) - m^x_0(r) \cdot g(r)}{\tau(0)}$$',...
'interpreter', 'latex');
xlabel('$r$', 'interpreter', 'latex');
hold off
