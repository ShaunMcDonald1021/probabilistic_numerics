function [post_mean, post_var, point_contribs]...
    = lap_diag(logf_inter, logf_at_mode, logdet_T, d, gam, alph,...
    Us, w_or_lambda, wce, s_star, options)

% Runs the LA diagnostic for a given function. Make sure you download FSKQ
% prior to running this (see README).
arguments
    logf_inter double
    logf_at_mode (1,1) double
    logdet_T (1,1) double
    d (1,1)
    gam (1,1) {mustBePositive}
    alph (1,1) {mustBePositive}
    Us cell
    w_or_lambda double
    wce double
    s_star double = cell2mat(Us)'
    options.should_plot logical = true
    options.is_w logical = true
    options.is_symm logical = true
end

% INPUTS:
% logf_inter: logs of interrogation values $\log f(\boldsymbol{s})$. Should
% be n*1 (column) vector
% logf_at_mode: $\log f(\hat{x})$: function value at mode on log scale
% logdet_T: log of determinant of T matrix, obtained from eigen-
% decomposition of the Hessian (see Section 4.1 of manuscript)
% d: dimensionality
% gam, alph: hyperparameters $\gamma$ and $\alpha$ respectively.
% w_or_lambda: either a vector of BQ weights (e.g. precomputed with
% diag_calib) or a single positive number giving the lambda
% hyperparameter, in which case weights will be calculated to get BQ
% is_w: logical indicating whether w_or_lambda should be interpreted as w
% or lambda.
% should_plot: logical indicating whether a plot of the integral diagnostic
% should be made.
% Us: The preliminary grid as a cellular array. Will usually be obtained by
% running fss_gen (from FSKQ) on a matrix of generator vectors (see demo.m
% in FSKQ).
% wce: the unscaled posterior integral standard deviation (presumably
% computed with diag_calib, along with w). If lambda is supplied instead of
% w, then wce can be computed alongside weights, so only supply this if
% is_w = true
% s_star: the preliminary grid in the form of an n*d matrix. Not necessary,
% but it'll save time if you precomputed it

% OUTPUTS:
% post_mean, post_var: posterior integral mean and variance, respectively


% logf_inter(1)
% logf_inter(144)
% Re-weight interrogations and subtract prior mean interrogation values
Y = exp(d*log(gam) + logdet_T + logf_inter - log(mvnpdf(s_star/gam)))-...
    exp(d*log(2*pi*gam) + logdet_T + logf_at_mode +...
    log(mvnpdf(sqrt(gam^2-1)*s_star/gam)));
%Y(1) = 0;
%Y(find(Y > 1)) = 0;
%find(Y == max(Y))
%exp(d*log(gam) + logdet_T + logf_inter(1) - log(mvnpdf(s_star(1,:)/gam)))-...
%exp(d*log(2*pi*gam) + logdet_T + logf_at_mode +...
%    log(mvnpdf(sqrt(gam^2-1)*s_star(1,:)/gam)))

%Y(2855) = 0;
%[~,ind] = sort(Y, 'ascend');
%Y(ind(1:10))
%max(Y)
%ind(1)
%Y = Y(ind(3:end));
%s_star = s_star(ind(3:end),:);
%Y = Y([1 3:144 146:end]);
%s_star = s_star([1 3:144 146:end], :);
% Calculate BQ weights and wce if not precomputed
if ~options.is_w
    % Kernel stuff
    k = @(r)exp(-r.^2/(2*w_or_lambda^2));
    kmean = @(x)(w_or_lambda^2/(gam^2+w_or_lambda^2))^(d/2)*...
        exp(-norm(x)^2/(2*(gam^2+w_or_lambda^2)));
    Ikmean = (w_or_lambda^2/(2*gam^2 + w_or_lambda^2))^(d/2);

    if options.is_symm
        [Q, wce, ~] = kq_fss(Y, Us, k, kmean, Ikmean);
    else
        [Q, wce, ~] = kq(Y, s_star', k, kmean, Ikmean);
    end
else
    % This part stolen from the kq_fss.m source code in FSKQ
    Q = 0;
    J = length(w_or_lambda);
    Ls = zeros(J,1);
    point_contribs = zeros(J, 1);
    ind = 0;
    for i = 1:J
        Ls(i) = size(Us{i}, 2);
        Q = Q + w_or_lambda(i) * sum(Y(ind+1:ind+Ls(i)));
        point_contribs(i) = w_or_lambda(i)*sum(Y(ind+1:ind+Ls(i)));
        ind = ind + Ls(i);
    end
end

lap_app = exp(logdet_T + logf_at_mode + d*log(2*pi)/2);
post_mean = Q + lap_app;
post_var = wce^2*exp(2*(logf_at_mode + logdet_T) - d*log(alph));

% Plot
if options.should_plot
    lower_bound = min([norminv(0.025, post_mean, sqrt(post_var)) lap_app]);
    upper_bound = max([norminv(0.975, post_mean, sqrt(post_var)) lap_app]);
    lap_dist = fplot(@(x) normpdf(x, post_mean, sqrt(post_var)),...
        [lower_bound - 0.35*(upper_bound - lower_bound) upper_bound+...
        0.35*(upper_bound - lower_bound)], '-k','LineWidth', 1.5);
    set(gca, 'FontSize', 18)
    set(gca, 'TickLabelInterpreter', 'latex')
    hold on
    lap_line = xline(lap_app, '--b', 'LineWidth', 1.75);
    quant_line = xline(norminv(0.025, post_mean, sqrt(post_var)), ':k',...
        'LineWidth', 1.5);
    quant_line = xline(norminv(0.975, post_mean, sqrt(post_var)), ':k',...
        'LineWidth', 1.5);
    legend([lap_dist, quant_line, lap_line],...
        {'Integral posterior', '95\% interval',...
        'Laplace approximation'},...
        'interpreter', 'latex', 'FontSize', 18);
    title('Posterior of $F$, $\mathcal{N}\left(m_1, C_1\right)$', ...
        'Interpreter', 'latex', 'FontSize', 22);
    hold off
end

