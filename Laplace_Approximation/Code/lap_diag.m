function [post_mean, post_var, point_contribs]...
    = lap_diag(logf_inter, logf_at_mode, logdet_T, d, gam, alph,...
    Us, w_or_lambda, wce, s_star, options)

% Runs the LA diagnostic for a given function. Make sure you download FSKQ
% prior to running this (see README).
arguments
    % n*1 vector of interrogation values on log scale
    logf_inter double
    % value of log(f) at mode
    logf_at_mode (1,1) double
    % log-determinant of T matrix (see Section 4.1 of manuscript)
    logdet_T (1,1) double
    % Dimensionality of domain
    d (1,1) {mustBePositive, mustBeInteger}
    % Hyperparameters gamma and alpha (see Sections 4-5 of manuscript)
    gam (1,1) {mustBePositive}
    alph (1,1) {mustBePositive}
    % The preliminary interrogation grid as a cellular array. Will
    % usually be obtained by running fss_gen (from FSKQ) on a matrix of
    % generator vectors (see demo.m in FSKQ)
    Us cell
    % Either a vector of BQ weights (e.g. precomputed with diag_calib) or 
    % a single positive number giving the lambda hyperparameter, in which
    % case weights will be calculated to get BQ
    w_or_lambda double
    % Unscaled posterior integral standard deviation (presumably computed 
    % w/ diag_calib alongside w). Only supply if is_w == TRUE, otherwise
    % wce is computed alongside weights
    wce double
    % Preliminary grid in the form of n*d array. Not necessary, but
    % saves time if precomputed
    s_star double = cell2mat(Us)'
    % Whether or not integral diagnostic should be plotted
    options.should_plot logical = true
    % Whether w_or_lambda should be interpreted as w or lambda.
    options.is_w logical = true
    % Whether or not to assume s_star is fully symmetric
    options.is_symm logical = true
end

% OUTPUTS:
% post_mean, post_var: Posterior integral mean and variance, respectively
% point_contribs: Contributions to integral mean by points at each distance
% from mode.

% Re-weight interrogations and subtract prior mean interrogation values
Y = exp(d*log(gam) + logdet_T + logf_inter - log(mvnpdf(s_star/gam)))-...
    exp(d*log(2*pi*gam) + logdet_T + logf_at_mode +...
    log(mvnpdf(sqrt(gam^2-1)*s_star/gam)));

% Calculate BQ weights and wce if not precomputed
if ~options.is_w
    % Kernel stuff
    k = @(r)exp(-r.^2/(2*w_or_lambda^2));
    kmean = @(x)(w_or_lambda^2/(gam^2+w_or_lambda^2))^(d/2)*...
        exp(-norm(x)^2/(2*(gam^2+w_or_lambda^2)));
    Ikmean = (w_or_lambda^2/(2*gam^2 + w_or_lambda^2))^(d/2);

    if options.is_symm
        % Fully-symmetric grid allows substantial speedup
        [Q, wce, ~] = kq_fss(Y, Us, k, kmean, Ikmean);
    else
        [Q, wce, ~] = kq(Y, s_star', k, kmean, Ikmean);
    end
else
    % This part taken from the kq_fss.m source code in FSKQ
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

% Laplace approximation (prior integral mean)
lap_app = exp(logdet_T + logf_at_mode + d*log(2*pi)/2);
% Posterior integral mean (LA + BQ calculation)
post_mean = Q + lap_app;
% Posterior integral variance
post_var = wce^2*exp(2*(logf_at_mode + logdet_T) - d*log(alph));

% Plot of diagnostic
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
