function [post_mean, post_var, lap_approx, true_int] ...
    = laplace_poc(f, a, b, d, s, lambda, alpha, gamma, input_true_int, make_plots)

% Conducts a diagnostic of the Laplace approximation using probabilistic
% numerics. The target function is given a GP prior with a decaying
% squared exponential covariance function (see documentation and
% accessory_code. We can then look at the resulting Normal posterior for
% the integral.
% NOTE: Internally, the domain of f is linearly transformed s.t.
% the approximating Gaussian function is proportional to a standard normal
% (i.e. centered at zero, and the negative Hessian of the log is the
% identity there). The interrogation points MUST be specified with respect
% to these transformed coordinates

% INPUTS:
% f: an anonymous function handle for the full function from R^d -> R.
% It must take a single array-valued argument: a vertical concatenation
% of d-dimensional row vectors.
% a, b: integration limits (possibly infinite) for integrated argument.
% Must be w.r.t. the same coordinates as s (see below).
% d: dimension of the domain.
% s: a set of n interrogation points of dimension d, stored in d*n matrix.
% They must be specified in terms of transformed coordinates (see above) -
% i.e. corresponding to numbers of standard deviations along the "principal
% axes". s determines whether or not we use the GPU - if it is a gpuArray,
% yes; if not, no.
% lambda, alpha: length-scale and precision hyperparameters for covariance,
% respectively.
% input_true_int: logical indicating if the user intends to input the true
% integral of f to check the diagnostic. Note that this is with respect to
% the UN-transformed scale.
% make_plots: logical indicating if diagnostic plots are to be made.

% OUTPUTS:
% post_mean, post_var: mean and variance, respectively, for the posterior
% distribution of the integral (w.r.t. the UNtransformed scale).
% lap_approx: the value of the Laplace approximation (the prior mean for
% the integral).
% true_int: the optional value for the true integral according to user
% input.

% Initialize (it stays this way if ~input_true_int)
true_int = NaN;

% Collect necessary components for Laplace approximation,
% and convert them from symbolic to numeric for faster computation
[mode, f_sym, f_at_mode, hess_at_mode, args] = laplace_ingredients(f, d,...
    true);
hess_at_mode = double(hess_at_mode);
mode = double(mode);
f_at_mode = double(f_at_mode);
%disp(f_sym);
% The coordinate transformation is based on an eigendecomposition of the
% Hessian.
if isa(s, 'gpuArray')
    [V, D] = eig(gpuArray(hess_at_mode));
    %f_at_mode = gpuArray(f_at_mode); % TODO: see if this speeds it up
    % Likewise for mode.
else
    [V, D] = eig(hess_at_mode);
end

% The matrix we use to rotate and scale the arguments, so that the negative 
% Hessian of the log is the identity at the mode.
rot_scale_mat = sqrt(-D)\V';
% NOTE: rot_scale_mat is not unique. In particular, the "banana" produces
% an "antidiagonal" matrix here, but it is obviously preferable to have a
% diagonal one. Consider changing it manually if need be.

% Need the absolute value of the determinant to put the integral posterior
% on the correct scale
rot_det = abs(det(rot_scale_mat));

% Rotate, scale, and translate the interrogation points to the original
% domain. The translation ensures the mode corresponds to the origin on
% the transformed scale.
s_trans = s'*rot_scale_mat + mode;

% Get value of Laplace approximation
lap_approx = rot_det*(f_at_mode*...
    (2^(1/2)*pi^(1/2)*(erf(2^(1/2)*b/2) - erf(2^(1/2)*a/2))/2)^d);

% Mean function for the GP prior on f (on un-transformed domain)
m0_t = @(x) f_at_mode*exp(sum(((x-mode)*hess_at_mode).*(x-mode), 2)/2);

% Function values at interrogation points
% NOTE: 2018 Shaun noted that function handle evaluations could be slower
% on the GPU than the CPU. May be worth investigating here
gauss_interr = m0_t(s_trans);
f_interr = f(s_trans);

% Covariance calculations
C0_t = se_decay_gram(s, lambda, gamma);
C0_cross = se_decay_cross_covar(s, lambda, gamma, a, b);

if gamma == Inf
    C0 = (4*lambda^3*pi^(1/2)*exp(-(a - b)^2/(4*lambda^2)) - ...
        4*lambda^3*pi^(1/2) + 2*lambda^2*pi*erf((a - b)/(2*lambda))*(a - b))^d;
else

    C0_multiple = ((4*gamma^2*lambda^2*pi^(3/2))/(2*gamma^2 + lambda^2)^(1/2))^d;

    % Calculate prior covariance for the integral. Save a bit of time
    % by noting it has a closed form if the integral is over all of R^d
    if a == -Inf && b == Inf
        C0 = C0_multiple;
    else
        C_sigma = 2*gamma^2*[lambda^2+gamma^2 gamma^2; gamma^2 lambda^2+gamma^2]/(lambda^2+2*gamma^2);
        C0 = C0_multiple*mvncdf([a a], [b b], [0 0], C_sigma)^d;
    end
end

weights = C0_cross/C0_t;
post_var = f_at_mode^2*rot_det^2*(C0 - weights*C0_cross')/alpha^d;

% Posterior mean
post_mean = lap_approx + rot_det*weights*(f_interr - gauss_interr);

% Prompt the user for an expression to calculate the true integral
% (e.g. by symbolic integration of f_sym w.r.t. args) on the
% original/untransformed scale
if input_true_int
    input_prompt = ["Input calculation of true integral (on the original"...
        "scale)." newline "You can use f_sym, a symbolic version of f "...
        "created automatically."];
        
    % If the bounds are finite, we need to transform them to the original
    % domain.
    if a > -Inf && b < Inf
        bounds= gather([repelem(a, d)*rot_scale_mat + mode;...
            repelem(b, d)*rot_scale_mat + mode]);
        upper_bounds = max(bounds);
        lower_bounds = min(bounds);
        input_prompt = [input_prompt newline "We've also taken the liberty "...
            "of providing transformed integral bounds in size-d row "...
            "vectors lower_bounds and upper_bounds."];
    end
    true_int = input(join(input_prompt, ''));

end

% Gather everything from GPU if applicable
if isa(s, 'gpuArray')
    [s_trans, post_mean, post_var, gauss_interr, f_interr] =...
        gather(s_trans, post_mean, post_var, gauss_interr, f_interr);
end
    
% Diagnostic plots
if make_plots
    % If d = 1 or 2, we can also show the full function along with its
    % Gaussian approximation and the interrogation points (on the
    % un-transformed scale).
    if d == 1
        figure(1)
        full_fun = fplot(f_sym, [min(s_trans)-0.5 max(s_trans)+0.5]);
        hold on
        % Temporarily suppress warning about irregular behaviour of m0_t
        % on arrays (which doesn't actually affect plotting behaviour)
        warning('off')
        gauss_approx = fplot(m0_t, [min(s_trans)-0.5 max(s_trans)+0.5]);
        warning('on');
        % Add in interrogation points
        for i = 1:size(s_trans, 1)
            point = s_trans(i);
            inter_line = plot([point(1) point(1)],...
                [gauss_interr(i) f_interr(i)], '-ok',...
                'LineWidth', 0.75, 'MarkerSize', 3, 'MarkerFaceColor', 'k');
        end
        
        legend([full_fun, gauss_approx, inter_line],...
        {'Full function', 'Gaussian approximation', 'Interrogation points'});
        title('True target function and Gaussian approximation');

        hold off
        figure(2)
        
    elseif d == 2
        % In 2-D, it's a surface plot
        figure(1)
        grid_bounds = [min(s_trans(:,1))-0.5 max(s_trans(:,1))+0.5...
            min(s_trans(:,2))-0.5 max(s_trans(:,2))+0.5];
    
        full_fun = fsurf(f_sym, grid_bounds, 'FaceAlpha', 0.65);
        
        % Set horizontal aspect ratio 1:1 for honesty's sake
        grid_ratio = (grid_bounds(4) - grid_bounds(3))/...
            (grid_bounds(2) - grid_bounds(1));
        pbaspect([1 grid_ratio geomean([grid_ratio 1/grid_ratio])])
        
        colormap winter
        axis vis3d
        hold on
        
        % This function handle returns NaN if m0_t < f_at_mode*10^-4, which
        % makes the plot look cleaner
        m0_plot = @(x,y) m0_t([x y]).*(m0_t([x y]) > f_at_mode*1e-4)./...
            (m0_t([x y]) > f_at_mode*1e-4);
        
        % Temporarily suppress warning about irregular behaviour of m0_t
        % on arrays (which doesn't actually affect plotting behaviour)
        warning('off')        
        gauss_approx = fsurf(m0_plot, grid_bounds,...
            'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'y');
        warning('on')
        
        % Add in interrogation points
        for i = 1:size(s_trans, 1)
            point = s_trans(i,:);
            inter_line = plot3([point(1) point(1)], [point(2) point(2)],...
                [gauss_interr(i) f_interr(i)], '-ok',...
                'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        end
        
        legend([full_fun, gauss_approx, inter_line],...
            {'Full function', 'Gaussian approximation', 'Interrogation points'});
        title('True target function and Gaussian approximation');
        xlabel('x'); ylabel('y');
        
        hold off
        figure(2)   
    end
    
    % Plot posterior of integral (can do for any d-value)
    lower_plot_bound = double(min([norminv(0.025, post_mean, sqrt(post_var))...
        true_int lap_approx])) - 0.5*sqrt(post_var);
    upper_plot_bound = double(max([norminv(0.975, post_mean, sqrt(post_var))...
        true_int lap_approx])) + 0.5*sqrt(post_var);
    lap_dist = fplot(@(x) normpdf(x, post_mean, sqrt(post_var)),...
        [lower_plot_bound upper_plot_bound], '-k');
    
    lap_line=xline(gather(lap_approx), '--b','LineWidth', 1);
    if ~isnan(true_int)
        true_line=xline(double(true_int), '--r','LineWidth', 1);
    end
    
    quant_line = xline(norminv(0.025, post_mean, sqrt(post_var)), ':k',...
       'LineWidth', 1.25);
    quant_line = xline(norminv(0.975, post_mean, sqrt(post_var)), ':k',...
       'LineWidth', 1.25);
   
    legend([lap_dist, quant_line, lap_line, true_line],...
        {'Posterior of integral', '95% region',...
        'Laplace approximation value', 'True integral value'});
    title('Laplace diagnostic plot');
end

end
