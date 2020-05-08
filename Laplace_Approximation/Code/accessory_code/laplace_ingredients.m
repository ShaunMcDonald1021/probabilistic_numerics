function [mode, f_sym, f_at_mode, hess_at_mode, args] =...
    laplace_ingredients(f, d, should_prompt)

% Uses symbolic toolkit to produce necessary ingredients for Laplace
% approximation of a function.

% INPUTS:
% f: an anonymous function handle for the full function from R^d -> R.
% It must take a single array-valued argument: a vertical concatenation
% of d-dimensional row vectors.
% d: dimension of the domain.
% should_prompt: logical indicating whether or not to prompt the user to
% input f_sym themselves. This may be a good idea in some applications (see
% below), but should be turned off for the interactive apps.

% OUTPUTS:
% mode: The location in the domain of the chosen mode (1xd sym row vector).
% f_sym: A symbolic version of f to be used for testing the diagnostic.
% f_at_mode: the value of f at mode, as a symbolic number.
% hess_at_mode: a dxd symbolic matrix for the Hessian of log(f) at mode.
% args: a d-dimensional sym row vector, to be used for calculations
% involving f_sym.

args = sym('t', [1 d]);
assumeAlso(args, 'real');

% Some target functions are inherently piecewise (e.g. things supported
% only on nonnegative numbers, like gamma pdf's). Conversion between
% symfun's and anonymous function handles doesn't work properly in this
% case, so we prompt the user to input their own f_sym.
if should_prompt
    f_sym = input(join(['Optional entry of symbolic version of f.' newline...
        'Press ENTER to skip, in which case it is created automatically.'...
        newline 'Manual entry is useful for piecewise f (e.g. a gamma'...
        "density) since those cases don't work when converting between sym"...
        'and matlab functions.' newline 'Manual entry MUST be defined in' ...
        'terms of symvar "args".'], ''));
end

% If f_sym entry is skipped, (attempt to) create it automatically
if ~should_prompt || size(f_sym, 1) == 0
    f_sym = f(args);
end

logf = log(f_sym);
hess = hessian(logf);

% Obtain all zeros of the gradient of logf
poss_modes = solve(gradient(logf));
% If there are multiple zeros, they are returned as a structure
% which we convert into a 2D-array
if isa(poss_modes, 'struct')
    poss_modes = struct2array(poss_modes);
end

num_poss = size(poss_modes, 1);

% First candidate for the mode
mode = poss_modes(1,:);
f_at_mode = subs(f_sym, args, mode);

% Cycle through the rest of the zeros. The one at which f
% attains the highest value is assumed to be a mode
% (this is presumably faster than checking negative-definiteness of
% the Hessian at each one)
for i = 2:num_poss
    cand_mode = poss_modes(i,:);
    cand_val = subs(f_sym, args, cand_mode);
    if isAlways(f_at_mode < cand_val)
        f_at_mode = cand_val;
        mode = cand_mode;
    end
end

hess_at_mode = subs(hess, args, mode);
end
