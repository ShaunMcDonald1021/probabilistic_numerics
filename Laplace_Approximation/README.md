Code and documentation for "A probabilistic diagnostic tool to assess Laplace approximations: proof of concept and non-asymptotic experimentation".

**DOCUMENTATION**

The main manuscript is `laplace_poc_draft.tex`. The other files are materials from various seminar/conference presentations/posters.

**CODE**

A combination of R and MATLAB code. **Note**: the MATLAB part requires the code for [fully symmetric kernel quadrature](https://github.com/tskarvone/fskq) (Karvonen and Särkkä, 2018), so make sure this has been cloned/downloaded and added to your MATLAB path before running the diagnostic.

Assuming you want to apply the diagnostic to, say, a TMB model object, use the R functions in `laplace_diag_functions.R` to extract the necessary components, calculate interrogation values (on the log scale), and save it all into a `.mat` file:

```
# Given: TMB object TMB_obj, preliminary interrogation grid s_star as n*d matrix

# Get necessary ingredients for diagnostic
# (various components used for Laplace approx., function to evaluate for interrogations, etc.)
diag_ingredients = lap_diag_from_tmb(TMB_obj)

# Calculate interrogation values and save these, along with ingredients calculated above
log_interrs = get_log_interrs(diag_ingredients$logf, diag_ingredients$T_mat, diag_ingredients$mode, s_star, write_mat = TRUE,
                              diag_ingredients$logf_at_mode, diag_ingredients$log_T_det, 'diag_stuff.mat')
```

Next, use MATLAB to calibrate the diagnostic (optionally) and run it on the function (`diag_calib.m` and `lap_diag.m`, respectively).
```
% Given: Interrogation values for function on log scale & all necessary ingredients for L.A. (saved in .mat file),
% values for hyperparameters lambda and gam, dimensionality d, degrees of freedom v for calibration function

% Load stuff calculated in R
load('diag_stuff.mat')

% Turn preliminary grid into cell array in order to use FSKQ
Us = fss_gen(s_star(:, [1 d + 1]));

% Need to transpose some things
s_star = s_star';
logf_interrs = logf_interrs';

% Calibrate diagnostic
[~, w, wce, alph] = diag_calib(gam, lambda, d, v, Us, s_star);

% Run diagnostic on function
[diag_mean, diag_var, ~] = lap_diag(logf_interrs, logf_at_mode, log_T_det, d, gam, alpha, Us, w, wce, s_star);
```

The other files in the main part of the `Code` directory are used for the various demos and figures in the manuscript. The `legacy_code` subfolder contains a few things that were relevant to earlier versions of the project, but are now obsolete.


