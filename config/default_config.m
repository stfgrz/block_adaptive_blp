function cfg = default_config()
% PURPOSE
% -------
% Central configuration for the block-adaptive Bayesian Local Projection
% (BLP) prototype. Every numerical choice of the prototype lives here so
% that experiments are reproducible and easy to modify. All downstream
% functions receive this struct; none of them hard-code these values.
%
% MODEL / EQUATIONS
% -----------------
% The prototype studies local projections (LP) at horizons h = 1,...,H,
%
%     y(t+h) = Beta_h' * z(t) + u(t+h),   z(t) = [1; y(t); ...; y(t-p+1)],
%
% with a prior that centres Beta_h on the LP coefficients implied by an
% estimated VAR(p), and block-specific "escape" scales tau_{g,h}:
%
%     beta_{g,h} - mu_{g,h}^VAR ~ N(0, lambda_h^2 * tau_{g,h}^2 * D_{g,h}).
%
% INPUTS
% ------
% (none)
%
% OUTPUTS
% -------
% cfg : struct with fields documented inline below.
%
% DIMENSIONS
% ----------
% K = number of variables, p = VAR lag order used by the ESTIMATED models,
% T = sample size, H = maximum LP horizon.
%
% NOTES
% -----
% This is a pedagogical prototype configuration. It is NOT the
% configuration of the published Ferreira/Miranda-Agrippino/Ricco BLP.

% ---------------------------------------------------------------------
% Core model dimensions
% ---------------------------------------------------------------------
cfg.K      = 3;      % number of endogenous variables
cfg.p      = 2;      % lag order of the ESTIMATED VAR / LP controls
cfg.T      = 200;    % sample size kept after burn-in
cfg.burnin = 100;    % burn-in periods discarded in DGP simulation
cfg.H      = 20;     % maximum LP horizon (IRFs reported for h = 0..H)

% Structural shock of interest: recursive (Cholesky) shock to variable 1.
% NOTE: recursive identification is used here ONLY as a transparent
% simulation device.  It is not the intended final empirical strategy.
cfg.shock_var = 1;

% Master random seed. Individual Monte Carlo replications derive their
% own seeds from this (see run_montecarlo.m) so runs are reproducible.
cfg.seed = 12345;

% ---------------------------------------------------------------------
% Coefficient blocks
% ---------------------------------------------------------------------
% 'per_variable': block g collects, within each LP equation, ALL lag
% coefficients that multiply variable g, i.e. the coefficients on
% y_g(t), y_g(t-1), ..., y_g(t-p+1).  This gives G = K blocks of p
% coefficients each.  The intercept is NOT shrunk by any tau (block 0).
cfg.blocks.scheme = 'per_variable';

% ---------------------------------------------------------------------
% Prior / global tightness (lambda_h)
% ---------------------------------------------------------------------
% lambda_mode:
%   'fixed' : use cfg.blp.lambda_fixed for every equation and horizon.
%   'grid'  : for each (equation, horizon) choose lambda on a grid by
%             maximising the Gaussian marginal likelihood of the LP
%             regression, conditional on sigma^2 = OLS estimate and all
%             tau_g = 1 (see priors/select_global_lambda.m).
% This module is deliberately isolated so it can later be replaced by
% the global-tightness procedure of the published BLP.
cfg.blp.lambda_mode  = 'grid';
cfg.blp.lambda_fixed = 0.30;
cfg.blp.lambda_grid  = [0.025 0.05 0.10 0.20 0.40 0.80 1.60];

% Loose prior standard deviation multiplier for the intercept:
% prior sd of intercept in equation i = intercept_scale * std(y_i).
cfg.blp.intercept_scale = 10;

% If true, return full posterior draws inside estimator outputs
% (memory heavy in Monte Carlo; keep false there).
cfg.blp.return_draws = false;

% ---------------------------------------------------------------------
% Gibbs sampler settings (single-dataset demonstrations)
% ---------------------------------------------------------------------
cfg.gibbs.n_burn  = 500;    % discarded warm-up draws
cfg.gibbs.n_keep  = 1500;   % retained posterior draws
cfg.gibbs.a0      = 0.01;   % weak Inverse-Gamma(a0, b0) prior on sigma^2
cfg.gibbs.b0      = 0.01;
% Numerical guard rails on tau_g^2 draws (wide truncation of the
% heavy-tailed horseshoe; documented in gibbs_block_horseshoe.m):
cfg.gibbs.tau2_min = 1e-10;
cfg.gibbs.tau2_max = 1e8;

% ---------------------------------------------------------------------
% Ordinary LP inference
% ---------------------------------------------------------------------
% 'nw'  : Newey-West HAC with truncation lag = h (simple prototype rule)
% 'ols' : homoskedastic iid covariance (understates uncertainty for h>1)
% The covariance computation is isolated in estimate_lp.m > lp_vcov()
% so that a more careful HAC / small-sample procedure can replace it.
cfg.lp.vcov = 'nw';

% ---------------------------------------------------------------------
% VAR interval settings
% ---------------------------------------------------------------------
% Intervals for VAR IRFs are computed by simulating coefficient draws
% from the asymptotic normal distribution of the OLS estimator
% (Sigma held at its point estimate; documented in estimate_var.m).
cfg.var.n_ci_sim = 500;

% ---------------------------------------------------------------------
% Interval nominal level (used by ALL estimators, frequentist and
% Bayesian, so coverage numbers are comparable)
% ---------------------------------------------------------------------
cfg.ci_level = 0.90;

% ---------------------------------------------------------------------
% Monte Carlo experiment
% ---------------------------------------------------------------------
cfg.mc.n_rep    = 100;                            % demonstration size
cfg.mc.dgp_list = {'correct', 'sparse', 'dense'}; % see dgp/ folder
% Cheaper sampler settings inside the Monte Carlo loop:
cfg.mc.gibbs_n_burn = 300;
cfg.mc.gibbs_n_keep = 700;
% Use parfor if Parallel Computing Toolbox present? (never required)
cfg.mc.use_parfor = false;

% ---------------------------------------------------------------------
% True-IRF computation for non-analytic DGPs
% ---------------------------------------------------------------------
% Number of simulated shocked/unshocked path pairs with common random
% numbers.  For the LINEAR DGPs in this prototype the analytic IRF is
% available and the simulation route is only used as a cross-check.
cfg.true_irf.n_paths = 1000;
cfg.true_irf.t_max_extra = 5;   % simulate a few periods beyond H

% ---------------------------------------------------------------------
% Methodology mode
% ---------------------------------------------------------------------
% 'prototype' : original semi-conjugate machinery (grid lambda, raw
%               data, posterior-quantile bands) -- kept as the verified
%               reference and for backward compatibility of all early
%               results.
% 'fmar'      : published Ferreira / Miranda-Agrippino / Ricco (REStat
%               2025) machinery: GLP Bayesian VAR, detrended horizon
%               regressions, NIW prior with Newey-West long-run scales,
%               marginal-likelihood tightness with horizon-dependent
%               Gamma hyperprior, quasi-Bayesian NW sandwich bands.
%               The block-adaptive estimator is then EXACTLY nested in
%               the FMAR baseline at tau = 1 (tests/test_fmar_nesting.m).
cfg.mode = 'prototype';

% ---------------------------------------------------------------------
% FMAR-mode settings (used only when cfg.mode = 'fmar'; values are the
% FMAR replication-code defaults)
% ---------------------------------------------------------------------
cfg.fmar.Vc         = 1e5;      % prior variance of the constant (their lambdaC)
cfg.fmar.isrw       = false;    % random-walk prior centre for the BVAR.
                                % Scalar (all K variables alike) or a 1 x K
                                % vector.  false = white-noise centre, right
                                % for the stationary simulated DGPs; true for
                                % all-levels data.  Mixed systems need the
                                % vector form -- the Ch. 7 empirical system
                                % uses [0 1 1 1 1], a white-noise centre for
                                % the policy surprise and random walks for
                                % the four levels (see empirical/).
cfg.fmar.lambda_min = 1e-4;     % FMAR/GLP optimisation bounds for lambda
cfg.fmar.lambda_max = 5;
cfg.fmar.hyper_mode = 0.4;      % Gamma hyperprior mode (all horizons)
% Horizon-dependent hyperprior sd (FMAR's logistic rule): loosens the
% hyperprior with the horizon, from ~0.11 at h = 1 to ~0.5 at long h.
cfg.fmar.hyper_sd_rule = @(h) 0.1 + 0.4 ./ (1 + exp(-0.3 * (h - 12)));
cfg.fmar.n_niw_draws = 500;     % NIW posterior draws for BVAR bands

% ---------------------------------------------------------------------
% Demonstration script switches (RUN_ME_FIRST.m)
% ---------------------------------------------------------------------
cfg.demo.run_montecarlo = true;       % run a small MC at the end
cfg.demo.mc_dgps        = {'sparse'}; % which DGP(s) in the demo MC
cfg.demo.save_figures   = true;       % print figures to results/

end
