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

% Probabilities at which the posterior of every tau_{i,g,h} is
% summarised (reported by both block-adaptive estimators and by
% montecarlo/tau_diagnostics.m).
cfg.blp.tau_probs = [0.05 0.25 0.50 0.75 0.95];

% How the POSTERIOR MEAN of the IRF is computed from the sampler output
% in the two block-adaptive estimators:
%   'rao_blackwell' (default) : average the CONDITIONAL means
%       E[beta | tau, sigma2, y] over the retained draws.  Same estimand
%       as averaging the draws, strictly smaller Monte Carlo error (the
%       exact Gaussian draw noise is integrated out instead of being
%       simulated), and it collapses to the exact closed form when tau
%       is fixed.  This is what makes the RMSE comparison against the
%       closed-form FMAR baseline fair: otherwise the adaptive
%       estimators are charged for avoidable simulation noise that the
%       baseline, having no sampler, never pays.
%   'draw_mean' : average the IRF draws (the behaviour of every result
%       produced before this option existed; set it to reproduce them).
cfg.blp.point_estimate = 'rao_blackwell';

% ---------------------------------------------------------------------
% Horizon pooling of the block scales (estimate_blp_blockpooled.m)
% ---------------------------------------------------------------------
% Prior on the tau path of each (equation, block):
%     log tau_{i,g,h} = phi * log tau_{i,g,h-1} + e,  e ~ N(0, kappa^2),
%     tau_{i,g,1} ~ half-Cauchy(0, 1),
%     kappa^2     ~ IG(a_kappa, b_kappa).
% phi = 1 is a random walk (a first-difference shrinkage prior on log
% tau); phi < 1 adds pull-back towards tau = 1.  The IG(2, 0.1) default
% has prior mean 0.1 for kappa^2, i.e. a typical one-horizon move of
% about 0.3 in log tau -- loose enough to track a genuine escape, tight
% enough to remove most of the horizon-to-horizon sampling noise.
cfg.blp.pool.phi         = 1.0;
cfg.blp.pool.a_kappa     = 2.0;
cfg.blp.pool.b_kappa     = 0.1;
cfg.blp.pool.mh_step     = 0.5;   % initial single-site step for log tau
cfg.blp.pool.mh_step_level = 0.3; % initial level-move step (whole path)
cfg.blp.pool.n_tau_sweeps  = 5;   % Metropolis sweeps of tau per iteration
cfg.blp.pool.kappa_fixed = [];    % non-empty = hold kappa there

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
cfg.mc.dgp_list = {'correct', 'sparse', 'intermediate', 'dense'}; % dgp/
% Which estimators run inside run_montecarlo.  'pooled' (the
% horizon-pooled block-adaptive estimator) is the most expensive one;
% drop it from this list for quick exploratory grids.  The baseline
% four are always present so every earlier result stays comparable.
cfg.mc.estimators = {'lp', 'var', 'global', 'block', 'pooled'};
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
% DGP design constants (see dgp/)
% ---------------------------------------------------------------------
% Every DGP reads these through a local helper with the documented
% default, so a configuration struct that predates this block still
% reproduces the original designs exactly.
%   sparse_a31   : the single misspecified entry A3(3,1) of the SPARSE
%                  DGP.  0 turns the misspecification off (the sparse
%                  DGP then coincides with the correct one); more
%                  negative = stronger prior-data conflict.
%   interm_scale : column-1 strength of the INTERMEDIATE DGP's A3;
%   interm_load  : how that strength is spread over the K equations
%                  (default puts the sparse DGP's -0.30 in equation 3).
%   dense_scale  : multiplies the dense DGP's MA matrix M (0 = no
%                  misspecification, 1 = the original design).
cfg.dgp.sparse_a31   = -0.30;
cfg.dgp.interm_scale = -0.30;
cfg.dgp.interm_load  = [0.5 0.7 1.0];
cfg.dgp.dense_scale  = 1.0;

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
% What the FMAR baseline reports at h = 1:
%   'bvar' : the Bayesian VAR response (the published FMAR convention,
%            and the DEFAULT so every earlier result is reproduced);
%   'lp'   : an h = 1 local projection with the same FMAR prior
%            machinery as h >= 2.  The block-adaptive estimators always
%            run an LP at h = 1, so 'lp' is the setting under which the
%            h = 1 comparison is like-for-like and the tau = 1 nesting
%            is exact at EVERY h >= 1.  See estimate_blp_fmar.m item 7.
cfg.fmar.h1_mode = 'bvar';

% ---------------------------------------------------------------------
% Demonstration script switches (RUN_ME_FIRST.m)
% ---------------------------------------------------------------------
cfg.demo.run_montecarlo = true;       % run a small MC at the end
cfg.demo.mc_dgps        = {'sparse'}; % which DGP(s) in the demo MC
cfg.demo.save_figures   = true;       % print figures to results/

end
