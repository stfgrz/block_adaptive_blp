% RUN_FMAR_DEMO  Demonstration of the FMAR-mode estimator stack.
%
% PURPOSE
% -------
% Companion to RUN_ME_FIRST.m (which demonstrates the original
% prototype machinery and is left untouched).  This script runs the
% pipeline with cfg.mode = 'fmar', i.e. the published methodology of
% Ferreira, Miranda-Agrippino & Ricco (REStat 2025) as the baseline and
% the block-adaptive layer nested inside it:
%   1.  set up path and configuration, switch to FMAR mode;
%   2.  single sparse-misspecification dataset: estimate
%         LP, BVAR (GLP), BLP-FMAR (global), BLP-block (adaptive,
%         tau independent by horizon) and BLP-pooled (adaptive, log tau
%         smoothed across horizons),
%       print the selected lambda path and the posterior block scales
%       tau (block 1 should escape in the y_3 equation at short h);
%   3.  verify the tau = 1 nesting on this dataset EXACTLY, using the
%       samplers' conditional posterior mean rather than the Gibbs
%       average, for both adaptive estimators;
%   4.  optionally run a small Monte Carlo (sparse + correct DGP) and
%       print the full report: both integrated-RMSE conventions, the
%       bias-variance decomposition, interval performance and the tau
%       diagnostics (montecarlo/report_montecarlo.m).
%
% The demonstration runs with cfg.fmar.h1_mode = 'lp', so the global
% baseline and the adaptive estimators are the same class of object at
% h = 1 and the nesting is exact from h = 1 upwards.  Set it back to
% 'bvar' for the published FMAR convention.
%
% QUICK MODE: QUICK_DEMO = true; RUN_FMAR_DEMO  shrinks everything.
%
% NOTES
% -----
% The FMAR-mode interval at h >= 2 is the quasi-Bayesian Newey-West
% sandwich band of FMAR (frequentist coverage target), while the BVAR
% band is a Bayesian credible band -- both at cfg.ci_level.

% ----------------------------------------------------------------------
% Step 0: path, configuration, seed
% ----------------------------------------------------------------------
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'config'), fullfile(this_dir, 'dgp'), ...
        fullfile(this_dir, 'estimators'), fullfile(this_dir, 'priors'), ...
        fullfile(this_dir, 'samplers'), fullfile(this_dir, 'montecarlo'), ...
        fullfile(this_dir, 'plots'), fullfile(this_dir, 'tests'), ...
        fullfile(this_dir, 'utils'));
results_dir = fullfile(this_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

cfg = default_config();
cfg.mode = 'fmar';
cfg.fmar.h1_mode = 'lp';        % like-for-like at h = 1 (see the header)

if exist('QUICK_DEMO', 'var') && QUICK_DEMO
    fprintf('*** QUICK_DEMO mode: reduced Monte Carlo and chains ***\n');
    cfg.mc.n_rep        = 6;
    cfg.mc.gibbs_n_burn = 150;
    cfg.mc.gibbs_n_keep = 350;
    cfg.gibbs.n_burn    = 250;
    cfg.gibbs.n_keep    = 600;
    cfg.fmar.n_niw_draws = 150;
end

rng(cfg.seed, 'twister');
fprintf('FMAR mode: K = %d, p = %d, T = %d, H = %d, seed = %d\n\n', ...
        cfg.K, cfg.p, cfg.T, cfg.H, cfg.seed);

% ----------------------------------------------------------------------
% Step 1: sparse-misspecification dataset, all four estimators
% ----------------------------------------------------------------------
fprintf('[1/3] Sparse-misspecification dataset, FMAR stack...\n');
dgp = simulate_sparse_misspec_dgp(cfg);

bvar  = estimate_bvar_niw(dgp.Y, cfg);
lp    = estimate_lp(dgp.Y, cfg, bvar.b1n);
blp_f = estimate_blp_fmar(dgp.Y, cfg, bvar);
blp_b = estimate_blp_blockadaptive(dgp.Y, cfg, bvar, blp_f.lambda);
blp_p = estimate_blp_blockpooled(dgp.Y, cfg, bvar, blp_f.lambda);

fprintf('  BVAR tightness lambda        = %.3f (max |eig| = %.3f)\n', ...
        bvar.lambda, bvar.max_eig);
fprintf('  BLP lambda_h (h = 2..%d)     = ', cfg.H);
fprintf('%.2f ', blp_f.lambda(1, 2:end));  fprintf('\n');

% Posterior block scales in the y_3 equation (the misspecified one):
fprintf('\n  Posterior mean tau, equation y_3 (rows = blocks, cols = h = 1..8):\n');
tau3 = squeeze(blp_b.tau_mean(3, :, 1:min(8, cfg.H)));
for g = 1:size(tau3, 1)
    fprintf('    block %d: ', g);  fprintf('%6.2f ', tau3(g, :));  fprintf('\n');
end
fprintf('  (block %d is the truly misspecified one)\n', dgp.misspec_block);

fprintf('\n  The SAME scales under horizon pooling (log tau smoothed across h):\n');
tau3p = squeeze(blp_p.tau_mean(3, :, 1:min(8, cfg.H)));
for g = 1:size(tau3p, 1)
    fprintf('    block %d: ', g);  fprintf('%6.2f ', tau3p(g, :));  fprintf('\n');
end
fprintf('  smoothing scale kappa (posterior mean, eq y_3): ');
fprintf('%.2f ', blp_p.kappa_mean(3, :));  fprintf('\n');
fprintf('  P(tau > 1) for block %d in eq y_3, h = 1..%d:\n    independent: ', ...
        dgp.misspec_block, min(8, cfg.H));
fprintf('%5.2f ', squeeze(blp_b.p_tau_gt1(3, dgp.misspec_block, 1:min(8, cfg.H))));
fprintf('\n    pooled     : ');
fprintf('%5.2f ', squeeze(blp_p.p_tau_gt1(3, dgp.misspec_block, 1:min(8, cfg.H))));
fprintf('\n\n');

% IRF comparison figure for all responses:
est_list = { ...
  struct('name', 'LP',        'theta', lp.theta,        'lo', lp.lo,   'hi', lp.hi,   'show_band', false), ...
  struct('name', 'BVAR',      'theta', bvar.theta,      'lo', bvar.theta_lo, 'hi', bvar.theta_hi, 'show_band', false), ...
  struct('name', 'BLP-FMAR',  'theta', blp_f.theta_mean,'lo', blp_f.lo, 'hi', blp_f.hi, 'show_band', false), ...
  struct('name', 'BLP-block', 'theta', blp_b.theta_mean,'lo', blp_b.lo, 'hi', blp_b.hi, 'show_band', true), ...
  struct('name', 'BLP-pooled','theta', blp_p.theta_mean,'lo', blp_p.lo, 'hi', blp_p.hi, 'show_band', false)};
fig1 = plot_irfs(dgp.theta_true, est_list, cfg, ...
                 'FMAR mode: IRFs, sparse misspecification');
if cfg.demo.save_figures
    print(fig1, fullfile(results_dir, 'fig_fmar_irf_sparse.png'), '-dpng', '-r120');
    fprintf('  Figure saved to results/fig_fmar_irf_sparse.png\n\n');
end

% ----------------------------------------------------------------------
% Step 2: nesting check on this dataset
% ----------------------------------------------------------------------
fprintf('[2/3] Nesting check (tau fixed at 1 vs FMAR posterior mean)...\n');
cfg_fix = cfg;  cfg_fix.blp.fix_tau = 1;
blp_fix  = estimate_blp_blockadaptive(dgp.Y, cfg_fix, bvar, blp_f.lambda);
blp_fixp = estimate_blp_blockpooled(dgp.Y, cfg_fix, bvar, blp_f.lambda);
% .theta_cond is the CONDITIONAL posterior mean, which with tau fixed is
% the exact closed form and carries no Monte Carlo error at all -- so
% this is a numerical identity, not a "close enough" check.
h_lo = 2;  if strcmp(cfg.fmar.h1_mode, 'lp'), h_lo = 1; end
cols = h_lo + 1:cfg.H + 1;
nest_err  = max(max(abs(blp_fix.theta_cond(:, cols)  - blp_f.theta_mean(:, cols))));
nest_errp = max(max(abs(blp_fixp.theta_cond(:, cols) - blp_f.theta_mean(:, cols))));
fprintf('  max |adaptive(tau=1) - BLP-FMAR| over h >= %d: %.2e (independent), %.2e (pooled)\n\n', ...
        h_lo, nest_err, nest_errp);

% ----------------------------------------------------------------------
% Step 3: small Monte Carlo (optional)
% ----------------------------------------------------------------------
if cfg.demo.run_montecarlo
    fprintf('[3/3] Monte Carlo in FMAR mode...\n');
    for dgps = {'sparse', 'correct'}
        dname = dgps{1};
        mc = run_montecarlo(cfg, dname);
        s  = summarize_montecarlo(mc);
        save(fullfile(results_dir, sprintf('mc_fmar_%s.mat', dname)), 'mc', 's');

        export_montecarlo_csv(s, fullfile(results_dir, sprintf('mc_fmar_%s', dname)));
        report_montecarlo(s);
    end
else
    fprintf('[3/3] Monte Carlo skipped (cfg.demo.run_montecarlo = false).\n');
end

fprintf('\nRUN_FMAR_DEMO finished.\n');
