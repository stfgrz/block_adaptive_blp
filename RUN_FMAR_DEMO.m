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
%         LP, BVAR (GLP), BLP-FMAR (global), BLP-block (adaptive),
%       print the selected lambda path and the posterior block scales
%       tau (block 1 should escape in the y_3 equation at short h);
%   3.  verify the tau = 1 nesting on this dataset (block estimator
%       with tau fixed at 1 reproduces the FMAR posterior mean);
%   4.  optionally run a small Monte Carlo (sparse + correct DGP) and
%       print bias / RMSE / coverage / detection summaries.
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
fprintf('  (block %d is the truly misspecified one)\n\n', dgp.misspec_block);

% IRF comparison figure for all responses:
est_list = { ...
  struct('name', 'LP',        'theta', lp.theta,        'lo', lp.lo,   'hi', lp.hi,   'show_band', false), ...
  struct('name', 'BVAR',      'theta', bvar.theta,      'lo', bvar.theta_lo, 'hi', bvar.theta_hi, 'show_band', false), ...
  struct('name', 'BLP-FMAR',  'theta', blp_f.theta_mean,'lo', blp_f.lo, 'hi', blp_f.hi, 'show_band', false), ...
  struct('name', 'BLP-block', 'theta', blp_b.theta_mean,'lo', blp_b.lo, 'hi', blp_b.hi, 'show_band', true)};
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
blp_fix = estimate_blp_blockadaptive(dgp.Y, cfg_fix, bvar, blp_f.lambda);
nest_err = max(max(abs(blp_fix.theta_mean(:, 3:end) - blp_f.theta_mean(:, 3:end))));
fprintf('  max |blockadaptive(tau=1) - BLP-FMAR| over h >= 2: %.4f\n\n', nest_err);

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

        fprintf('\n--- FMAR-mode Monte Carlo summary, DGP = %s (R = %d) ---\n', dname, s.R);
        fprintf('Integrated RMSE (mean over horizons), per response variable:\n');
        fprintf('  %-10s', 'estimator');
        for i = 1:cfg.K, fprintf('     y_%d ', i); end
        fprintf('\n');
        for e = 1:4
            fprintf('  %-10s', s.est_names{e});
            fprintf('  %7.3f', s.irmse(e, :));
            fprintf('\n');
        end
        fprintf('Average coverage of %d%% intervals (over horizons, response y_%d):\n', ...
                round(100 * cfg.ci_level), cfg.K);
        for e = 1:4
            fprintf('  %-10s  %5.2f\n', s.est_names{e}, ...
                    mean(squeeze(s.coverage(e, cfg.K, :))));
        end
        if ~isnan(s.detect_prob)
            fprintf('P(truly misspecified block %d has largest scale) = %.2f\n', ...
                    s.misspec_block, s.detect_prob);
        end
    end
else
    fprintf('[3/3] Monte Carlo skipped (cfg.demo.run_montecarlo = false).\n');
end

fprintf('\nRUN_FMAR_DEMO finished.\n');
