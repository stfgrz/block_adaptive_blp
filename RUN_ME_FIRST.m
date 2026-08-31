% RUN_ME_FIRST  Master demonstration of the block-adaptive BLP prototype.
%
% PURPOSE
% -------
% Reading and running this one script shows the entire workflow:
%   1.  set up the path and configuration;
%   2.  simulate a correctly specified VAR(2) dataset;
%   3.  estimate all four estimators (LP, VAR, global BLP, block-
%       adaptive BLP) and plot the IRF comparison;
%   4.  repeat on a SPARSE-misspecification dataset (omitted delayed
%       effect of the shock variable y_1 on y_3) and display the
%       posterior block scales: block 1 should escape the VAR prior,
%       above all in the y_3 equation; blocks 2-3 should not;
%   5.  optionally run a small Monte Carlo demonstration and summarise
%       bias / RMSE / coverage / block-detection metrics.
%
% Figures are saved to results/ when cfg.demo.save_figures is true.
%
% QUICK MODE: to smoke-test everything fast, set in the base workspace
%     QUICK_DEMO = true; RUN_ME_FIRST
% which shrinks the Monte Carlo and the Gibbs chains. Default is the
% full demonstration (cfg.mc.n_rep = 100 replications).
%
% NOTES
% -----
% This script does not call clear (it would delete QUICK_DEMO); start
% from a fresh workspace if in doubt.

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

if exist('QUICK_DEMO', 'var') && QUICK_DEMO
    fprintf('*** QUICK_DEMO mode: reduced Monte Carlo and chains ***\n');
    cfg.mc.n_rep        = 6;
    cfg.mc.gibbs_n_burn = 150;
    cfg.mc.gibbs_n_keep = 350;
    cfg.gibbs.n_burn    = 250;
    cfg.gibbs.n_keep    = 600;
    cfg.var.n_ci_sim    = 150;
end

rng(cfg.seed, 'twister');
fprintf('Configuration: K = %d, p = %d, T = %d, H = %d, seed = %d\n\n', ...
        cfg.K, cfg.p, cfg.T, cfg.H, cfg.seed);

% ----------------------------------------------------------------------
% Step 1: correctly specified DGP -- simulate and estimate
% ----------------------------------------------------------------------
fprintf('[1/5] Simulating correctly specified VAR(2) data...\n');
dgp1 = simulate_var_dgp(cfg);
fprintf('      %s\n', dgp1.description);

fprintf('[2/5] Estimating all four estimators on the correct DGP...\n');
var1 = estimate_var(dgp1.Y, cfg);
fprintf('      VAR stability: max|eig| = %.3f (stable = %d)\n', ...
        var1.max_eig, var1.is_stable);
lp1  = estimate_lp(dgp1.Y, cfg, var1.b1n);
bg1  = estimate_blp_global(dgp1.Y, cfg, var1);
bb1  = estimate_blp_blockadaptive(dgp1.Y, cfg, var1, bg1.lambda);

est_list1 = { ...
  struct('name', 'LP',        'theta', lp1.theta,      'lo', lp1.lo, 'hi', lp1.hi, 'show_band', false), ...
  struct('name', 'VAR',       'theta', var1.theta,     'lo', var1.theta_lo, 'hi', var1.theta_hi, 'show_band', false), ...
  struct('name', 'BLP-glob',  'theta', bg1.theta_mean, 'lo', bg1.lo, 'hi', bg1.hi, 'show_band', false), ...
  struct('name', 'BLP-block', 'theta', bb1.theta_mean, 'lo', bb1.lo, 'hi', bb1.hi, 'show_band', true)};
fig1 = plot_irfs(dgp1.theta_true, est_list1, cfg, ...
                 'IRFs, correctly specified VAR(2) DGP');
if cfg.demo.save_figures
    print(fig1, fullfile(results_dir, 'fig1_irf_correct.png'), '-dpng', '-r120');
end

% ----------------------------------------------------------------------
% Step 2: sparse-misspecification DGP -- simulate and estimate
% ----------------------------------------------------------------------
fprintf('\n[3/5] Simulating sparse-misspecification data...\n');
dgp2 = simulate_sparse_misspec_dgp(cfg);
fprintf('      %s\n', dgp2.description);

fprintf('[4/5] Estimating all four estimators on the sparse DGP...\n');
var2 = estimate_var(dgp2.Y, cfg);
fprintf('      Fitted VAR(2) stability: max|eig| = %.3f (stable = %d)\n', ...
        var2.max_eig, var2.is_stable);
lp2  = estimate_lp(dgp2.Y, cfg, var2.b1n);
bg2  = estimate_blp_global(dgp2.Y, cfg, var2);
bb2  = estimate_blp_blockadaptive(dgp2.Y, cfg, var2, bg2.lambda);

est_list2 = { ...
  struct('name', 'LP',        'theta', lp2.theta,      'lo', lp2.lo, 'hi', lp2.hi, 'show_band', false), ...
  struct('name', 'VAR',       'theta', var2.theta,     'lo', var2.theta_lo, 'hi', var2.theta_hi, 'show_band', false), ...
  struct('name', 'BLP-glob',  'theta', bg2.theta_mean, 'lo', bg2.lo, 'hi', bg2.hi, 'show_band', false), ...
  struct('name', 'BLP-block', 'theta', bb2.theta_mean, 'lo', bb2.lo, 'hi', bb2.hi, 'show_band', true)};
fig2 = plot_irfs(dgp2.theta_true, est_list2, cfg, ...
                 'IRFs, sparse misspecification (omitted lag-3 effect of y_1 on y_3)');
if cfg.demo.save_figures
    print(fig2, fullfile(results_dir, 'fig2_irf_sparse.png'), '-dpng', '-r120');
end

% ----------------------------------------------------------------------
% Step 3: posterior block scales on the sparse DGP
% ----------------------------------------------------------------------
fprintf('\n[5/5] Posterior block scales tau_{g,h} (sparse DGP):\n');
fprintf('      average posterior-mean tau by (equation, block), over horizons:\n');
fprintf('      %-14s', 'equation \ block');
for g = 1:cfg.K, fprintf('  block %d ', g); end
fprintf('\n');
for i = 1:cfg.K
    fprintf('      eq for y_%d    ', i);
    for g = 1:cfg.K
        fprintf('  %7.3f', mean(squeeze(bb2.tau_mean(i, g, :))));
    end
    fprintf('\n');
end
fprintf(['      Expectation under sparse misspecification: block %d should\n' ...
         '      carry the largest scales, especially in the equation for y_3.\n'], ...
         dgp2.misspec_block);

fig3 = plot_block_scales(bb2.tau_mean, cfg, ...
    'Posterior block scales, sparse DGP (thick line = truly misspecified block)', ...
    dgp2.misspec_block);
if cfg.demo.save_figures
    print(fig3, fullfile(results_dir, 'fig3_block_scales_sparse.png'), '-dpng', '-r120');
end

% ----------------------------------------------------------------------
% Step 4 (optional): small Monte Carlo demonstration
% ----------------------------------------------------------------------
if cfg.demo.run_montecarlo
    fprintf('\nMonte Carlo demonstration (this is the slow part)...\n');
    for id = 1:numel(cfg.demo.mc_dgps)
        dname = cfg.demo.mc_dgps{id};
        mc = run_montecarlo(cfg, dname);
        s  = summarize_montecarlo(mc);
        save(fullfile(results_dir, sprintf('mc_%s.mat', dname)), 'mc', 's');

        fprintf('\n--- Monte Carlo summary, DGP = %s (R = %d) ---\n', dname, s.R);
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
        fprintf('Share of replications with stable fitted VAR: %.2f\n', ...
                s.var_stable_share);

        fig4 = plot_rmse(s, cfg);
        fig5 = plot_block_scales(s.tau_bar, cfg, ...
            sprintf('MC-average posterior block scales, DGP: %s', dname), ...
            s.misspec_block);
        if cfg.demo.save_figures
            print(fig4, fullfile(results_dir, sprintf('fig4_rmse_%s.png', dname)), '-dpng', '-r120');
            print(fig5, fullfile(results_dir, sprintf('fig5_scales_%s.png', dname)), '-dpng', '-r120');
        end
    end
else
    fprintf('\n(Monte Carlo demonstration skipped; set cfg.demo.run_montecarlo = true.)\n');
end

fprintf('\nDone. Figures and .mat files are in results/.\n');
fprintf('Run tests with:  cd tests; run_all_tests\n');
