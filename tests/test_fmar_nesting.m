function test_fmar_nesting()
% PURPOSE
% -------
% Verify the EXACT NESTING of the block-adaptive estimator in the FMAR
% baseline: with cfg.mode = 'fmar' and tau fixed at 1, the conditional
% posterior mean of beta in the conjugate-scaled sampler is
%     (Z'Z + W^{-1})^{-1} (Z'y + W^{-1} mu)
% for EVERY sigma2 draw (sigma2 cancels; see conditionals (1') in
% samplers/gibbs_block_horseshoe.m).  The Gibbs average of beta is
% therefore this constant vector up to pure Monte Carlo noise from the
% Gaussian draw step, and the implied IRFs must match
% estimate_blp_fmar's closed-form posterior mean tightly.
%
% Horizons h >= 2 are compared: at h = 1 the FMAR baseline reports the
% Bayesian VAR itself, while the block-adaptive estimator runs an
% h = 1 local projection (documented design difference).
%
% The comparison uses the SAME lambdas (blp_fmar.lambda passed into the
% block estimator), the same detrended data, the same prior centre and
% the same scales -- so any disagreement beyond Monte Carlo noise in
% the posterior-mean map would flag a port bug.
%
% TOLERANCE
% ---------
% The Gibbs mean of n_keep iid... (weakly dependent) draws around the
% fixed conditional mean has standard error ~ posterior sd / sqrt(ESS).
% With n_keep = 1500 draws a tolerance of 0.02 on IRF units (shock of
% size 1) is conservative for the sparse DGP's scales; the test also
% checks the MAXIMUM deviation across all variables/horizons.

fprintf('test_fmar_nesting: ');

cfg = default_config();
cfg.mode = 'fmar';
cfg.gibbs.n_burn = 500;
cfg.gibbs.n_keep = 1500;
cfg.blp.return_draws = false;

rng(777, 'twister');
dgp = simulate_sparse_misspec_dgp(cfg);
Y = dgp.Y;

bvar = estimate_bvar_niw(Y, cfg);
blp_f = estimate_blp_fmar(Y, cfg, bvar);

cfg_fix = cfg;
cfg_fix.blp.fix_tau = 1;
rng(778, 'twister');
blp_b = estimate_blp_blockadaptive(Y, cfg_fix, bvar, blp_f.lambda);

K = cfg.K;  H = cfg.H;
err = abs(blp_b.theta_mean(:, 3:H + 1) - blp_f.theta_mean(:, 3:H + 1)); % h >= 2
max_err = max(err(:));

assert(max_err < 0.02, ...
    'nesting violated: max |blockadaptive(tau=1) - FMAR| = %.4f', max_err);

% The prior centres must agree exactly (same BVAR, same mapping).
pc_err = max(max(abs(blp_b.prior_theta - blp_f.prior_theta)));
assert(pc_err < 1e-10, 'prior centres differ: %.3g', pc_err);

% And with tau sampled, the estimator must move AWAY from the global
% solution somewhere (sanity that the switch actually does something).
rng(779, 'twister');
blp_free = estimate_blp_blockadaptive(Y, cfg, bvar, blp_f.lambda);
dev = max(max(abs(blp_free.theta_mean(:, 3:end) - blp_f.theta_mean(:, 3:end))));
assert(dev > 1e-4, 'tau sampling appears inert (max dev %.2g)', dev);

fprintf('PASSED (max nesting error %.4f at h>=2; free-tau deviation %.3f)\n', ...
    max_err, dev);
end
