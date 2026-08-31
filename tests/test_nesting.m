function test_nesting()
% PURPOSE
% -------
% Verify the NESTING property: when every local block scale is FIXED at
% tau_g = 1, the block-adaptive prior
%     beta - mu ~ N(0, lambda^2 * tau_g^2 * D)
% collapses to the global prior N(0, lambda^2 D), so the two estimators
% target the SAME posterior.  Run with INDEPENDENT random seeds, their
% posterior summaries must agree up to Monte Carlo (sampling) error.
%
% MODEL / EQUATIONS
% -----------------
% Both estimators use the same fixed lambda, same data, same prior
% centre; the block-adaptive run sets cfg.blp.fix_tau = 1, which routes
% through the SAME Gibbs code with sample_tau = false, but is seeded
% differently, so agreement is a Monte-Carlo-error statement, not a
% trivial identity of the random stream.
%
% Tolerances: (a) a HIGH-PRECISION single-regression comparison with
% 8000 kept draws must agree within 6 Monte Carlo standard errors
% (MCSE approximated by sd/sqrt(n_keep/10), a conservative effective-
% sample-size deflation); (b) the full IRF surfaces with the default
% draw count must agree within a loose absolute tolerance.
%
% INPUTS / OUTPUTS
% ----------------
% None; prints PASS or raises an assertion error.
%
% NOTES
% -----
% Any systematic disagreement here means the two estimators do not
% implement the same prior family -- a bug, not sampling noise.

fprintf('--- test_nesting ---\n');
cfg = default_config();
cfg.blp.lambda_mode  = 'fixed';
cfg.blp.lambda_fixed = 0.30;
cfg.blp.return_draws = false;

rng(cfg.seed, 'twister');
dgp = simulate_var_dgp(cfg);
var_est = estimate_var(dgp.Y, cfg);

% ---- (a) High-precision comparison on one equation/horizon -----------
h = 6;  i = 2;
[T, K] = size(dgp.Y);  p = cfg.p;
Zall = build_lp_regressors(dgp.Y, p);
Zh = Zall(1:end - h, :);
yh = dgp.Y(p + h:T, i);
bp = build_block_prior(dgp.Y, cfg);
Mu = var_implied_lp_prior(var_est, h);

prior.mu = Mu(:, i);  prior.d = bp.d_all(:, i);
prior.block_id = bp.block_id;
prior.lambda = cfg.blp.lambda_fixed;
prior.a0 = cfg.gibbs.a0;  prior.b0 = cfg.gibbs.b0;

opts_g = struct('n_burn', 1000, 'n_keep', 8000, 'sample_tau', false, ...
    'tau_fixed', 1, 'tau2_min', cfg.gibbs.tau2_min, ...
    'tau2_max', cfg.gibbs.tau2_max, 'seed', 111);
opts_b = opts_g;  opts_b.seed = 222;    % independent stream

out_g = gibbs_block_horseshoe(yh, Zh, prior, opts_g);
out_b = gibbs_block_horseshoe(yh, Zh, prior, opts_b);

th_g = out_g.beta_draws(:, 2:1 + K) * var_est.b1n;
th_b = out_b.beta_draws(:, 2:1 + K) * var_est.b1n;
mcse = std(th_g) / sqrt(opts_g.n_keep / 10);   % conservative ESS/10
gap  = abs(mean(th_g) - mean(th_b));
fprintf('  single-regression gap %.4f vs 6*MCSE %.4f\n', gap, 6 * mcse);
assert(gap < 6 * mcse, 'Nesting failed (high precision): gap %.4f.', gap);

% ---- (b) Full IRF surface with default draws --------------------------
rng(777, 'twister');
blp_glob = estimate_blp_global(dgp.Y, cfg, var_est);
rng(888, 'twister');
cfg_fix = cfg;  cfg_fix.blp.fix_tau = 1;
blp_fix = estimate_blp_blockadaptive(dgp.Y, cfg_fix, var_est, blp_glob.lambda);

d_all = abs(blp_glob.theta_mean(:, 2:end) - blp_fix.theta_mean(:, 2:end));
fprintf('  full-surface max gap %.4f, mean gap %.4f\n', ...
        max(d_all(:)), mean(d_all(:)));
assert(max(d_all(:)) < 0.05, 'Nesting failed: max IRF gap %.4f.', max(d_all(:)));
assert(mean(d_all(:)) < 0.01, 'Nesting failed: mean IRF gap %.4f.', mean(d_all(:)));

fprintf('PASS: test_nesting\n\n');
end
