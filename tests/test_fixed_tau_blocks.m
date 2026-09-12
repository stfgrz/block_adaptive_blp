function test_fixed_tau_blocks()
% PURPOSE
% -------
% Pin the behaviour of cfg.blocks.fixed_tau, the option that holds the
% escape scale of chosen blocks at 1 (the global FMAR prior) while the
% other blocks adapt.  The euro-area application uses it for the block of
% the policy surprise (see default_config.m for why).
%
% WHAT IS CHECKED, for BOTH adaptive estimators
% ---------------------------------------------
% 1. DEFAULT IS THE ORIGINAL SAMPLER.  With cfg.blocks.fixed_tau = [] the
%    output is bit-identical to a configuration that lacks the field.
% 2. HOLDING EVERY BLOCK equals the tau = 1 nesting case: the conditional
%    posterior mean is the FMAR closed form to machine precision, and
%    every reported tau is exactly 1.
% 3. HOLDING ONE BLOCK: that block's tau_mean is exactly 1 and its
%    P(tau > 1) exactly 0 at every equation and horizon, while the other
%    blocks still move (their tau is not identically 1), and the IRF
%    differs from both the all-free and the all-held runs.
% 4. INPUT VALIDATION: an out-of-range block index is refused.

fprintf('test_fixed_tau_blocks:\n');

cfg = default_config();
cfg.mode = 'fmar';
cfg.fmar.h1_mode = 'lp';
cfg.H = 8;  cfg.T = 180;
cfg.gibbs.n_burn = 80;  cfg.gibbs.n_keep = 160;
cfg.fmar.n_niw_draws = 40;

rng(6101, 'twister');
dgp  = simulate_sparse_misspec_dgp(cfg);
bvar = estimate_bvar_niw(dgp.Y, cfg);
blpf = estimate_blp_fmar(dgp.Y, cfg, bvar);
G = cfg.K;

ests = {'independent', @(c, s) run_ind(dgp.Y, c, bvar, blpf.lambda, s); ...
        'pooled',      @(c, s) run_pool(dgp.Y, c, bvar, blpf.lambda, s)};
for e = 1:size(ests, 1)
    nm = ests{e, 1};  f = ests{e, 2};

    % --- 1. default unchanged -------------------------------------------
    c_nofield = cfg;  c_nofield.blocks = rmfield(c_nofield.blocks, 'fixed_tau');
    c_empty   = cfg;  c_empty.blocks.fixed_tau = [];
    a = f(c_nofield, 6102);  b = f(c_empty, 6102);
    assert(isequal(a.theta_mean, b.theta_mean) && isequal(a.tau_mean, b.tau_mean), ...
        '%s: an empty fixed_tau changed the output', nm);

    % --- 2. all blocks held == tau = 1 nesting --------------------------
    c_all = cfg;  c_all.blocks.fixed_tau = 1:G;
    r_all = f(c_all, 6103);
    d = max(max(abs(r_all.theta_cond(:, 2:end) - blpf.theta_mean(:, 2:end))));
    assert(d < 1e-10, '%s: holding every block is not the FMAR closed form (%.2e)', nm, d);
    assert(all(r_all.tau_mean(:) == 1) && all(r_all.p_tau_gt1(:) == 0), ...
        '%s: held blocks do not report tau = 1 / P(tau>1) = 0', nm);

    % --- 3. one block held -----------------------------------------------
    c_one = cfg;  c_one.blocks.fixed_tau = 1;
    r_one = f(c_one, 6104);
    assert(all(all(r_one.tau_mean(:, 1, :) == 1)) && all(all(r_one.p_tau_gt1(:, 1, :) == 0)), ...
        '%s: the held block moved', nm);
    others = r_one.tau_mean(:, 2:end, :);
    assert(any(others(:) ~= 1), '%s: the free blocks did not move', nm);
    assert(isequal(r_one.fixed_tau_blocks, 1), '%s: fixed_tau_blocks not recorded', nm);
    d1 = max(max(abs(r_one.theta_mean(:, 2:end) - r_all.theta_mean(:, 2:end))));
    d2 = max(max(abs(r_one.theta_mean(:, 2:end) - a.theta_mean(:, 2:end))));
    assert(d1 > 1e-6 && d2 > 1e-6, ...
        '%s: holding one block did not produce a distinct estimator (%.1e, %.1e)', nm, d1, d2);
    fprintf('  %-11s: default unchanged; all held = FMAR to %.1e; one held: tau = 1 there, free elsewhere: OK\n', nm, d);
end

% --- 4. validation ---------------------------------------------------------
threw = false;
try
    c_bad = cfg;  c_bad.blocks.fixed_tau = G + 1;
    run_ind(dgp.Y, c_bad, bvar, blpf.lambda, 1);
catch
    threw = true;
end
assert(threw, 'an out-of-range fixed_tau index was accepted');
fprintf('  out-of-range block index refused: OK\n');

fprintf('PASS: test_fixed_tau_blocks\n\n');
end

% =====================================================================
function r = run_ind(Y, cfg, bvar, lam, seed)
rng(seed, 'twister');
r = estimate_blp_blockadaptive(Y, cfg, bvar, lam);
end

function r = run_pool(Y, cfg, bvar, lam, seed)
rng(seed, 'twister');
r = estimate_blp_blockpooled(Y, cfg, bvar, lam);
end
