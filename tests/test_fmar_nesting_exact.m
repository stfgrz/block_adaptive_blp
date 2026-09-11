function test_fmar_nesting_exact()
% PURPOSE
% -------
% TOLERANCE-FREE verification that both block-adaptive estimators are
% EXACTLY nested in the FMAR baseline at tau = 1, across DGPs, seeds,
% EVERY equation and EVERY horizon.
%
% WHY A SECOND NESTING TEST
% -------------------------
% tests/test_fmar_nesting.m compares the GIBBS AVERAGE of the block
% estimator with the FMAR closed form.  That average carries Monte Carlo
% error, so the test needs a tolerance (0.02 in IRF units) and can only
% ever say "close enough".  A tolerance is exactly the wrong instrument
% for a nesting claim: it hides small systematic port bugs, and the
% temptation when it fails is to loosen it.
%
% With tau FIXED and the conjugate (sigma2-scaled) prior, conditional
% (1') of samplers/gibbs_block_horseshoe.m makes the conditional
% posterior mean of beta
%       mp = (Z'Z + W^{-1})^{-1} (Z'y + W^{-1} mu)
% INDEPENDENT of sigma2 and of the iteration, i.e. a deterministic
% function of the data and the prior.  The samplers now return it as
% .beta_cond_mean, and the estimators expose the implied IRF as
% .theta_cond.  Comparing THAT with estimate_blp_fmar's closed form is a
% comparison of two deterministic quantities: it should agree to machine
% precision, and this test asserts 1e-10, which is ~1e8 times tighter
% than the sampling-based test.  A genuine mismatch cannot hide.
%
% WHAT IS COVERED
% ---------------
%   * 4 datasets: correct, sparse, intermediate and dense DGPs;
%   * 2 seeds each (so 8 datasets in total);
%   * both cfg.fmar.h1_mode settings:
%       'bvar' -> nesting is claimed and checked for h = 2..H,
%       'lp'   -> nesting is claimed and checked for h = 1..H,
%     which is the point of the harmonised setting;
%   * all K equations and all H horizons, not an aggregate;
%   * BOTH the independent block-adaptive estimator and the
%     horizon-pooled one (a pooling prior on tau must not disturb the
%     tau = 1 special case).
%
% It also re-checks two things a pure equality test would miss:
%   * the prior centres agree exactly (same BVAR, same companion map);
%   * with tau sampled the estimators actually MOVE away from the global
%     solution, so "nesting holds" is not achieved by inertness.
%
% TOLERANCE
% ---------
% 1e-10 in IRF units, on the maximum over (equation, horizon).  Nothing
% here is statistical, so this is a numerical-linear-algebra tolerance,
% not a judgement call.

fprintf('test_fmar_nesting_exact:\n');

cfg0 = default_config();
cfg0.mode = 'fmar';
cfg0.H = 10;                 % full-length H adds cost, not coverage
cfg0.T = 180;
cfg0.gibbs.n_burn = 40;      % the conditional mean does not depend on
cfg0.gibbs.n_keep = 60;      % the chain length -- keep it short
cfg0.blp.return_draws = false;

dgps = {'correct', 'sparse', 'intermediate', 'dense'};
sim  = {@simulate_var_dgp, @simulate_sparse_misspec_dgp, ...
        @simulate_intermediate_misspec_dgp, @simulate_dense_misspec_dgp};
seeds = [4001 4002];
tol = 1e-10;

worst = 0;  worst_where = '';
n_checked = 0;

for h1 = {'bvar', 'lp'}
    h1_mode = h1{1};
    cfg = cfg0;
    cfg.fmar.h1_mode = h1_mode;
    if strcmp(h1_mode, 'lp'), h_lo = 1; else, h_lo = 2; end

    for d = 1:numel(dgps)
        for sd = seeds
            rng(sd, 'twister');
            dgp = sim{d}(cfg);
            Y = dgp.Y;

            bvar  = estimate_bvar_niw(Y, cfg);
            blp_f = estimate_blp_fmar(Y, cfg, bvar);

            cfg_fix = cfg;  cfg_fix.blp.fix_tau = 1;
            rng(sd + 500, 'twister');
            b_ind = estimate_blp_blockadaptive(Y, cfg_fix, bvar, blp_f.lambda);
            rng(sd + 900, 'twister');
            b_pool = estimate_blp_blockpooled(Y, cfg_fix, bvar, blp_f.lambda);

            cols = h_lo + 1 : cfg.H + 1;         % theta columns for h_lo..H
            for which = 1:2
                if which == 1
                    est = b_ind;  nm = 'independent';
                else
                    est = b_pool; nm = 'pooled';
                end
                E = abs(est.theta_cond(:, cols) - blp_f.theta_mean(:, cols));
                % per equation AND per horizon, not an aggregate
                [me, lin] = max(E(:));
                n_checked = n_checked + numel(E);
                if me > worst
                    worst = me;
                    [ii, hh] = ind2sub(size(E), lin);
                    worst_where = sprintf('%s/%s seed %d h1=%s eq %d h %d', ...
                        dgps{d}, nm, sd, h1_mode, ii, hh + h_lo - 1);
                end
                assert(me < tol, ...
                    ['nesting violated (%s, %s, seed %d, h1_mode=%s): ' ...
                     'max |theta_cond(tau=1) - FMAR| = %.3e over h >= %d'], ...
                    dgps{d}, nm, sd, h1_mode, me, h_lo);
            end

            pc = max(max(abs(b_ind.prior_theta - blp_f.prior_theta)));
            assert(pc < 1e-12, 'prior centres differ (%s, seed %d): %.3g', ...
                   dgps{d}, sd, pc);
            pc2 = max(max(abs(b_pool.prior_theta - blp_f.prior_theta)));
            assert(pc2 < 1e-12, 'pooled prior centres differ (%s): %.3g', ...
                   dgps{d}, pc2);
        end
    end
    fprintf('  h1_mode = %-4s : exact nesting verified for h = %d..%d\n', ...
            h1_mode, h_lo, cfg.H);
end

% --- the switch must not be inert --------------------------------------
cfg = cfg0;
cfg.gibbs.n_burn = 200;  cfg.gibbs.n_keep = 400;
rng(4321, 'twister');
dgp = simulate_sparse_misspec_dgp(cfg);
bvar  = estimate_bvar_niw(dgp.Y, cfg);
blp_f = estimate_blp_fmar(dgp.Y, cfg, bvar);
rng(4322, 'twister');
free_i = estimate_blp_blockadaptive(dgp.Y, cfg, bvar, blp_f.lambda);
rng(4323, 'twister');
free_p = estimate_blp_blockpooled(dgp.Y, cfg, bvar, blp_f.lambda);
dev_i = max(max(abs(free_i.theta_mean(:, 3:end) - blp_f.theta_mean(:, 3:end))));
dev_p = max(max(abs(free_p.theta_mean(:, 3:end) - blp_f.theta_mean(:, 3:end))));
assert(dev_i > 1e-4, 'independent tau sampling appears inert (%.2g)', dev_i);
assert(dev_p > 1e-4, 'pooled tau sampling appears inert (%.2g)', dev_p);

fprintf(['  %d (equation, horizon) comparisons; worst deviation %.3e ' ...
         '(%s), tol %.0e\n'], n_checked, worst, worst_where, tol);
fprintf('  free-tau deviation from the global solution: %.3f (independent), %.3f (pooled)\n', ...
        dev_i, dev_p);
fprintf('PASS: test_fmar_nesting_exact\n\n');
end
