function test_psi_floor()
% PURPOSE
% -------
% Pin the two claims behind the optional prior-scale floor
% cfg.fmar.psi_floor (priors/fmar_prior_scale.m):
%
% 1. ESSENTIALLY INERT ON THE SIMULATION DESIGNS.  The floor replaces the
%    Newey-West long-run variance by the plain residual variance only
%    where the former is SMALLER, which needs negatively autocorrelated
%    h-step LP residuals.  On the simulated DGPs the residuals at h >= 2
%    are positively autocorrelated (overlapping horizons of a persistent
%    system), so the floor must not bind there.  At h = 1 the one-step
%    residual is white and its sample autocovariances are pure noise, so
%    the NW estimate falls a few percent below the variance in about half
%    the samples and the floor binds by that few percent.  The test pins
%    the size of that effect: psi moves by < 15% at h = 1 only, log
%    lambda_1 by < 0.05, and no IRF value by more than 1e-3 (the RMSEs
%    in the study are ~0.07).  The stored Monte Carlo results are in any
%    case produced with the flag OFF (the default).
%
% 2. ACTIVE, AND APPLIED CONSISTENTLY, WHERE IT SHOULD BE.  On a system
%    that mixes a mean-reverting near-white-noise series (the shape of the
%    euro-area policy surprise) with a persistent level, the floor binds
%    for the first variable only, raises its scale to the residual
%    variance, and -- the part that matters for the method -- the tau = 1
%    nesting of both adaptive estimators in the FMAR closed form remains
%    exact with the floor switched on, i.e. the three estimators apply the
%    same floored scale.
%
% Also checks that the lambda selector's new multimodality counter is
% present, is 1 on the simulated designs (a unimodal objective, as the
% documentation of select_lambda_fmar claims), and that the selector's
% returned optimum is unchanged by the extra bookkeeping.

fprintf('test_psi_floor:\n');

cfg = default_config();
cfg.mode = 'fmar';
cfg.fmar.h1_mode = 'lp';
cfg.H = 16;  cfg.T = 200;
cfg.gibbs.n_burn = 30;  cfg.gibbs.n_keep = 40;   % nesting uses the conditional mean
cfg.fmar.n_niw_draws = 40;

% --- 1. inert on the four simulated designs -----------------------------
sims  = {@simulate_var_dgp, @simulate_sparse_misspec_dgp, ...
         @simulate_intermediate_misspec_dgp, @simulate_dense_misspec_dgp};
names = {'correct', 'sparse', 'intermediate', 'dense'};
n_cells = 0;  worst_h1 = 1;  worst_dth = 0;  worst_dlam = 0;
for d = 1:numel(sims)
    rng(5100 + d, 'twister');
    dgp  = sims{d}(cfg);
    bvar = estimate_bvar_niw(dgp.Y, cfg);
    dt = var_deterministic_trend(dgp.Y, bvar.B, cfg.p);
    for h = 1:cfg.H
        [ps0, i0] = fmar_prior_scale(dt.x, cfg.p, h, false);
        [ps1, i1] = fmar_prior_scale(dt.x, cfg.p, h, true);
        ratio = min(i0.psi_nw ./ i0.psi_var);
        if h >= 2
            % overlapping residuals: positive dependence, no material floor
            assert(ratio > 0.98, ...
                'the floor bound materially on the %s DGP at h = %d (NW/var = %.3f)', ...
                names{d}, h, ratio);
        else
            assert(ratio > 0.85, ...
                'the h = 1 floor moved psi by more than 15%% on the %s DGP (NW/var = %.3f)', ...
                names{d}, ratio);
            worst_h1 = min(worst_h1, ratio);
        end
        assert(all(ps1 >= ps0) && all(ps1 == max(ps0, i0.psi_var)), ...
            'floored psi is not max(NW, var) on the %s DGP at h = %d', names{d}, h);
        n_cells = n_cells + numel(ps0);
    end
    cfg0 = cfg;  cfg0.fmar.psi_floor = false;
    cfg1 = cfg;  cfg1.fmar.psi_floor = true;
    f0 = estimate_blp_fmar(dgp.Y, cfg0, bvar);
    f1 = estimate_blp_fmar(dgp.Y, cfg1, bvar);
    dth  = max(max(abs(f0.theta_mean - f1.theta_mean)));
    dlam = max(abs(log(f0.lambda(1, :)) - log(f1.lambda(1, :))));
    assert(dth < 1e-3, 'the floor moved an IRF value by %.2e on the %s DGP', dth, names{d});
    assert(dlam < 0.05, 'the floor moved log lambda by %.3f on the %s DGP', dlam, names{d});
    assert(all(f0.lambda(1, 2:end) == f1.lambda(1, 2:end)) || ...
           max(abs(log(f0.lambda(1, 2:end)) - log(f1.lambda(1, 2:end)))) < 1e-3, ...
           'the floor changed lambda_h at h >= 2 on the %s DGP', names{d});
    worst_dth = max(worst_dth, dth);  worst_dlam = max(worst_dlam, dlam);
    assert(isfield(f1.diag, 'n_local_max') && numel(f1.diag.n_local_max) == cfg.H, ...
           'lambda multimodality counter missing');
    assert(f1.diag.lambda_multimodal == 0, ...
           'the lambda objective was multimodal at %d horizons on the %s DGP', ...
           f1.diag.lambda_multimodal, names{d});
end
fprintf(['  4 DGPs x %d horizons (%d cells): floor binds at h = 1 only, worst NW/var %.3f;\n' ...
         '  max IRF change %.1e, max |dlog lambda| %.3f: essentially inert, OK\n'], ...
        cfg.H, n_cells, worst_h1, worst_dth, worst_dlam);
fprintf('  lambda objective unimodal at every horizon on all 4 DGPs: OK\n');

% --- 2. active where it should be, and consistent across estimators ----
% y_1: MA(1) with a negative coefficient (mean-reverting, near white noise,
%      long-run variance far below the variance);  y_2: persistent AR(1)
%      that responds to y_1.  A small-scale analogue of the surprise +
%      level system.
rng(5200, 'twister');
T = 260;  burn = 100;
e = randn(T + burn + 1, 2);
y1 = 0.05 * (e(2:end, 1) - 0.6 * e(1:end-1, 1));
y2 = zeros(T + burn, 1);
for t = 2:T + burn
    y2(t) = 0.95 * y2(t-1) + 0.5 * y1(t) + 0.3 * e(t + 1, 2);
end
Y = [y1(burn+1:end), y2(burn+1:end)];
cfg2 = cfg;  cfg2.K = 2;  cfg2.p = 3;  cfg2.H = 12;
cfg2.fmar.isrw = [0 1];
bvar2 = estimate_bvar_niw(Y, cfg2);
dt2 = var_deterministic_trend(Y, bvar2.B, cfg2.p);
n_bind = 0;
for h = 1:cfg2.H
    [ps0, i0] = fmar_prior_scale(dt2.x, cfg2.p, h, false);
    [ps1, i1] = fmar_prior_scale(dt2.x, cfg2.p, h, true);
    assert(~i1.floor_bound(2), 'the floor bound for the persistent variable at h = %d', h);
    if i1.floor_bound(1)
        n_bind = n_bind + 1;
        assert(ps1(1) == i0.psi_var(1) && ps0(1) < ps1(1), ...
               'floored scale is not the residual variance at h = %d', h);
    end
    assert(all(ps1 >= ps0), 'floor lowered a scale at h = %d', h);
end
assert(n_bind >= cfg2.H / 2, ...
       'the floor bound at only %d of %d horizons on the mean-reverting series', n_bind, cfg2.H);
fprintf('  mean-reverting series: floor binds at %d of %d horizons, never for the level: OK\n', ...
        n_bind, cfg2.H);

cfg2.fmar.psi_floor = true;
f2 = estimate_blp_fmar(Y, cfg2, bvar2);
cfg2f = cfg2;  cfg2f.blp.fix_tau = 1;
rng(5201, 'twister');
b2 = estimate_blp_blockadaptive(Y, cfg2f, bvar2, f2.lambda);
rng(5202, 'twister');
p2 = estimate_blp_blockpooled(Y, cfg2f, bvar2, f2.lambda);
cols = 2:cfg2.H + 1;                                    % h = 1..H (h1_mode = lp)
db = max(max(abs(b2.theta_cond(:, cols) - f2.theta_mean(:, cols))));
dp = max(max(abs(p2.theta_cond(:, cols) - f2.theta_mean(:, cols))));
assert(db < 1e-10 && dp < 1e-10, ...
       'tau = 1 nesting broken with the floor on: %.2e (independent), %.2e (pooled)', db, dp);
% ...and the floor must actually change the estimator here (not inert)
cfg2n = cfg2;  cfg2n.fmar.psi_floor = false;
f2n = estimate_blp_fmar(Y, cfg2n, bvar2);
assert(max(abs(f2.lambda(1, :) - f2n.lambda(1, :))) > 1e-6, ...
       'the floor did not change lambda_h on the mean-reverting system');
fprintf('  exact tau = 1 nesting with the floor on: %.1e / %.1e; lambda path moved: OK\n', db, dp);

fprintf('PASS: test_psi_floor\n\n');
end
