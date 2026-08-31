function test_var_prior_mapping()
% PURPOSE
% -------
% Verify numerically that priors/var_implied_lp_prior.m maps the
% estimated VAR into the LP coefficients it truly implies.
%
% MODEL / EQUATIONS
% -----------------
% Three independent checks (tolerances at machine-precision scale):
%
% CHECK 1 (forecast equivalence): for every t in the sample and several
%   horizons h, the "LP prediction" z(t)' * Mu_h must equal the h-step
%   iterated forecast of the estimated VAR started from the state
%   (y(t), ..., y(t-p+1)) -- computed here by an EXPLICIT loop over the
%   VAR recursion, independent of any companion-matrix algebra.
%
% CHECK 2 (IRF block): the coefficient block on y(t) in Mu_h must equal
%   the reduced-form IRF matrix Psi_h -- recomputed here by the MA
%   recursion Psi_h = sum_j A_j Psi_{h-j}, again independent of the
%   companion-power code inside var_implied_lp_prior.
%
% CHECK 3 (h = 1 exactness): Mu_1 must equal the OLS VAR coefficient
%   matrix B itself, because the h = 1 "LP" IS the VAR regression.
%
% INPUTS / OUTPUTS
% ----------------
% None; prints PASS or raises an assertion error.
%
% DIMENSIONS
% ----------
% Uses default configuration (K = 3, p = 2, T = 200).
%
% NOTES
% -----
% Uses rng(cfg.seed) locally so the test is reproducible.

fprintf('--- test_var_prior_mapping ---\n');
cfg = default_config();
rng(cfg.seed, 'twister');

dgp = simulate_var_dgp(cfg);
var_est = estimate_var(dgp.Y, cfg);
Y = dgp.Y;
[T, K] = size(Y);
p = cfg.p;
Zall = build_lp_regressors(Y, p);

% ---- CHECK 3: h = 1 equals the OLS coefficients ----------------------
Mu1 = var_implied_lp_prior(var_est, 1);
err3 = max(abs(Mu1(:) - var_est.B(:)));
assert(err3 < 1e-10, 'CHECK 3 failed: |Mu_1 - B| = %.2e', err3);
fprintf('  CHECK 3 (h=1 equals OLS B): max err %.2e  OK\n', err3);

% ---- CHECK 1: iterated forecasts ---------------------------------------
for h = [2, 5, 10, cfg.H]
    Mu = var_implied_lp_prior(var_est, h);
    pred_lp = Zall * Mu;                        % ((T-p+1) x K)

    % Explicit iterated VAR forecast from each state, no companion algebra.
    n_states = size(Zall, 1);
    pred_iter = zeros(n_states, K);
    for r = 1:n_states
        t = p + r - 1;
        hist = Y(t - p + 1:t, :);               % (p x K), oldest first
        for step = 1:h
            f = var_est.c;
            for j = 1:p
                f = f + var_est.A(:, :, j) * hist(end - j + 1, :)';
            end
            hist = [hist; f'];                  % append forecast
        end
        pred_iter(r, :) = hist(end, :);
    end
    err1 = max(abs(pred_lp(:) - pred_iter(:)));
    assert(err1 < 1e-8, 'CHECK 1 failed at h=%d: err %.2e', h, err1);
    fprintf('  CHECK 1 (iterated forecast, h=%2d): max err %.2e  OK\n', h, err1);
end

% ---- CHECK 2: y(t) block equals Psi_h (independent MA recursion) ------
Psi = zeros(K, K, cfg.H + 1);
Psi(:, :, 1) = eye(K);
for h = 1:cfg.H
    acc = zeros(K);
    for j = 1:min(h, p)
        acc = acc + var_est.A(:, :, j) * Psi(:, :, h - j + 1);
    end
    Psi(:, :, h + 1) = acc;
end
for h = [1, 3, 8, cfg.H]
    Mu = var_implied_lp_prior(var_est, h);
    block_on_yt = Mu(2:1 + K, :)';              % (K x K)
    err2 = max(max(abs(block_on_yt - Psi(:, :, h + 1))));
    assert(err2 < 1e-10, 'CHECK 2 failed at h=%d: err %.2e', h, err2);
end
fprintf('  CHECK 2 (y(t) block equals Psi_h): OK\n');
fprintf('PASS: test_var_prior_mapping\n\n');
end
