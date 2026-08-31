function lp = estimate_lp_lagaug(Y, cfg, var_est)
% ESTIMATE_LP_LAGAUG  Lag-augmented local projection (Montiel Olea &
% Plagborg-Moller, Econometrica 2021, "Local projection inference is
% simpler and more robust than you think").
%
% PURPOSE
% -------
% Frequentist benchmark aimed at an LP-oriented audience.  Two differences
% from the plain LP in the repo: (i) LAG AUGMENTATION -- the horizon
% regressions include p + 1 lag blocks instead of p, and (ii) plain
% Eicker-Huber-White standard errors instead of HAC.  MOP-M show that with
% the extra lag the score is serially uncorrelated, so EHW inference is
% asymptotically valid WITHOUT Newey-West truncation choices, uniformly
% over stationary and (near-)unit-root persistence -- attractive here since
% the euro-area levels are highly persistent.
%
% MODEL / EQUATIONS
% -----------------
% For each equation i and horizon h >= 1, OLS on raw (non-detrended) data:
%     y_{i,t+h} = z_t' beta_{i,h} + u_{i,t+h},
%     z_t = [1; y_t; y_{t-1}; ...; y_{t-p}]          ((p+1) lag blocks).
% The K x K matrix Gamma_h formed by stacking each equation's coefficient
% block on y_t consistently estimates the reduced-form Wold coefficient
% Psi_h; structural IRFs are theta(:, h+1) = Gamma_h * b1n, with b1n the
% recursive impact vector held FIXED at var_est.b1n (the same convention
% as the FMAR-mode bands in the repo: impact-normalisation uncertainty is
% not propagated, keeping bands comparable across estimators).
%     se(theta_i(h)) = sqrt(b1n' * V_i * b1n),
% where V_i is the y_t-block of the EHW covariance for equation i.
% h = 0: theta = b1n, se = 0 (impact identity).  Bands: theta +/- z * se
% at cfg.ci_level.
%
% INPUTS
% ------
%   Y       : (T x K) raw data.
%   cfg     : uses cfg.p, cfg.H, cfg.ci_level.
%   var_est : estimate_bvar_niw (or estimate_var) output; only .b1n used.
%
% OUTPUT
% ------
%   lp.theta / .theta_mean (K x (H+1)), lp.se, lp.lo, lp.hi, lp.name.
%
% NOTES
% -----
% * Uses build_lp_regressors(Y, p+1); its layout z = [1; y(t); y(t-1); ...]
%   puts the y_t block at positions 2..K+1 (verified in tests of the port).
% * Effective sample at horizon h: T - (p+1) - h + 1 + ... (one obs less at
%   the start than the p-lag LP; the honest price of augmentation).

p = cfg.p;  H = cfg.H;
[T, K] = size(Y);
alpha = 1 - cfg.ci_level;
zc = normal_quantile(1 - alpha / 2);
b1n = var_est.b1n;

pa = p + 1;                                   % augmented lag order
Zall = build_lp_regressors(Y, pa);            % rows t = pa..T-?; [1, y_t, ...]
x1 = Y(pa:T, :);                              % dependent base: y_t aligned to Zall
% Zall has one row per t = pa..T with z_t = [1, y_t', ..., y_{t-pa+1}'];
% regressing y_{t+h} on z_t uses rows 1..end-h of Zall and y rows pa+h..T.

m = 1 + K * pa;
assert(size(Zall, 2) == m, 'estimate_lp_lagaug: unexpected regressor layout.');

lp.theta = zeros(K, H + 1);
lp.se    = zeros(K, H + 1);
lp.theta(:, 1) = b1n;                          % h = 0

idx = 2:K + 1;                                 % y_t block
for h = 1:H
    Zh = Zall(1:end - h, :);
    Yh = x1(1 + h:end, :);
    n  = size(Zh, 1);
    assert(n > m + 5, 'estimate_lp_lagaug: sample too short at h = %d.', h);
    ZtZi = inv(Zh' * Zh);                      %#ok<MINV>
    Gam  = zeros(K, K);
    for i = 1:K
        bi = ZtZi * (Zh' * Yh(:, i));
        ui = Yh(:, i) - Zh * bi;
        meat = Zh' * (Zh .* (ui.^2 * ones(1, m)));
        Vi   = ZtZi * meat * ZtZi;             % EHW sandwich
        Gam(i, :) = bi(idx)';
        lp.se(i, h + 1) = sqrt(max(b1n' * Vi(idx, idx) * b1n, 0));
    end
    lp.theta(:, h + 1) = Gam * b1n;
end

lp.theta_mean = lp.theta;
lp.lo = lp.theta - zc * lp.se;
lp.hi = lp.theta + zc * lp.se;
lp.name = 'LP-LA (MOP)';
lp.p_augmented = pa;
end
