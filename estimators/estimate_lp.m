function lp_est = estimate_lp(Y, cfg, b1n)
% PURPOSE
% -------
% Ordinary (frequentist) local projections, estimated by OLS horizon by
% horizon.  This is the unregularised benchmark: unbiased-ish but noisy
% at long horizons.
%
% MODEL / EQUATIONS
% -----------------
% For each horizon h = 1,...,H run the SYSTEM regression
%
%     y(t+h)' = z(t)' * Beta_h + u(t+h)',    z(t) = [1; y(t); ...; y(t-p+1)],
%
% with Beta_h (m x K), m = 1 + K*p, estimated by OLS (Z_h \ Y_h).
%
% Shock variable and controls:
%   * the "shock" regressor block is y(t) (columns 2 : 1+K of z);
%   * lag controls are y(t-1), ..., y(t-p+1) plus the intercept.
% Under recursive identification with the shocked variable ordered
% first, the structural response of variable i at horizon h is the
% linear combination of the coefficients on y(t) with the (fixed,
% VAR-estimated) impact vector b1n:
%
%     theta_i(h) = b1n' * Beta_h(2:1+K, i),      h >= 1,
%     theta(:, 1) = b1n                          (h = 0, impact).
%
% Inference: Var(theta_i(h)) = w' * Vcov(beta_i) * w with the fixed
% weight vector w = [0; b1n; 0; ...; 0] (uncertainty in b1n itself is
% IGNORED in this prototype -- documented simplification shared by the
% Bayesian estimators, so comparisons focus on the dynamic
% coefficients).  Vcov(beta_i) is computed in the ISOLATED local
% function lp_vcov() below:
%   'ols' : sigma_i^2 (Z'Z)^{-1}     (iid; wrong for h > 1 because
%           u(t+h) is MA(h-1) by construction)
%   'nw'  : Newey-West HAC, Bartlett kernel, truncation lag L = h.
% Replace lp_vcov() with a better HAC / small-sample procedure later
% without touching anything else.
%
% INPUTS
% ------
% Y   : (T x K) data.
% cfg : configuration struct (cfg.p, cfg.H, cfg.ci_level, cfg.lp.vcov).
% b1n : (K x 1) normalised impact vector from estimate_var (shared
%       identification device across all estimators).
%
% OUTPUTS
% -------
% lp_est : struct with fields
%   .theta (K x (H+1)) point estimates (column h+1 = horizon h)
%   .se    (K x (H+1)) standard errors (NaN at h = 0)
%   .lo/.hi (K x (H+1)) confidence bands at level cfg.ci_level
%   .Beta  {H x 1} cell of (m x K) coefficient matrices (debugging)
%
% DIMENSIONS
% ----------
% Horizon-h regression: N_h = T - p + 1 - h observations, m regressors.
%
% NOTES
% -----
% h = 0 is not estimated: it is the identification step itself and is
% fixed at b1n for every estimator (see README).

[T, K] = size(Y);
p = cfg.p;  H = cfg.H;
zcrit = normal_quantile(1 - (1 - cfg.ci_level) / 2);

Zall = build_lp_regressors(Y, p);       % ((T-p+1) x m)
m    = size(Zall, 2);
w    = [0; b1n; zeros(K * (p - 1), 1)]; % weights picking the y(t) block
assert(numel(w) == m, 'estimate_lp: weight vector dimension mismatch.');

theta = zeros(K, H + 1);  se = nan(K, H + 1);
theta(:, 1) = b1n;                      % h = 0: impact normalisation
Beta_all = cell(H, 1);

for h = 1:H
    Zh = Zall(1:end - h, :);            % (N_h x m), t = p .. T-h
    Yh = Y(p + h:T, :);                 % (N_h x K), y(t+h)
    Nh = size(Zh, 1);
    assert(Nh > m + 5, 'estimate_lp: too few observations at h = %d.', h);

    Bh = Zh \ Yh;                       % (m x K) OLS
    Uh = Yh - Zh * Bh;                  % (N_h x K) residuals
    Beta_all{h} = Bh;

    theta(:, h + 1) = Bh(2:1 + K, :)' * b1n;    % (K x 1)

    ZtZ = Zh' * Zh;
    for i = 1:K
        V = lp_vcov(Zh, Uh(:, i), ZtZ, h, cfg);  % (m x m)
        se(i, h + 1) = sqrt(max(w' * V * w, 0));
    end
end

lp_est.theta = theta;
lp_est.se    = se;
lp_est.lo    = theta - zcrit * se;
lp_est.hi    = theta + zcrit * se;
lp_est.Beta  = Beta_all;

assert(all(isfinite(theta(:))), 'estimate_lp: non-finite estimates.');
end

% =====================================================================
function V = lp_vcov(Z, u, ZtZ, h, cfg)
% ISOLATED covariance computation for one LP equation.  Swap this
% function to change LP inference project-wide.
%
% 'ols': V = sigma2 * (Z'Z)^{-1},  sigma2 = u'u / (N - m).
% 'nw' : sandwich V = (Z'Z)^{-1} * S * (Z'Z)^{-1} with the Newey-West
%        long-run score covariance
%        S = G_0 + sum_{l=1}^{L} w_l (G_l + G_l'),  w_l = 1 - l/(L+1),
%        G_l = sum_t g(t) g(t-l)',  g(t) = z(t) * u(t),  L = h.
% The truncation-lag rule L = h reflects that u(t+h) is MA(h-1) under
% correct specification; it is a simple prototype rule, replace later.
[N, m] = size(Z);
switch cfg.lp.vcov
    case 'ols'
        sigma2 = (u' * u) / (N - m);
        V = sigma2 * inv_from_chol(ZtZ);
    case 'nw'
        G = Z .* u;                     % (N x m), rows g(t)'
        S = G' * G;                     % G_0
        L = h;
        for l = 1:min(L, N - 1)
            wl = 1 - l / (L + 1);
            Gl = G(l + 1:end, :)' * G(1:end - l, :);
            S  = S + wl * (Gl + Gl');
        end
        ZtZinv = inv_from_chol(ZtZ);
        V = ZtZinv * S * ZtZinv;
    otherwise
        error('estimate_lp: unknown cfg.lp.vcov "%s".', cfg.lp.vcov);
end
end

function Ainv = inv_from_chol(A)
% Inverse of a symmetric positive definite matrix via Cholesky solves,
% avoiding inv() on a possibly ill-conditioned raw matrix.
R = chol((A + A') / 2);                 % A = R'R
Ainv = R \ (R' \ eye(size(A, 1)));
end
