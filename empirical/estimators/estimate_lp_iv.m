function lp = estimate_lp_iv(Y, z, cfg, opts)
% ESTIMATE_LP_IV  Local projection with an external instrument (LP-IV,
% Stock & Watson, Economic Journal 2018), with HAC Wald bands, a
% lag-augmented EHW variant (Montiel Olea & Plagborg-Moller 2021) and an
% Anderson-Rubin confidence set that is valid under weak identification
% (Montiel Olea, Stock & Watson 2021).
%
% PURPOSE
% -------
% Frequentist benchmark for the proxy-identified Chapter 7 results.  The
% Bayesian estimators in the repo hold the proxy impact vector b_z fixed
% (Gamma_h * b_z, docs/APPROXIMATIONS.md); this estimator instead carries
% the instrument uncertainty in full, horizon by horizon, and reports how
% much the bands widen once weak-instrument-robust inference is used.
%
% MODEL / EQUATIONS
% -----------------
% For h = 0..H and each equation i (s = cfg.shock_var, the policy
% indicator), the horizon-h structural equation is
%     y_{i,t+h} = theta_i(h) * y_{s,t} + gamma' w_t + u_{i,t+h},
%     w_t = [1; y_{t-1}; ...; y_{t-p}]                 (p lags of ALL variables),
% with y_{s,t} instrumented by z_t.  The controls deliberately EXCLUDE the
% contemporaneous y_t of the other variables: the instrument, not a
% recursive ordering, does the identification.  theta_s(0) = 1 exactly
% (unit-effect normalisation), so theta_i(h) is the response to a shock
% that raises the indicator by one unit on impact.
%
% 2SLS BY PARTIALLING (Frisch-Waugh-Lovell).  Let a = M_w y_{i,t+h},
% b = M_w y_{s,t}, zt = M_w z_t be the residuals from OLS on w_t.  Then
%     theta_hat = (zt' a) / (zt' b),
% which is EXACTLY the first element of the full just-identified 2SLS
% estimator with instruments [z_t; w_t] and regressors [y_{s,t}; w_t]:
% e_1'(Q'X)^{-1} Q' u = (zt' b)^{-1} zt' u by the partitioned inverse,
% so the sandwich variance of theta_hat is
%     V = (zt' b)^{-2} * S,   S = long-run variance of g_t = zt_t * uhat_t,
% with uhat = a - theta_hat b the 2SLS residual (again exact: FWL).  S is
% Newey-West with Bartlett weights w_l = 1 - l/(L+1), L = h + 1 (same
% truncation as the FMAR sandwich in estimate_blp_fmar.m).
%
% LAG-AUGMENTED VARIANT (opts.lag_augment).  Same 2SLS with w_t carrying
% p + 1 lags and EHW (L = 0) errors.  MOP-M 2021 show the extra lag makes
% the LP score serially uncorrelated (also for LP-IV, their Section 5), so
% EHW inference is valid without a truncation choice and uniformly over
% persistence.  Returned as .theta_la / .se_la / .lo_la / .hi_la.
%
% ANDERSON-RUBIN SET.  For a null theta0, regress r(theta0) = a - theta0 b
% on zt and test the slope c(theta0) = 0 with the same HAC.  Algebra that
% makes the grid O(n_grid) after an O(N L) setup per (i, h):
%     c(theta0)   = (Sza - theta0 Szb) / Szz,     Sxy = zt' x,
%     e_t(theta0) = a_t - theta0 b_t - c(theta0) zt_t,
%     g_t(theta0) = zt_t e_t = A_t - theta0 B_t,
%         A_t = zt_t (a_t - (Sza/Szz) zt_t),  B_t = zt_t (b_t - (Szb/Szz) zt_t),
%     LRV(g)      = Q_AA - 2 theta0 Q_AB + theta0^2 Q_BB,
%         Q_XY = X' Omega Y,  Omega_{tr} = w_{|t-r|} (Bartlett, symmetric),
%     t(theta0)   = (Sza - theta0 Szb) / sqrt(LRV(g)).
% The set is {theta0 : |t(theta0)| <= zcrit} on a grid of opts.ar_n_grid
% points spanning theta_hat +/- opts.ar_grid_halfwidth * se.  Because the
% grid has an odd number of points its centre is theta_hat itself, where
% t = 0, so the 2SLS point always belongs to the set.  ar_unbounded flags
% that a grid endpoint is accepted (the set runs off the grid: with a
% single instrument |t(theta0)| tends to the first-stage HAC t as
% |theta0| -> inf, so an insignificant first stage means an unbounded
% set); ar_empty flags that nothing is accepted (only possible in the
% degenerate se = 0 cell, handled separately).
%
% FIRST STAGE (h = 0 sample): b = pi zt + v.  Reports pi, HAC t (L = 1),
% EHW t, the effective F = HAC t^2 (with one instrument the Montiel
% Olea-Pflueger 2013 effective F is the robust F; critical values
% 12.05 / 15.06 / 23.11 / 37.42 for 30 / 20 / 10 / 5 % Nagar bias),
% partial R^2 = Szb^2 / (Szz Sbb), N and the number of nonzero z_t.
%
% INPUTS
% ------
%   Y    : (T x K) raw data.
%   z    : (T x 1) instrument aligned with the rows of Y.  NaN allowed:
%          a row with NaN z is dropped from EVERY regression it would
%          enter (the controls and leads for the other rows are unaffected).
%   cfg  : uses cfg.p, cfg.H, cfg.ci_level, cfg.shock_var (default 1).
%   opts : optional struct with fields (defaults in brackets)
%          .lag_augment        [true]  compute the LA/EHW variant
%          .ar_grid_halfwidth  [10]    AR grid half-width in HAC se units
%          .ar_n_grid          [401]   AR grid points (made odd if even)
%
% OUTPUTS
% -------
%   lp.theta, .se, .lo, .hi        (K x (H+1)) 2SLS point, HAC se, Wald band
%   lp.theta_mean                  = lp.theta (interface parity with the repo)
%   lp.ar_lo, .ar_hi               (K x (H+1)) AR set endpoints on the grid
%   lp.ar_unbounded, .ar_empty     (K x (H+1)) logical flags
%   lp.theta_la, .se_la, .lo_la, .hi_la   lag-augmented EHW variant ([] if off)
%   lp.first_stage                 struct: .pi, .t_hac, .t_ehw, .F_eff,
%                                  .partial_R2, .N, .n_nonzero_z, .crit_F
%   lp.N_h                         (1 x (H+1)) usable observations per h
%   lp.hac_L                       (1 x (H+1)) Bartlett truncation per h
%   lp.name = 'LP-IV', plus .p, .H, .shock_var, .ci_level, .opts
%
% NOTES
% -----
% * Regressor layout: build_lp_regressors(Y, p) row r holds
%   [1, y_t, ..., y_{t-p+1}] for t = p + r - 1.  Using rows 1..end-1 as
%   the controls of the NEXT period gives w_t = [1; y_{t-1}; ...; y_{t-p}]
%   for t = p+1..T without a second design matrix.
% * Sample at horizon h: t = p+1..T-h intersected with {z_t finite}, so
%   N_h = #{t in p+1..T-h : z_t finite}.  The LA variant starts at p+2.
% * No degrees-of-freedom correction anywhere (asymptotic inference).
% * The (s, h = 0) cell is the identity theta = 1, se = 0: the AR set is
%   reported as the point {1} and no grid is built.

% ---------------------------------------------------------------- setup
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'lag_augment'),       opts.lag_augment       = true; end
if ~isfield(opts, 'ar_grid_halfwidth'), opts.ar_grid_halfwidth = 10;   end
if ~isfield(opts, 'ar_n_grid'),         opts.ar_n_grid         = 401;  end
if mod(opts.ar_n_grid, 2) == 0, opts.ar_n_grid = opts.ar_n_grid + 1; end

p = cfg.p;  H = cfg.H;
s = 1; if isfield(cfg, 'shock_var') && ~isempty(cfg.shock_var), s = cfg.shock_var; end
[T, K] = size(Y);
z = z(:);
assert(numel(z) == T, 'estimate_lp_iv: z must have one entry per row of Y.');
assert(s >= 1 && s <= K, 'estimate_lp_iv: cfg.shock_var out of range.');
zc = normal_quantile(1 - (1 - cfg.ci_level) / 2);

% Controls w_t = [1; y_{t-1}; ...; y_{t-p}] for t = p+1..T (see NOTES).
Zall  = build_lp_regressors(Y, p);
Wall  = Zall(1:end - 1, :);                 % rows t = p+1..T
t0    = p + 1;                              % first usable t
N0    = T - t0 + 1;
ys_all = Y(t0:T, s);
z_all  = z(t0:T);
zok    = isfinite(z_all);

lp.theta = zeros(K, H + 1);   lp.se = zeros(K, H + 1);
lp.ar_lo = zeros(K, H + 1);   lp.ar_hi = zeros(K, H + 1);
lp.ar_unbounded = false(K, H + 1);  lp.ar_empty = false(K, H + 1);
lp.N_h   = zeros(1, H + 1);   lp.hac_L = (0:H) + 1;
first_stage = [];

for h = 0:H
    L    = h + 1;
    keep = find(zok(1:N0 - h));             % rows of Wall usable at h
    n    = numel(keep);
    assert(n > size(Wall, 2) + 5, 'estimate_lp_iv: sample too short at h = %d.', h);
    lp.N_h(h + 1) = n;

    W   = Wall(keep, :);
    Yh  = Y(t0 + h + keep - 1, :);          % y_{t+h}, all K equations
    R   = [ys_all(keep), z_all(keep), Yh];
    R   = R - W * (W \ R);                  % partial out w_t (FWL)
    b   = R(:, 1);  zt = R(:, 2);  A = R(:, 3:end);

    Szz = zt' * zt;  Szb = zt' * b;
    if h == 0
        first_stage = first_stage_stats(zt, b, Szz, Szb, z_all(keep));
    end

    for i = 1:K
        a   = A(:, i);
        Sza = zt' * a;
        th  = Sza / Szb;
        if h == 0 && i == s
            lp.theta(i, 1) = 1;  lp.se(i, 1) = 0;   % identity cell
            lp.ar_lo(i, 1) = 1;  lp.ar_hi(i, 1) = 1;
            continue
        end
        % ---- AR / HAC setup: A_t, B_t and the three kernel cross-products
        At  = zt .* (a - (Sza / Szz) * zt);
        Bt  = zt .* (b - (Szb / Szz) * zt);
        Qaa = bartlett_xprod(At, At, L);
        Qab = bartlett_xprod(At, Bt, L);
        Qbb = bartlett_xprod(Bt, Bt, L);
        % HAC se of theta_hat: LRV(zt .* uhat) = Q(theta_hat) since
        % c(theta_hat) = 0 makes g_t(theta_hat) = zt_t * uhat_t exactly.
        S   = max(Qaa - 2 * th * Qab + th^2 * Qbb, 0);
        se  = sqrt(S) / abs(Szb);
        lp.theta(i, h + 1) = th;
        lp.se(i, h + 1)    = se;
        % ---- Anderson-Rubin set on the grid
        if se > 0
            grid = th + linspace(-1, 1, opts.ar_n_grid) * opts.ar_grid_halfwidth * se;
            num  = Sza - grid * Szb;
            lrv  = max(Qaa - 2 * grid * Qab + grid.^2 * Qbb, 0);
            tst  = num ./ sqrt(lrv);
            tst(lrv == 0 & num == 0) = 0;   % 0/0 only at theta_hat in degenerate cases
            acc  = abs(tst) <= zc;
            if any(acc)
                lp.ar_lo(i, h + 1) = min(grid(acc));
                lp.ar_hi(i, h + 1) = max(grid(acc));
                lp.ar_unbounded(i, h + 1) = acc(1) || acc(end);
            else
                lp.ar_lo(i, h + 1) = NaN;  lp.ar_hi(i, h + 1) = NaN;
                lp.ar_empty(i, h + 1) = true;
            end
        else
            lp.ar_lo(i, h + 1) = th;  lp.ar_hi(i, h + 1) = th;
        end
    end
end

lp.theta_mean = lp.theta;
lp.lo = lp.theta - zc * lp.se;
lp.hi = lp.theta + zc * lp.se;

% ------------------------------------------- lag-augmented EHW variant
lp.theta_la = [];  lp.se_la = [];  lp.lo_la = [];  lp.hi_la = [];
if opts.lag_augment
    pa     = p + 1;
    Wla    = build_lp_regressors(Y, pa);  Wla = Wla(1:end - 1, :);   % t = pa+1..T
    ta     = pa + 1;  Na = T - ta + 1;
    ys_la  = Y(ta:T, s);  z_la = z(ta:T);  zok_la = isfinite(z_la);
    lp.theta_la = zeros(K, H + 1);  lp.se_la = zeros(K, H + 1);
    for h = 0:H
        keep = find(zok_la(1:Na - h));
        W  = Wla(keep, :);
        R  = [ys_la(keep), z_la(keep), Y(ta + h + keep - 1, :)];
        R  = R - W * (W \ R);
        b  = R(:, 1);  zt = R(:, 2);
        Szb = zt' * b;
        for i = 1:K
            a  = R(:, 2 + i);
            th = (zt' * a) / Szb;
            g  = zt .* (a - th * b);
            lp.theta_la(i, h + 1) = th;
            lp.se_la(i, h + 1)    = sqrt(g' * g) / abs(Szb);       % EHW, L = 0
        end
    end
    lp.theta_la(s, 1) = 1;  lp.se_la(s, 1) = 0;
    lp.lo_la = lp.theta_la - zc * lp.se_la;
    lp.hi_la = lp.theta_la + zc * lp.se_la;
end

lp.first_stage = first_stage;
lp.name = 'LP-IV';
lp.p = p;  lp.H = H;  lp.shock_var = s;  lp.ci_level = cfg.ci_level;
lp.opts = opts;
end

% =====================================================================
function Q = bartlett_xprod(x, y, L)
% Q = x' Omega y with Omega_{tr} = max(1 - |t-r|/(L+1), 0): the Newey-West
% (Bartlett) weighted sum of cross-products at lags -L..L.  O(N L).
Q = x' * y;
for l = 1:L
    w = 1 - l / (L + 1);
    Q = Q + w * (x(1 + l:end)' * y(1:end - l) + x(1:end - l)' * y(1 + l:end));
end
end

% =====================================================================
function fs = first_stage_stats(zt, b, Szz, Szb, zraw)
% First stage b = pi zt + v on the partialled h = 0 sample.
pi_hat = Szb / Szz;
v  = b - pi_hat * zt;
g  = zt .* v;
S_hac = bartlett_xprod(g, g, 1);
S_ehw = g' * g;
fs.pi          = pi_hat;
fs.se_hac      = sqrt(S_hac) / Szz;
fs.se_ehw      = sqrt(S_ehw) / Szz;
fs.t_hac       = pi_hat / fs.se_hac;
fs.t_ehw       = pi_hat / fs.se_ehw;
fs.F_eff       = fs.t_hac^2;
fs.F_ehw       = fs.t_ehw^2;
fs.partial_R2  = Szb^2 / (Szz * (b' * b));
fs.N           = numel(zt);
fs.n_nonzero_z = sum(zraw ~= 0);
fs.crit_F      = struct('nagar_bias_30', 12.05, 'nagar_bias_20', 15.06, ...
                        'nagar_bias_10', 23.11, 'nagar_bias_05', 37.42, ...
                        'rule_of_thumb', 10);
fs.hac_L       = 1;
end
