function bz = ea_identify_proxy(bvar, Y, z, cfg, opts)
% EA_IDENTIFY_PROXY  External-instrument (proxy) identification of the
% impact vector, dropped into the BVAR struct so every estimator uses it.
%
% PURPOSE
% -------
% v2 identification (docs/CH7_REDESIGN.md Sec. 1).  The monthly instrument
% z_t never enters the VAR or the local projections; it is used once, here:
%
%     u_t  = y_t - B' z_{t-1}          VAR(p) innovations at the posterior mean
%     c    = Cov(u_t, z_t)             (K x 1), months with finite z
%     b_z  = c / c(s)                  unit effect on the policy indicator s
%                                      (Stock & Watson 2018)
%     theta(h) = Psi_h b_z             BVAR IRF under proxy identification
%
% The returned struct is a copy of bvar with .b1n := b_z (and .theta,
% .theta_lo/.theta_hi, .b1n_draws recomputed), so estimate_blp_fmar,
% estimate_blp_blockadaptive, estimate_blp_blockpooled and estimate_lp_lagaug
% -- which all read var_est.b1n and form theta_i(h) = b1n' * Gamma_{h,i} --
% are proxy-identified without any change.  Under the maintained VAR(p)
% structure Gamma_h b_z is the LP-IV estimand (Plagborg-Moller & Wolf
% 2021), which the frequentist benchmark estimate_lp_iv estimates directly.
%
% BANDS.  The NIW posterior draws of the BVAR re-identify b_z draw by draw
% (residuals at the drawn coefficients, covariance with the same z), so the
% BVAR bands carry parameter uncertainty in b_z but NOT the sampling error
% of the covariance with z given the coefficients (nor do FMAR's bands for
% the recursive b1n; see docs/APPROXIMATIONS.md).  That sampling error is
% reported separately: a HAC delta-method standard error for each element
% of b_z, the first-stage statistics, and -- fully -- in estimate_lp_iv's
% Anderson-Rubin sets.
%
% INPUTS
% ------
%   bvar : estimate_bvar_niw output (fresh: with .F_draws for bands; a
%          slimmed struct without them gets point estimates and NaN bands).
%   Y    : (T x K) the data the BVAR was estimated on.
%   z    : (T x 1) monthly instrument aligned with the rows of Y (any unit;
%          NaN = not available that month, the row is dropped).
%   cfg  : uses cfg.p, cfg.H, cfg.shock_var (policy indicator, default 1),
%          cfg.ci_level.
%   opts : optional .z_name (string), .hac_lag (default 1), .verbose.
%
% OUTPUT
% ------
%   bz : bvar with .b1n = b_z, .theta, .theta_lo, .theta_hi, .b1n_draws
%        replaced; .b1n_chol = the recursive vector it replaced; .ident:
%          .type 'proxy', .z_name, .s, .N, .cov_uz (K x 1), .b_z, .b_z_se
%          (HAC delta method), .first_stage (b, se_ehw, t_ehw, se_hac,
%          t_hac, F_eff, R2, N, n_nonzero, sd_z, sd_us), .sign_ok,
%          .theta_chol (the recursive IRF, for comparison).
%
% NOTES
% -----
% * Sign: a positive surprise should raise the indicator innovation
%   (c(s) > 0).  If c(s) <= 0 a warning is issued; the normalisation still
%   divides by c(s), so every IRF is then sign-flipped relative to a
%   tightening -- read .sign_ok before interpreting anything.
% * The effective F for a single instrument is the HAC-robust first-stage
%   F (Montiel Olea & Pflueger 2013; critical values 12.05 / 15.06 /
%   23.11 / 37.42 for 30/20/10/5 % Nagar bias).

if nargin < 5, opts = struct(); end
if ~isfield(opts, 'z_name'), opts.z_name = 'z'; end
if ~isfield(opts, 'hac_lag'), opts.hac_lag = 1; end
if ~isfield(opts, 'verbose'), opts.verbose = true; end
[T, K] = size(Y);
p = cfg.p;  H = cfg.H;
s = 1;  if isfield(cfg, 'shock_var') && ~isempty(cfg.shock_var), s = cfg.shock_var; end
alpha_ci = 1 - cfg.ci_level;
assert(numel(z) == T, 'ea_identify_proxy: z must have one entry per row of Y (%d), got %d.', T, numel(z));
z = z(:);

Zall = build_lp_regressors(Y, p);
Zr = Zall(1:end - 1, :);  Y1 = Y(p + 1:T, :);
U  = Y1 - Zr * bvar.B;                         % (N x K) innovations
zz = z(p + 1:T);
ok = isfinite(zz) & all(isfinite(U), 2);
N  = sum(ok);
assert(N > 10 * K, 'ea_identify_proxy: only %d months with a usable instrument.', N);
Uc = U(ok, :) - mean(U(ok, :), 1);
zc = zz(ok) - mean(zz(ok));
c  = (Uc' * zc) / N;                           % Cov(u, z), K x 1
sign_ok = c(s) > 0;
if ~sign_ok
    warning('ea_identify_proxy:sign', ...
            'Cov(u_%d, %s) = %.4g <= 0: a positive surprise LOWERS the indicator innovation; IRFs are sign-flipped.', ...
            s, opts.z_name, c(s));
end
b_z = c / c(s);

% HAC delta-method standard errors of b_z(i) = c_i / c_s
M  = Uc .* zc - ones(N, 1) * c';               % (N x K) moment deviations
S  = M' * M;
L  = opts.hac_lag;
for l = 1:min(L, N - 1)
    G = M(l + 1:end, :)' * M(1:end - l, :);
    S = S + (1 - l / (L + 1)) * (G + G');
end
Vc = S / N^2;                                  % Var(c_hat)
b_z_se = zeros(K, 1);
for i = 1:K
    g = zeros(K, 1);  g(i) = g(i) + 1 / c(s);  g(s) = g(s) - c(i) / c(s)^2;
    b_z_se(i) = sqrt(max(g' * Vc * g, 0));
end
b_z_se(s) = 0;

% first stage: u_s on z (the same slope up to var(z))
fs = slope_stats(zc, Uc(:, s), L);
fs.F_eff = fs.t_hac^2;
fs.N = N;  fs.n_nonzero = sum(zz(ok) ~= 0);
fs.sd_z = std(zz(ok));  fs.sd_us = std(U(ok, s));
fs.corr = corr_(zc, Uc(:, s));

% point IRF and bands
theta = zeros(K, H + 1);
for h = 0:H, theta(:, h + 1) = bvar.Psi(:, :, h + 1) * b_z; end
theta_lo = nan(K, H + 1);  theta_hi = nan(K, H + 1);  b1n_draws = [];
if isfield(bvar, 'F_draws') && ~isempty(bvar.F_draws)
    nd = size(bvar.F_draws, 3);
    J = [eye(K), zeros(K, K * (p - 1))];
    th_d = zeros(K, H + 1, nd);  b1n_draws = zeros(K, nd);
    Zl = Zr(:, 2:end);                          % lag block only: the constant drops out after demeaning
    for d = 1:nd
        Fd = bvar.F_draws(:, :, d);
        Ud = Y1 - Zl * Fd(1:K, :)';             % residuals up to a constant
        Ud = Ud(ok, :);  Ud = Ud - mean(Ud, 1);
        cd = (Ud' * zc) / N;
        bd = cd / cd(s);
        b1n_draws(:, d) = bd;
        Fh = eye(K * p);
        for h = 0:H
            th_d(:, h + 1, d) = (J * Fh * J') * bd;
            Fh = Fh * Fd;
        end
    end
    for i = 1:K
        for h = 1:H + 1
            q = empirical_quantile(squeeze(th_d(i, h, :)), [alpha_ci / 2, 1 - alpha_ci / 2]);
            theta_lo(i, h) = q(1);  theta_hi(i, h) = q(2);
        end
    end
end

bz = bvar;
bz.b1n_chol = bvar.b1n;
bz.theta_chol = bvar.theta;
bz.b1n = b_z;
bz.theta = theta;  bz.theta_lo = theta_lo;  bz.theta_hi = theta_hi;
bz.b1n_draws = b1n_draws;
bz.ident = struct('type', 'proxy', 'z_name', opts.z_name, 's', s, 'N', N, 'cov_uz', c, ...
                  'b_z', b_z, 'b_z_se', b_z_se, 'first_stage', fs, 'sign_ok', sign_ok, ...
                  'hac_lag', L, 'theta_chol', bvar.theta);
if opts.verbose
    fprintf('ea_identify_proxy: %s on %d months; first stage b = %.4f (HAC t = %.2f, F_eff = %.2f, R2 = %.3f); sign %s\n', ...
            opts.z_name, N, fs.b, fs.t_hac, fs.F_eff, fs.R2, tern(sign_ok, 'ok', 'WRONG'));
    fprintf('  b_z (unit effect on var %d):', s);  fprintf(' %.3f', b_z);  fprintf('\n  se(b_z):                  ');
    fprintf(' %.3f', b_z_se);  fprintf('\n');
end
end

% =============================================================================
function st = slope_stats(x, y, L)
X = [ones(numel(x), 1), x(:)];
b = X \ y(:);  u = y(:) - X * b;
XtXi = inv(X' * X);                                                      %#ok<MINV>
Gs = X .* u;  S0 = Gs' * Gs;
V_ehw = XtXi * S0 * XtXi;
S = S0;  n = numel(u);
for l = 1:min(L, n - 1)
    G = Gs(l + 1:end, :)' * Gs(1:end - l, :);
    S = S + (1 - l / (L + 1)) * (G + G');
end
V_hac = XtXi * S * XtXi;
st.b = b(2);  st.se_ehw = sqrt(V_ehw(2, 2));  st.t_ehw = b(2) / st.se_ehw;
st.se_hac = sqrt(V_hac(2, 2));  st.t_hac = b(2) / st.se_hac;
st.R2 = 1 - sum(u.^2) / sum((y - mean(y)).^2);
end

function r = corr_(a, b)
c = corrcoef(a, b);  r = c(1, 2);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
