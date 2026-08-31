function bvar = estimate_bvar_niw(Y, cfg)
% PURPOSE
% -------
% Bayesian VAR(p) with a Minnesota Normal - Inverse-Wishart prior and
% global tightness lambda chosen by maximising the marginal likelihood
% with a Gamma hyperprior -- the Giannone-Lenza-Primiceri (2015)
% procedure exactly as wired in the FMAR replication code
% (maxMLikelihoodVAR.m + the VAR block of IRFbayesianLocalProj.m).
% In cfg.mode = 'fmar' this object replaces the plain OLS VAR as
%   (i)   the generator of the VAR-implied LP prior CENTRE,
%   (ii)  the identification device (b1n from its residual covariance),
%   (iii) the "VAR" comparator in the Monte Carlo, and
%   (iv)  the source of the deterministic trend removed before the
%         horizon regressions.
%
% MODEL / EQUATIONS
% -----------------
% Regression (shared regressor convention of this project):
%     y(t+1)' = z(t)' B + e(t+1)',   z(t) = [1; y(t); ...; y(t-p+1)].
% Prior (FMAR/GLP, with their d = K + 2 so d - K - 1 = 1):
%     Sigma ~ IW(diag(psi), d),  psi_j = univariate AR(1) residual
%             variance of variable j (computed on the SAME estimation
%             sample, as in FMAR);
%     vec(B) | Sigma ~ N(vec(b), Sigma (x) Omega),
%     b = 0 except first own lag = 1 for variables flagged as random
%         walks (cfg.fmar.isrw; FALSE for the stationary simulations,
%         matching GLP's white-noise centring for stationary data);
%     omega_1 = Vc (loose constant),
%     omega for lag-l of variable v = lambda^2 / (l^2 * psi_v)
%         (Minnesota lag decay with exponent alpha = 2, ACTIVE for the
%         VAR prior -- unlike the BLP horizons, where FMAR switch the
%         decay off).
% lambda maximises  logML(lambda) + log Gamma(lambda; mode .4, sd .2)
% over [1e-4, 5] (GLP hyperprior; golden-section as in
% select_lambda_fmar.m).
%
% Posterior (Kadiyala & Karlsson 1997 conjugate formulas):
%     Omega_end = (Omega^{-1} + Z'Z)^{-1},
%     B_end     = Omega_end (Omega^{-1} b + Z'Y1),
%     Sigma | Y ~ IW(S_end, a_end),
%     S_end = Ehat'Ehat + diag(psi) + (B_end-b)' Omega^{-1} (B_end-b),
%     a_end = d + N,
%     B | Sigma, Y ~ N(B_end, Sigma (x) Omega_end).
%
% IRFs: point estimate from companion powers of B_end with impact
% vector b1n = normalised Cholesky column of cov(Y1 - Z B_end)
% (identical normalisation to estimate_var.m; FMAR's row-normalised
% upper Cholesky is the same object transposed).  Credible bands from
% NIW posterior draws: for each draw, Sigma ~ IW, B | Sigma matrix-
% normal, b1n recomputed from the residuals at the drawn B (as FMAR
% do), IRFs by companion powers, pointwise quantiles.  Draw stability
% is RECORDED but draws are not rejected (FMAR comment the stability
% check out; reproduced faithfully and reported in .diag).
%
% INPUTS
% ------
% Y   : (T x K) data.
% cfg : configuration struct (cfg.p, cfg.H, cfg.shock_var, cfg.ci_level,
%       cfg.fmar.* : Vc, isrw, lambda_min/max, n_niw_draws).
%
% OUTPUTS
% -------
% bvar : struct compatible with estimate_var output where it matters:
%   .B, .c, .A, .Sigma, .F, .max_eig, .is_stable, .Psi, .b1n,
%   .theta, .theta_lo, .theta_hi, .N, .m      (same meanings), plus
%   .lambda      selected tightness
%   .psi         (K x 1) AR(1) prior scales
%   .Omega_diag  (m x 1) prior coefficient variances at lambda
%   .S_end, .a_end, .Omega_end   NIW posterior parameters
%   .diag.share_unstable_draws
%
% DIMENSIONS
% ----------
% N = T - p observations, m = 1 + K*p regressors, K equations.
%
% NOTES
% -----
% .Sigma is set to the posterior MODE sigmahat (used for b1n and for
% downstream code expecting a covariance); .theta uses B_end.

[T, K] = size(Y);
p = cfg.p;  H = cfg.H;  s = cfg.shock_var;
alpha_ci = 1 - cfg.ci_level;

Zall = build_lp_regressors(Y, p);
Z  = Zall(1:end-1, :);
Y1 = Y(p+1:T, :);
[N, m] = size(Z);
d = K + 2;

% --- psi: univariate AR(1) residual variances (FMAR's sigmaj) ---------
psi = zeros(K, 1);
for j = 1:K
    Xj = [ones(N, 1), Z(:, 1 + j)];     % constant + first own lag y_j(t)
    uj = Y1(:, j) - Xj * (Xj \ Y1(:, j));
    psi(j) = std(uj)^2;                 % FMAR use std(u)^2 (biased dof), kept
end

% --- prior mean --------------------------------------------------------
b = zeros(m, K);
if cfg.fmar.isrw
    b(1 + (1:K) + 0, :) = eye(K);       % first own lag = 1 (rows 2..K+1)
end

% --- lambda by logML + Gamma(mode .4, sd .2) hyperprior ----------------
[gk, gtheta] = gamma_coef(0.4, 0.2);
lo = log(cfg.fmar.lambda_min);  hi = log(cfg.fmar.lambda_max);

    function omega = make_omega(lam)
        omega = zeros(m, 1);
        omega(1) = cfg.fmar.Vc;
        for lag = 1:p
            omega(1 + (lag-1)*K + (1:K)) = lam^2 ./ (lag^2 .* psi);
        end
    end
    function Q = objective(loglam)
        lam = exp(loglam);
        lml = niw_logml(Y1, Z, b, make_omega(lam), psi, d);
        Q = lml + (gk - 1) * log(lam) - lam / gtheta ...
                - gk * log(gtheta) - gammaln(gk);
    end

n_grid = 25;  grid = linspace(lo, hi, n_grid);  Qg = zeros(n_grid, 1);
for i = 1:n_grid, Qg(i) = objective(grid(i)); end
[~, i0] = max(Qg);
a = grid(max(i0 - 1, 1));  c = grid(min(i0 + 1, n_grid));
gr = (sqrt(5) - 1) / 2;
x1 = c - gr * (c - a);  x2 = a + gr * (c - a);
f1 = objective(x1);     f2 = objective(x2);
while (c - a) > 1e-6
    if f1 > f2
        c = x2;  x2 = x1;  f2 = f1;  x1 = c - gr * (c - a);  f1 = objective(x1);
    else
        a = x1;  x1 = x2;  f1 = f2;  x2 = a + gr * (c - a);  f2 = objective(x2);
    end
end
lambda = exp((a + c) / 2);
omega  = make_omega(lambda);

% --- conjugate posterior ----------------------------------------------
[~, B_end, sigmahat] = niw_logml(Y1, Z, b, omega, psi, d);
Omega_end = (diag(1 ./ omega) + Z' * Z) \ eye(m);
Ehat  = Y1 - Z * B_end;
S_end = Ehat' * Ehat + diag(psi) + (B_end - b)' * diag(1 ./ omega) * (B_end - b);
a_end = d + N;

% --- unpack, companion, IRFs at the posterior mean ---------------------
c_vec = B_end(1, :)';
A = zeros(K, K, p);
for j = 1:p
    A(:, :, j) = B_end(1 + (j-1)*K + (1:K), :)';
end
F = zeros(K * p);
for j = 1:p
    F(1:K, (j-1)*K + (1:K)) = A(:, :, j);
end
if p > 1, F(K+1:end, 1:K*(p-1)) = eye(K * (p - 1)); end
max_eig = max(abs(eig(F)));

J = [eye(K), zeros(K, K * (p - 1))];
Psi = zeros(K, K, H + 1);
Fh = eye(K * p);
for h = 0:H
    Psi(:, :, h + 1) = J * Fh * J';
    Fh = Fh * F;
end

Sig_resid = cov(Ehat);                  % FMAR use cov(residuals) for Bzero
B0 = safe_chol_lower((Sig_resid + Sig_resid') / 2);
b1n = B0(:, s) / B0(s, s);
theta = zeros(K, H + 1);
for h = 0:H
    theta(:, h + 1) = Psi(:, :, h + 1) * b1n;
end

% --- credible bands from NIW posterior draws ---------------------------
n_draw = cfg.fmar.n_niw_draws;
Lo = safe_chol_lower((Omega_end + Omega_end') / 2);   % row covariance factor
theta_sims = zeros(K, H + 1, n_draw);
n_unstable = 0;
for isim = 1:n_draw
    Sig_d = draw_iw(S_end, a_end);
    Rs = chol((Sig_d + Sig_d') / 2);                  % upper: Sig = Rs'Rs
    B_d = B_end + Lo * randn(m, K) * Rs;              % matrix-normal draw
    Fd = zeros(K * p);
    for j = 1:p
        Fd(1:K, (j-1)*K + (1:K)) = B_d(1 + (j-1)*K + (1:K), :)';
    end
    if p > 1, Fd(K+1:end, 1:K*(p-1)) = eye(K * (p - 1)); end
    if max(abs(eig(Fd))) >= 1, n_unstable = n_unstable + 1; end
    E_d = Y1 - Z * B_d;
    Sig_rd = cov(E_d);
    B0d = safe_chol_lower((Sig_rd + Sig_rd') / 2);
    b1d = B0d(:, s) / B0d(s, s);
    Fh = eye(K * p);
    for h = 0:H
        theta_sims(:, h + 1, isim) = (J * Fh * J') * b1d;
        Fh = Fh * Fd;
    end
end
theta_lo = zeros(K, H + 1);  theta_hi = zeros(K, H + 1);
for i = 1:K
    for h = 1:H + 1
        q = empirical_quantile(squeeze(theta_sims(i, h, :)), ...
                               [alpha_ci/2, 1 - alpha_ci/2]);
        theta_lo(i, h) = q(1);  theta_hi(i, h) = q(2);
    end
end

% --- pack ---------------------------------------------------------------
bvar.B = B_end;      bvar.c = c_vec;    bvar.A = A;
bvar.Sigma = sigmahat;
bvar.F = F;          bvar.max_eig = max_eig;  bvar.is_stable = max_eig < 1;
bvar.Psi = Psi;      bvar.b1n = b1n;
bvar.theta = theta;  bvar.theta_lo = theta_lo;  bvar.theta_hi = theta_hi;
bvar.N = N;          bvar.m = m;
bvar.lambda = lambda;
bvar.psi = psi;
bvar.Omega_diag = omega;
bvar.S_end = S_end;  bvar.a_end = a_end;  bvar.Omega_end = Omega_end;
bvar.diag.share_unstable_draws = n_unstable / n_draw;

assert(all(isfinite(theta(:))), 'estimate_bvar_niw: non-finite IRFs.');
end
