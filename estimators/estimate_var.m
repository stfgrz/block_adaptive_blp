function var_est = estimate_var(Y, cfg)
% PURPOSE
% -------
% Estimate a VAR(p) by equation-by-equation OLS (identical to matrix
% OLS), fully in the open: design matrix, coefficient extraction,
% companion matrix, stability diagnostics, impulse responses and simple
% simulated confidence bands.  No toolbox calls.
%
% MODEL / EQUATIONS
% -----------------
% Written as the h = 1 local projection on the shared regressor vector
%     z(t) = [1; y(t); y(t-1); ...; y(t-p+1)],   m = 1 + K*p:
%
%     y(t+1)' = z(t)' * B + e(t+1)',    B is (m x K),
%
% so column i of B holds equation i's coefficients and
%     B = [c, A_1, ..., A_p]'  blockwise:  c' = B(1,:),  A_j = B(2+(j-1)K : 1+jK, :)'.
%
% OLS:      B_hat = (Z'Z)^{-1} Z'Y1  computed as Z \ Y1 (stable solve).
% Sigma:    Sigma_hat = U'U / (N - m), U = OLS residuals (N x K).
% Companion (Kp x Kp):
%     F = [A_1 ... A_p; I_{K(p-1)}  0],   stability: max|eig(F)| < 1.
% Reduced-form IRFs: Psi_h = J F^h J', J = [I_K 0], h = 0..H.
% Structural IRF (recursive identification, simulation device only):
%     B0_hat = chol(Sigma_hat,'lower'),  b1n = B0_hat(:,s)/B0_hat(s,s),
%     theta(:,h+1) = Psi_h * b1n         (unit impact on variable s).
%
% Confidence bands (prototype): draw n_ci_sim coefficient matrices from
% the asymptotic OLS distribution
%     vec(B) ~ N(vec(B_hat), Sigma_hat (x) (Z'Z)^{-1}),
% holding Sigma_hat (hence b1n) at its point estimate, map each draw to
% IRFs, and take empirical quantiles.  This ignores uncertainty in
% Sigma_hat -- documented simplification of the prototype.
%
% INPUTS
% ------
% Y   : (T x K) data.
% cfg : configuration struct (cfg.p, cfg.H, cfg.shock_var, cfg.ci_level,
%       cfg.var.n_ci_sim).
%
% OUTPUTS
% -------
% var_est : struct with fields
%   .B        (m x K)   stacked OLS coefficients [c; A_1'; ...; A_p']-wise
%   .c        (K x 1)   intercept
%   .A        (K x K x p) AR matrices
%   .Sigma    (K x K)   residual covariance (dof-corrected)
%   .F        (Kp x Kp) companion matrix
%   .max_eig  scalar    max |eigenvalue| of F
%   .is_stable logical  max_eig < 1
%   .Psi      (K x K x (H+1)) reduced-form IRFs
%   .b1n      (K x 1)   normalised impact vector of the shock
%   .theta    (K x (H+1)) structural IRF point estimates
%   .theta_lo/.theta_hi (K x (H+1)) simulated CI bounds
%   .resid    (N x K),  .N, .m, .ZtZ
%
% DIMENSIONS
% ----------
% N = T - p observations; m = 1 + K*p regressors.
%
% NOTES
% -----
% Cholesky ordering with the shocked variable first (default s = 1) is
% purely a simulation-device identification choice.

[T, K] = size(Y);
p = cfg.p;  H = cfg.H;  s = cfg.shock_var;

% --- Design and OLS ---------------------------------------------------
Zall = build_lp_regressors(Y, p);       % ((T-p+1) x m), row r <-> t=p+r-1
Z    = Zall(1:end-1, :);                % t = p .. T-1
Y1   = Y(p+1:T, :);                     % y(t+1)
[N, m] = size(Z);
assert(N > m + 10, 'estimate_var: too few observations for OLS.');

B = Z \ Y1;                             % (m x K), QR-based stable solve
U = Y1 - Z * B;                         % (N x K) residuals
Sigma = (U' * U) / (N - m);             % dof-corrected covariance

% --- Unpack coefficients ---------------------------------------------
c = B(1, :)';                           % (K x 1)
A = zeros(K, K, p);
for j = 1:p
    A(:, :, j) = B(1 + (j-1)*K + (1:K), :)';   % (K x K)
end

% --- Companion matrix and stability -----------------------------------
F = zeros(K * p);
for j = 1:p
    F(1:K, (j-1)*K + (1:K)) = A(:, :, j);
end
if p > 1
    F(K+1:end, 1:K*(p-1)) = eye(K * (p - 1));
end
max_eig   = max(abs(eig(F)));
is_stable = max_eig < 1;

% --- Reduced-form and structural IRFs ---------------------------------
J   = [eye(K), zeros(K, K * (p - 1))];
Psi = zeros(K, K, H + 1);
Fh  = eye(K * p);
for h = 0:H
    Psi(:, :, h + 1) = J * Fh * J';
    Fh = Fh * F;
end
B0_hat = safe_chol_lower(Sigma);
b1n    = B0_hat(:, s) / B0_hat(s, s);
theta  = zeros(K, H + 1);
for h = 0:H
    theta(:, h + 1) = Psi(:, :, h + 1) * b1n;
end

% --- Simulated confidence bands ---------------------------------------
[theta_lo, theta_hi] = var_irf_ci(B, Z, Sigma, b1n, K, p, H, cfg);

% --- Pack outputs ------------------------------------------------------
var_est.B         = B;
var_est.c         = c;
var_est.A         = A;
var_est.Sigma     = Sigma;
var_est.F         = F;
var_est.max_eig   = max_eig;
var_est.is_stable = is_stable;
var_est.Psi       = Psi;
var_est.b1n       = b1n;
var_est.theta     = theta;
var_est.theta_lo  = theta_lo;
var_est.theta_hi  = theta_hi;
var_est.resid     = U;
var_est.N         = N;
var_est.m         = m;
var_est.ZtZ       = Z' * Z;

assert(all(isfinite(theta(:))), 'estimate_var: non-finite IRFs.');
end

% =====================================================================
function [lo, hi] = var_irf_ci(B, Z, Sigma, b1n, K, p, H, cfg)
% Simulate from vec(B) ~ N(vec(B_hat), Sigma (x) (Z'Z)^{-1}) using the
% matrix-normal identity: if G is (m x K) iid N(0,1) then
%     B_hat + Rz \ G * Rs   has the required covariance,
% where Z'Z = Rz'Rz (Rz upper triangular from chol) and Sigma = Rs'Rs.
% Each draw is mapped to IRFs; empirical quantiles give the bands.
n_sim = cfg.var.n_ci_sim;
alpha = 1 - cfg.ci_level;
Rz = chol(Z' * Z);                      % upper: Z'Z = Rz'Rz
Rs = chol((Sigma + Sigma') / 2);        % upper: Sigma = Rs'Rs
m  = size(B, 1);
J  = [eye(K), zeros(K, K * (p - 1))];

theta_sims = zeros(K, H + 1, n_sim);
for isim = 1:n_sim
    Bd = B + (Rz \ randn(m, K)) * Rs;   % coefficient draw (m x K)
    Fd = zeros(K * p);
    for j = 1:p
        Fd(1:K, (j-1)*K + (1:K)) = Bd(1 + (j-1)*K + (1:K), :)';
    end
    if p > 1
        Fd(K+1:end, 1:K*(p-1)) = eye(K * (p - 1));
    end
    Fh = eye(K * p);
    for h = 0:H
        theta_sims(:, h + 1, isim) = (J * Fh * J') * b1n;
        Fh = Fh * Fd;
    end
end

lo = zeros(K, H + 1);  hi = zeros(K, H + 1);
for i = 1:K
    for h = 1:H + 1
        q = empirical_quantile(squeeze(theta_sims(i, h, :)), ...
                               [alpha/2, 1 - alpha/2]);
        lo(i, h) = q(1);  hi(i, h) = q(2);
    end
end
end
