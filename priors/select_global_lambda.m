function [lambda_star, logml] = select_global_lambda(y, Z, mu, d, block_id, cfg)
% PURPOSE
% -------
% Choose the global prior tightness lambda_h for ONE LP equation and
% ONE horizon by maximising the Gaussian marginal likelihood on a grid,
% conditional on (i) sigma^2 fixed at its OLS estimate and (ii) all
% local scales tau_g = 1.  This is a deliberately simple, transparent
% empirical-Bayes device.  It is MODULAR: the eventual thesis should
% swap this file for the global-tightness procedure of the published
% BLP (Ferreira / Miranda-Agrippino / Ricco) without touching the
% estimators.
%
% MODEL / EQUATIONS
% -----------------
% Under the prior beta ~ N(mu, W(lambda)) with
%     W = diag(w),  w_j = lambda^2 d_j  (blocks),  w_1 = d_1 (intercept),
% and the Gaussian likelihood y | beta ~ N(Z beta, sigma2 I), the
% marginal distribution of the data is
%
%     y ~ N( Z mu,  C ),   C = sigma2 I_N + Z W Z',
%
% so with r = y - Z mu the log marginal likelihood is
%
%     log p(y | lambda) = -N/2 log(2 pi) - 1/2 log det C - 1/2 r' C^{-1} r,
%
% evaluated via the Cholesky factor of C (log det C = 2 sum log diag L;
% r' C^{-1} r = || L \ r ||^2 ).
%
% INPUTS
% ------
% y        : (N x 1) dependent variable of this LP regression.
% Z        : (N x m) regressors.
% mu       : (m x 1) VAR-implied prior centre.
% d        : (m x 1) diagonal of D for this equation.
% block_id : (m x 1) block labels (0 = intercept, untouched by lambda).
% cfg      : uses cfg.blp.lambda_grid.
%
% OUTPUTS
% -------
% lambda_star : scalar, grid value with the highest marginal likelihood.
% logml       : (n_grid x 1) log marginal likelihoods (diagnostics).
%
% DIMENSIONS
% ----------
% N x N Cholesky per grid point; N <= 200 in the prototype, cheap.
%
% NOTES
% -----
% Fixing sigma2 at OLS and tau = 1 makes this an approximation to a
% full empirical-Bayes treatment; documented prototype simplification.

[N, m] = size(Z);
grid = cfg.blp.lambda_grid(:);
n_grid = numel(grid);

% OLS residual variance used as the fixed sigma2.
beta_ols = Z \ y;
u = y - Z * beta_ols;
sigma2 = (u' * u) / max(N - m, 1);

r = y - Z * mu;
logml = zeros(n_grid, 1);
for k = 1:n_grid
    lam = grid(k);
    w = lam^2 .* d;
    w(block_id == 0) = d(block_id == 0);        % intercept unscaled
    C = sigma2 * eye(N) + Z * diag(w) * Z';
    L = safe_chol_lower(C);
    v = L \ r;
    logml(k) = -0.5 * N * log(2*pi) - sum(log(diag(L))) - 0.5 * (v' * v);
end

[~, kstar] = max(logml);
lambda_star = grid(kstar);
end
