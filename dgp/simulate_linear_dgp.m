function [Y, eta] = simulate_linear_dgp(dgp, cfg)
% PURPOSE
% -------
% Shared simulation engine for the linear DGPs ('var' and 'varma').
% Kept in ONE place so all three DGP files generate data with identical
% mechanics and differ only in their true parameters.
%
% MODEL / EQUATIONS
% -----------------
%   y(t) = c + sum_{j=1}^{p*} A_j y(t-j) + e(t) [+ M e(t-1) if 'varma'],
%   e(t) = B0 eta(t),  eta(t) ~ N(0, I_K).
% The process is started at zeros and cfg.burnin initial periods are
% discarded so the retained sample is (approximately) a draw from the
% stationary distribution.
%
% INPUTS
% ------
% dgp : struct with .type, .A (K x K x p*), .M, .c, .B0.
% cfg : uses cfg.T (kept sample size) and cfg.burnin.
%       RNG state is controlled by the caller.
%
% OUTPUTS
% -------
% Y   : (T x K) simulated data, burn-in removed, rows = time.
% eta : (K x T) structural shocks for the KEPT sample (useful for
%       debugging; estimators never see these).
%
% DIMENSIONS
% ----------
% Internal simulation length Tsim = cfg.burnin + cfg.T.
%
% NOTES
% -----
% Innovations are drawn as B0 * eta so the recursive/Cholesky structure
% used for identification is true by construction in the DGP.

K    = size(dgp.A, 1);
ps   = size(dgp.A, 3);
Tsim = cfg.burnin + cfg.T;

eta_all = randn(K, Tsim);
E       = dgp.B0 * eta_all;             % reduced-form innovations

y = zeros(K, Tsim);
for t = (ps + 1):Tsim
    v = dgp.c + E(:, t);
    for j = 1:ps
        v = v + dgp.A(:, :, j) * y(:, t - j);
    end
    if strcmp(dgp.type, 'varma')
        v = v + dgp.M * E(:, t - 1);
    end
    y(:, t) = v;
end

Y   = y(:, cfg.burnin + 1:end)';        % (T x K)
eta = eta_all(:, cfg.burnin + 1:end);   % (K x T)

assert(size(Y, 1) == cfg.T && size(Y, 2) == K, ...
    'simulate_linear_dgp: unexpected output dimensions.');
assert(all(isfinite(Y(:))), 'simulate_linear_dgp: non-finite data.');
end
