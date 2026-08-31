function [Zall, t_idx] = build_lp_regressors(Y, p)
% PURPOSE
% -------
% Build the common regressor matrix used by EVERY estimator in this
% project (VAR, LP, global BLP, block-adaptive BLP), so that all of them
% condition on exactly the same information set
%
%     z(t) = [1; y(t); y(t-1); ...; y(t-p+1)]   (m x 1, m = 1 + K*p).
%
% MODEL / EQUATIONS
% -----------------
% For horizon h >= 1 the local projection is
%     y(t+h) = Beta_h' z(t) + u(t+h),
% and the VAR(p) is exactly the h = 1 case:
%     y(t+1) = [c, A_1, ..., A_p] z(t) + e(t+1).
% Using one shared design matrix guarantees that the VAR-implied prior
% centre (priors/var_implied_lp_prior.m) refers to the same regressor
% ordering as the LP design matrix.
%
% INPUTS
% ------
% Y : (T x K) data matrix, rows = time, columns = variables.
% p : scalar, number of lags entering z(t) (y(t) counts as "lag 0").
%
% OUTPUTS
% -------
% Zall  : ((T-p+1) x m) matrix.  Row r corresponds to time t = p + r - 1
%         and contains [1, y(t)', y(t-1)', ..., y(t-p+1)'].
% t_idx : ((T-p+1) x 1) vector of the time indices t for each row.
%
% DIMENSIONS
% ----------
% T = size(Y,1), K = size(Y,2), m = 1 + K*p.
% For a horizon-h regression use rows 1:(T-p+1-h) of Zall paired with
% dependent observations Y(p+h : T, :); there are N_h = T - p + 1 - h
% usable observations.
%
% NOTES
% -----
% Column ordering (intercept first, then y(t), then deeper lags) matches
% the companion-form state ordering used in var_implied_lp_prior.m.
% Do not change one without the other.

[T, K] = size(Y);
assert(p >= 1 && p < T, 'build_lp_regressors: need 1 <= p < T.');

m     = 1 + K * p;
nrow  = T - p + 1;
Zall  = zeros(nrow, m);
t_idx = (p:T)';

Zall(:, 1) = 1;                                  % intercept
for lag = 0:(p-1)
    cols = 1 + lag*K + (1:K);                    % columns for y(t-lag)
    Zall(:, cols) = Y(t_idx - lag, :);
end

assert(all(isfinite(Zall(:))), 'build_lp_regressors: non-finite entries.');
end
