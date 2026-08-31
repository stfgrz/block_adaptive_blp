function res = var_deterministic_trend(Y, B, p)
% PURPOSE
% -------
% Remove the VAR-implied DETERMINISTIC COMPONENT from the data, as done
% by the FMAR Bayesian Local Projection before running the horizon
% regressions (their child function handleVARtrend).  The h > 1 BLP
% regressions of FMAR are run on the detrended data and the prior on
% the LP constant is then centred at ZERO (with a very loose variance).
%
% MODEL / EQUATIONS
% -----------------
% Given fitted VAR coefficients B = [c'; A_1'; ...; A_p'] (m x K, the
% column-per-equation stacking used everywhere in this project), the
% deterministic component is the path the VAR would follow WITHOUT
% shocks, started from the first p actual observations:
%     trend(t) = Y(t)                                   t = 1..p,
%     trend(t) = c + sum_j A_j trend(t-j)               t = p+1..T.
% As t grows, trend(t) converges to the unconditional mean of the
% fitted VAR; near the start it also absorbs initial-condition
% dynamics.  The detrended series is x(t) = Y(t) - trend(t).
%
% *** Deviation from the FMAR code, documented ***
% FMAR pad the first p detrended observations with ZEROS to preserve
% the sample length (their `res.detrended=[zeros(nL,n);x]`), so a few
% early LP regressor rows contain zeros.  Here the first p rows are
% DROPPED instead (x has T - p rows): cleaner, at the cost of p
% observations.  Since detrended data have mean ~ 0, the two choices
% are numerically very close; neither affects the estimand.
%
% INPUTS
% ------
% Y : (T x K) raw data.
% B : (m x K) VAR coefficients, m = 1 + K*p (intercept in row 1).
% p : VAR lag order.
%
% OUTPUTS
% -------
% res : struct with fields
%   .x     ((T-p) x K) detrended data (rows correspond to t = p+1..T)
%   .trend (T x K)     the deterministic component (first p rows = Y)
%
% DIMENSIONS
% ----------
% T = size(Y,1); K = size(Y,2); m = 1 + K*p.
%
% NOTES
% -----
% The trend uses the POINT ESTIMATE of B (FMAR use the BVAR posterior
% mean); trend uncertainty is ignored, exactly as in FMAR.

[T, K] = size(Y);
m = 1 + K * p;
assert(size(B, 1) == m && size(B, 2) == K, ...
    'var_deterministic_trend: B must be (1+K*p) x K.');
assert(T > p, 'var_deterministic_trend: T must exceed p.');

c = B(1, :)';                          % (K x 1)
A = zeros(K, K, p);
for j = 1:p
    A(:, :, j) = B(1 + (j-1)*K + (1:K), :)';
end

trend = zeros(T, K);
trend(1:p, :) = Y(1:p, :);             % start from actual history
for t = p+1:T
    v = c;
    for j = 1:p
        v = v + A(:, :, j) * trend(t - j, :)';
    end
    trend(t, :) = v';
end

res.trend = trend;
res.x     = Y(p+1:end, :) - trend(p+1:end, :);

assert(all(isfinite(res.x(:))), 'var_deterministic_trend: non-finite output.');
end
