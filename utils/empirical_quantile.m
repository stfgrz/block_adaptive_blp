function q = empirical_quantile(x, probs)
% PURPOSE
% -------
% Empirical quantiles of a sample by linear interpolation of the order
% statistics.  Implemented locally because MATLAB's quantile/prctile
% live in the Statistics Toolbox, which this project avoids.
%
% MODEL / EQUATIONS
% -----------------
% Sort x(1) <= ... <= x(n).  For probability p define the (continuous)
% index  i = 1 + p*(n-1)  and linearly interpolate between the two
% neighbouring order statistics (MATLAB "method 7" / default type).
%
% INPUTS
% ------
% x     : (n x 1) or (1 x n) numeric vector (NaNs are removed).
% probs : vector of probabilities in [0,1].
%
% OUTPUTS
% -------
% q : vector, same shape as probs, of empirical quantiles.
%
% DIMENSIONS
% ----------
% n = numel(x) after NaN removal; requires n >= 1.
%
% NOTES
% -----
% Deterministic; no randomness.

x = x(:);
x = x(~isnan(x));
n = numel(x);
assert(n >= 1, 'empirical_quantile: empty sample.');
assert(all(probs >= 0 & probs <= 1), 'empirical_quantile: probs in [0,1].');

xs = sort(x);
q  = zeros(size(probs));
for k = 1:numel(probs)
    idx = 1 + probs(k) * (n - 1);       % continuous index in [1, n]
    lo  = floor(idx);
    hi  = ceil(idx);
    w   = idx - lo;                     % interpolation weight
    q(k) = (1 - w) * xs(lo) + w * xs(hi);
end
end
