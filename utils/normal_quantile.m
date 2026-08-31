function z = normal_quantile(p)
% PURPOSE
% -------
% Quantile (inverse CDF) of the standard normal distribution, used for
% frequentist confidence intervals.  Implemented with base-MATLAB erfinv
% so the Statistics Toolbox (norminv) is not required.
%
% MODEL / EQUATIONS
% -----------------
% If Phi is the N(0,1) CDF then
%     Phi^{-1}(p) = sqrt(2) * erfinv(2p - 1).
%
% INPUTS
% ------
% p : scalar or array of probabilities in (0,1).
%
% OUTPUTS
% -------
% z : same size as p, standard normal quantiles.
%
% DIMENSIONS
% ----------
% Elementwise.
%
% NOTES
% -----
% Deterministic.

assert(all(p(:) > 0 & p(:) < 1), 'normal_quantile: p must be in (0,1).');
z = sqrt(2) .* erfinv(2 .* p - 1);
end
