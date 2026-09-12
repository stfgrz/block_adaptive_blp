function rel = ea_relevance(Y, bvar, cfg)
% EA_RELEVANCE  Instrument-relevance diagnostics for the euro-area system.
%
% PURPOSE
% -------
% The surprise is ordered first and identified recursively, i.e. as an
% INTERNAL INSTRUMENT (Plagborg-Moller & Wolf 2021).  Its relevance for
% the policy indicator is therefore NOT the naive regression of the
% monthly rate change on the surprise (that is contaminated by every
% other piece of news in the month, and by the aggregation of a mid-month
% event into a monthly AVERAGE rate); it is the regression of the rate
% INNOVATION on the surprise INNOVATION, both taken from the VAR that
% conditions on the lags.  That coefficient is exactly the impact
% response the 25 bp normalisation divides by:
%
%     e_i1y(t) = a + b e_mps(t) + v(t),     b = Sigma(2,1) / Sigma(1,1)
%                                             = bvar.b1n(2) = bvar.theta(2,1).
%
% Both are reported, with heteroskedasticity-robust t statistics, so that
% the strength of the identification is a NUMBER in the chapter and not
% an assumption.  With T = 240 monthly observations and a meetings-only
% surprise, the design document budgets an F of about 4 (AGKL, 1M OIS on
% both sides); a smaller value here is expected because the indicator is
% the monthly average of the 12-month Euribor, not an end-of-window OIS
% rate.  A value near zero means the IRFs are weakly identified and the
% 25 bp scale factor is a noisy number; the tau diagnostic, which is a
% statement about the DYNAMIC specification of the whole system, does not
% depend on it.
%
% INPUTS
% ------
%   Y    : (T x K) the dataset, surprise in column 1, rate in column 2.
%   bvar : estimate_bvar_niw(Y, cfg) output (uses .B, .b1n, .theta).
%   cfg  : uses cfg.p.
%
% OUTPUT
% ------
%   rel : struct
%     .impact_b, .impact_t, .impact_F   VAR-innovation regression
%     .sd_innov_mps_bp, .sd_innov_i1y_bp
%     .naive_b, .naive_t, .naive_F      d i1y_t on mps_t (both in pp)
%     .k25                              0.25 / impact response of i1y
%     .p, .T, .N

[T, K] = size(Y);  %#ok<ASGLU>
p = cfg.p;
Zall = build_lp_regressors(Y, p);
Z  = Zall(1:end - 1, :);
Y1 = Y(p + 1:T, :);
E  = Y1 - Z * bvar.B;                     % one-step residuals at the posterior mean
[b, t] = ehw_slope(E(:, 1), E(:, 2));
rel.impact_b = b;  rel.impact_t = t;  rel.impact_F = t^2;
rel.sd_innov_mps_bp = 100 * std(E(:, 1));
rel.sd_innov_i1y_bp = 100 * std(E(:, 2));
% b1n(2) is Sigma(2,1)/Sigma(1,1) from cov(E): the same slope.
assert(abs(rel.impact_b - bvar.b1n(2)) < 1e-8, ...
       'ea_relevance: impact slope %.6f differs from bvar.b1n(2) = %.6f', ...
       rel.impact_b, bvar.b1n(2));

[bn, tn] = ehw_slope(Y(2:end, 1), Y(2:end, 2) - Y(1:end - 1, 2));
rel.naive_b = bn;  rel.naive_t = tn;  rel.naive_F = tn^2;

rel.k25 = 0.25 / bvar.theta(2, 1);
rel.p = p;  rel.T = T;  rel.N = size(E, 1);
end

% -------------------------------------------------------------------------
function [b, t] = ehw_slope(x, y)
X = [ones(numel(x), 1), x(:)];
bb = X \ y(:);
u  = y(:) - X * bb;
XtXi = inv(X' * X);                                        %#ok<MINV>
V = XtXi * (X' * (X .* (u.^2 * ones(1, 2)))) * XtXi;
b = bb(2);  t = bb(2) / sqrt(V(2, 2));
end
