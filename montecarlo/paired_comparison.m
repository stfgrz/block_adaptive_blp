function P = paired_comparison(mc, e_test, e_base, h_range)
% PURPOSE
% -------
% Compare two estimators on the SAME simulated samples and attach a
% standard error to the difference, so that a 1-2% gap in integrated
% RMSE can be called real or not called anything at all.
%
% WHY A PAIRED TEST IS THE RIGHT ONE HERE
% ---------------------------------------
% The Monte Carlo standard error of a SINGLE estimated RMSE is about
% rmse/sqrt(2R) -- 3.2% at R = 500.  Read naively, that would make every
% difference in this study statistically invisible.  But the estimators
% are not run on independent data: every replication feeds the SAME
% simulated sample to all of them (common random numbers), so their
% errors are strongly positively correlated and the DIFFERENCE is
% estimated far more precisely than either level.  The paired standard
% error below is typically an order of magnitude smaller than the
% naive one, and it is what any claim about which estimator wins must
% be judged against.
%
% MODEL / EQUATIONS
% -----------------
% For replication r, with theta0 the truth:
%
%   d_r = mean over (i, h in h_range) of
%           [ (theta_test - theta0)^2 - (theta_base - theta0)^2 ]
%
% i.e. the per-replication difference in mean squared error.  Then
%   dMSE  = mean_r d_r,
%   se    = std(d_r) / sqrt(R),
%   t     = dMSE / se.
% Negative dMSE means the test estimator has the LOWER MSE.  The ratio
% of integrated MSEs and the implied RMSE ratio are reported alongside,
% because those are the numbers the tables quote.
%
% INPUTS
% ------
% mc      : run_montecarlo output.
% e_test  : name or index of the estimator being judged.
% e_base  : name or index of the baseline (default 'BLP-FMAR').
% h_range : horizons to integrate over (default 2:H, the preferred
%           metric).
%
% OUTPUTS
% -------
% P : struct with .dMSE, .se, .t, .mse_test, .mse_base, .mse_ratio,
%     .rmse_ratio, .R, .h_range, .names.
%
% NOTES
% -----
% Deterministic given mc.  The t-statistic is a Monte Carlo precision
% statement about THIS design, not a test about the population of DGPs.

if nargin < 3 || isempty(e_base), e_base = 'BLP-FMAR'; end
[nE, K, Hp1, R] = size(mc.theta);
H = Hp1 - 1;
if nargin < 4 || isempty(h_range), h_range = 2:H; end

it = resolve(mc, e_test);
ib = resolve(mc, e_base);

theta0 = mc.theta_true;
mse_t = zeros(R, 1);  mse_b = zeros(R, 1);
for r = 1:R
    at = 0;  ab = 0;  cnt = 0;
    for i = 1:K
        for h = h_range
            et = mc.theta(it, i, h + 1, r) - theta0(i, h + 1);
            eb = mc.theta(ib, i, h + 1, r) - theta0(i, h + 1);
            at = at + et^2;  ab = ab + eb^2;  cnt = cnt + 1;
        end
    end
    mse_t(r) = at / cnt;  mse_b(r) = ab / cnt;
end
d = mse_t - mse_b;

P.names      = {mc.est_names{it}, mc.est_names{ib}};
P.R          = R;
P.h_range    = h_range;
P.dMSE       = mean(d);
P.se         = std(d) / sqrt(R);
P.t          = P.dMSE / max(P.se, realmin);
P.mse_test   = mean(mse_t);
P.mse_base   = mean(mse_b);
P.mse_ratio  = P.mse_test / P.mse_base;
P.rmse_ratio = sqrt(P.mse_ratio);
% What the se of the same difference would have been if the two
% estimators had been run on INDEPENDENT samples.  The ratio of the two
% is how much the common-random-numbers design buys.
P.se_unpaired = sqrt(var(mse_t) + var(mse_b)) / sqrt(R);
P.pairing_gain = P.se_unpaired / max(P.se, realmin);
% The same difference expressed on the RMSE-ratio scale, with a
% one-standard-error band, since that is the scale the tables use.
lo_mse = (P.mse_base + P.dMSE - P.se) / P.mse_base;
hi_mse = (P.mse_base + P.dMSE + P.se) / P.mse_base;
P.rmse_ratio_lo = sqrt(max(lo_mse, 0));
P.rmse_ratio_hi = sqrt(max(hi_mse, 0));
P.correlation = corr_(mse_t, mse_b);
end

% =====================================================================
function r = corr_(x, y)
x = x - mean(x);  y = y - mean(y);
r = (x' * y) / max(sqrt((x' * x) * (y' * y)), realmin);
end

% =====================================================================
function idx = resolve(mc, e)
if ischar(e)
    idx = find(strcmp(mc.est_names, e), 1);
    assert(~isempty(idx), 'paired_comparison: no estimator named "%s".', e);
else
    idx = e;
end
end
