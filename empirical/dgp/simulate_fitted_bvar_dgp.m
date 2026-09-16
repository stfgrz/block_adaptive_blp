function dgp = simulate_fitted_bvar_dgp(bvar, Ydata, cfg, opts)
% SIMULATE_FITTED_BVAR_DGP  Bootstrap the null "the VAR prior is correct".
%
% PURPOSE
% -------
% Empirical analogue of the Monte Carlo's correct-DGP benchmark: generate
% artificial datasets from the VAR(p) whose posterior mean was estimated on
% the REAL data (estimate_bvar_niw).  Under this null the VAR-centred prior
% is correct by construction, so the distribution of the block scales
% tau_{g,h} across replications is exactly the false-positive / null
% distribution against which the real-data tau map must be read
% (see montecarlo/run_null_calibration.m and DESIGN.md Sec. 7).
%
% MODEL
% -----
%   y_t = c + A_1 y_{t-1} + ... + A_p y_{t-p} + e_t,     t = p+1..T,
% with (c, A_j) the BVAR posterior-mean coefficients (bvar.c, bvar.A) and
% e_t drawn by one of four schemes (opts.method):
%   'resample'  iid resampling of the demeaned in-sample one-step residual
%               ROWS (default; preserves fat tails and the cross-equation
%               covariance, but makes the innovations HOMOSKEDASTIC over
%               time: every simulated month draws from the same pool);
%   'wild'      wild bootstrap, e_t = eta_t * Ehat_t with eta_t Rademacher
%               (+1/-1 with probability 1/2) and Ehat_t the residual of the
%               SAME calendar month t.  Keeps the time profile of the
%               volatility (the euro-area short-rate innovation has a
%               four-times smaller standard deviation after 2012 than in
%               2001-2011) and the contemporaneous covariance at each t,
%               so the null then describes a correctly specified VAR with
%               the DATA's heteroskedasticity.  Needs opts.T = T (the
%               default);
%   'block'     moving-block bootstrap of residual rows, block length
%               opts.block_len (default 12): preserves volatility
%               clustering and within-block serial dependence of the
%               residuals without tying the simulation to the calendar;
%   'gaussian'  N(0, cov(Ehat)).
% The simulation is initialised at the ACTUAL first p observations of the
% real data (standard conditional parametric bootstrap for persistent
% systems; the fitted VAR may have roots near unity under the random-walk
% prior mean, which is fine over a finite T but rules out simulating from
% an unconditional distribution).
%
% INPUTS
% ------
%   bvar  : output of estimate_bvar_niw(Ydata, cfg).  Fields used: c
%           (K x 1 constant), A (K x K x p), B (m x K; constant row first,
%           then p lag blocks -- used only to recompute residuals), theta
%           (true-under-null IRF, passed through for coverage checks).
%   Ydata : (T x K) the real data the BVAR was estimated on.
%   cfg   : uses cfg.p only.
%   opts  : optional struct:
%           .method     'resample' (default) | 'wild' | 'block' | 'gaussian'
%           .T          length of the simulated sample (default size(Ydata,1))
%           .block_len  block length for 'block' (default 12)
%
% OUTPUT
% ------
%   dgp.Y      (T x K) simulated data (first p rows = real data).
%   dgp.theta  (K x (H+1)) the IRF that is true under this null
%              (= bvar.theta), for optional coverage-under-null checks.
%   dgp.name   'fitted_bvar'
%   dgp.p_true cfg.p
%
% NOTES
% -----
% * RNG is controlled by the caller (rng(seed) before each replication).
% * Residuals are recomputed here from (Ydata, bvar.B) rather than stored
%   in the bvar struct; the regressor layout matches estimate_bvar_niw
%   (constant first, lag blocks in order), which is asserted numerically:
%   the two residual constructions (via B and via c/A) must agree.

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'method'), opts.method = 'resample'; end
if ~isfield(opts, 'block_len'), opts.block_len = 12; end
[Td, K] = size(Ydata);
if ~isfield(opts, 'T'), opts.T = Td; end
p = cfg.p;
T = opts.T;
assert(T > p + 10, 'simulate_fitted_bvar_dgp: T too small.');
assert(any(strcmp(opts.method, {'resample', 'wild', 'block', 'gaussian'})), ...
       'simulate_fitted_bvar_dgp: unknown method ''%s''.', opts.method);

% --- recompute in-sample residuals ----------------------------------------
N = Td - p;
Z = ones(N, 1 + K * p);
for j = 1:p
    Z(:, 1 + (j - 1) * K + (1:K)) = Ydata(p + 1 - j:Td - j, :);
end
Y1   = Ydata(p + 1:Td, :);
Ehat = Y1 - Z * bvar.B;

% cross-check the (c, A) representation against B (layout assertion)
E2 = Y1 - ones(N, 1) * bvar.c(:)';
for j = 1:p
    E2 = E2 - Ydata(p + 1 - j:Td - j, :) * bvar.A(:, :, j)';
end
assert(max(abs(Ehat(:) - E2(:))) < 1e-8, ...
       'simulate_fitted_bvar_dgp: bvar.B vs bvar.c/A layout mismatch.');

Ec = Ehat - ones(N, 1) * mean(Ehat, 1);          % centred residuals

switch opts.method
    case 'gaussian'
        Sig = (Ec' * Ec) / (N - 1);
        Lch = safe_chol_lower((Sig + Sig') / 2);
        E = (Lch * randn(K, T - p))';
    case 'resample'
        E = Ec(ceil(rand(T - p, 1) * N), :);
    case 'wild'
        assert(T - p == N, ...
               'simulate_fitted_bvar_dgp: ''wild'' needs opts.T equal to the data length (%d).', Td);
        eta = 2 * (rand(N, 1) > 0.5) - 1;                  % Rademacher
        E = Ec .* eta;                                     % same calendar month, random sign
    case 'block'
        L = max(1, min(round(opts.block_len), N));
        nb = ceil((T - p) / L);
        starts = ceil(rand(nb, 1) * (N - L + 1));
        idx = reshape(bsxfun(@plus, starts, 0:L - 1)', [], 1);
        E = Ec(idx(1:T - p), :);
end

% --- simulate ---------------------------------------------------------------
Y = zeros(T, K);
Y(1:p, :) = Ydata(1:p, :);
for t = p + 1:T
    m = bvar.c(:)';
    for j = 1:p
        m = m + Y(t - j, :) * bvar.A(:, :, j)';
    end
    Y(t, :) = m + E(t - p, :);
end
assert(all(isfinite(Y(:))) && max(abs(Y(:))) < 1e8, ...
       'simulate_fitted_bvar_dgp: simulation exploded (unstable fitted VAR?).');

dgp.Y      = Y;
dgp.theta  = bvar.theta;
dgp.name   = 'fitted_bvar';
dgp.method = opts.method;
dgp.p_true = p;
end
