function blp = estimate_blp_fmar(Y, cfg, bvar)
% PURPOSE
% -------
% GLOBAL VAR-centred Bayesian Local Projection following the published
% methodology of Ferreira, Miranda-Agrippino & Ricco (REStat 2025), as
% implemented in their replication code (IRFbayesianLocalProj.m with
% priorType = 'VAR', identification = 'CHOL', sur = noc = 0,
% MNpsi = 0).  This replaces the prototype estimate_blp_global.m as the
% baseline when cfg.mode = 'fmar'.  Every FMAR ingredient is used:
%
%   1. DETRENDING: the horizon regressions are run on data with the
%      BVAR-implied deterministic component removed
%      (utils/var_deterministic_trend.m), and the prior on the LP
%      constant is centred at 0 with loose variance Vc.
%   2. PRIOR CENTRE: at horizon h, the coefficient block on the lags is
%      the companion power A^h of the BVAR posterior-mean coefficients
%      (FMAR's setPriorMean_VAR; identical to the (J F^h)' mapping of
%      priors/var_implied_lp_prior.m, which tests verify).
%   3. PRIOR SCALE: NIW with IW scale psi(h) = Newey-West LONG-RUN
%      variances of univariate own-lag LP residuals
%      (priors/fmar_prior_scale.m); coefficient prior variance
%      Sigma_ii * lambda_h^2 / psi_v(h), NO Minnesota lag decay.
%   4. GLOBAL TIGHTNESS lambda_h: maximises the NIW marginal likelihood
%      plus FMAR's horizon-dependent Gamma hyperprior
%      (priors/select_lambda_fmar.m).
%   5. POINT ESTIMATE: conjugate posterior mean
%      Bhat_h = (X'X + Omega^{-1})^{-1} (X'y + Omega^{-1} b).
%   6. UNCERTAINTY, h >= 2 ("quasi-Bayesian"): frequentist Newey-West
%      SANDWICH bands around the posterior-mean IRF,
%          V_i = (X'X)^{-1} S_i (X'X)^{-1},
%      with scores z(t) * uhat_i(t) evaluated at the POSTERIOR MEAN
%      residuals (demeaned), Bartlett truncation L = h + 1, and the
%      impact vector b1n held FIXED at the BVAR point estimate.  This
%      is FMAR's ex-post correction for the MA(h-1) serial correlation
%      that the per-horizon Gaussian likelihood ignores.
%   7. h = 0 is always the shared identification step (theta = b1n).
%      What happens at h = 1 is controlled by cfg.fmar.h1_mode:
%        'bvar' (DEFAULT, the FMAR convention): h = 1 is taken from the
%           Bayesian VAR itself (point = posterior mean IRF, bands = NIW
%           posterior-draw quantiles), and lambda(1) is the BVAR
%           tightness -- exactly as in FMAR, where the h <= 1 responses
%           come from the VAR block.
%        'lp' (HARMONISED comparison): h = 1 is estimated as a local
%           projection with the same FMAR machinery as h >= 2 (prior
%           centre A^1, psi(1), lambda_1 from select_lambda_fmar, NW
%           sandwich band with L = 2).
%      WHY THE OPTION EXISTS.  The block-adaptive estimator has always
%      run an LP at h = 1, so under 'bvar' the two estimators are not
%      even the same class of object at h = 1 and any h = 1 RMSE gap
%      confounds "LP vs VAR" with "adaptive vs global".  Under 'lp' the
%      two agree exactly at tau = 1 for EVERY h >= 1, so h = 1 enters
%      the comparison on the same footing as the other horizons.  The
%      Monte Carlo reports an integrated RMSE over h = 2..H (which is
%      unaffected by this switch) as the PREFERRED metric, and the
%      legacy h = 1..H version alongside it.
%
% MODEL / EQUATIONS
% -----------------
% See items 1-7; full formulas live in the called functions.  The
% structural response of variable i at horizon h is
%     theta_i(h) = b1n' * Bhat_h(2:1+K, i),
% the same convention as everywhere in this project.
%
% INPUTS
% ------
% Y    : (T x K) raw data.
% cfg  : configuration struct (cfg.p, cfg.H, cfg.ci_level, cfg.fmar.*,
%        incl. the optional cfg.fmar.h1_mode described above).
% bvar : output of estimate_bvar_niw (prior centre, b1n, detrending,
%        h <= 1 responses and bands, lambda at h = 1).
%
% OUTPUTS
% -------
% blp : struct with the same interface as estimate_blp_global:
%   .theta_mean/.theta_med (K x (H+1)) point IRFs (identical here: the
%                          conjugate posterior is symmetric in beta)
%   .lo/.hi                (K x (H+1)) bands (NIW quantiles at h <= 1,
%                          HAC sandwich at h >= 2)
%   .lambda                (K x H) tightness (constant across equations
%                          at each horizon; replicated for interface
%                          compatibility)
%   .prior_theta           (K x (H+1)) IRF implied by the prior centre
%   .diag                  struct (.lambda_at_bound count, .h1_mode,
%                          .h_start = first horizon estimated as an LP)
%
% DIMENSIONS
% ----------
% Detrended sample Tx = T - p; horizon-h regression has
% N_h = Tx - p + 1 - h observations and m = 1 + K*p regressors.
%
% NOTES
% -----
% No coefficient sampling happens at h >= 2 (FMAR summarise the
% posterior by its mean and attach sandwich bands); the block-adaptive
% extension adds the sampling layer in estimate_blp_blockadaptive.m.

[T, K] = size(Y);
p = cfg.p;  H = cfg.H;
zcrit = normal_quantile(1 - (1 - cfg.ci_level) / 2);
b1n = bvar.b1n;
m = 1 + K * p;
h1_mode = fmar_h1_mode(cfg);

% --- detrend and build shared regressors -------------------------------
dt = var_deterministic_trend(Y, bvar.B, p);
x  = dt.x;                              % ((T-p) x K)
Tx = size(x, 1);
Zall = build_lp_regressors(x, p);       % ((Tx-p+1) x m)

% --- containers ---------------------------------------------------------
theta = zeros(K, H + 1);
lo = nan(K, H + 1);  hi = nan(K, H + 1);
lambda_used = zeros(K, H);
prior_theta = zeros(K, H + 1);
n_at_bound = 0;

% --- h = 0: shared identification step ---------------------------------
theta(:, 1) = bvar.theta(:, 1);   lo(:, 1) = bvar.theta_lo(:, 1);   hi(:, 1) = bvar.theta_hi(:, 1);
prior_theta(:, 1) = b1n;

% --- h = 1: BVAR (FMAR convention) or LP (harmonised comparison) -------
switch h1_mode
    case 'bvar'
        theta(:, 2) = bvar.theta(:, 2);
        lo(:, 2) = bvar.theta_lo(:, 2);   hi(:, 2) = bvar.theta_hi(:, 2);
        lambda_used(:, 1) = bvar.lambda;
        prior_theta(:, 2) = bvar.Psi(:, :, 2) * b1n;
        h_start = 2;
    case 'lp'
        h_start = 1;                  % h = 1 goes through the LP loop
    otherwise
        error('estimate_blp_fmar: unknown cfg.fmar.h1_mode ''%s''.', h1_mode);
end

% --- horizons h >= h_start ----------------------------------------------
J = [eye(K), zeros(K, K * (p - 1))];
for h = h_start:H
    Zh = Zall(1:end - h, :);            % (N_h x m)
    Yh = x(p + h:Tx, :);                % (N_h x K)
    Nh = size(Zh, 1);
    assert(Nh > m + 5, 'estimate_blp_fmar: too few observations at h = %d.', h);

    % prior centre: zero constant + companion power of the BVAR
    Fh = bvar.F^h;
    b_center = [zeros(1, K); (J * Fh)'];        % (m x K)
    prior_theta(:, h + 1) = b_center(2:1 + K, :)' * b1n;

    % prior scales and tightness
    psi_h = fmar_prior_scale(x, p, h);
    [lam, Bhat, ~, info] = select_lambda_fmar(Yh, Zh, b_center, psi_h, h, cfg);
    lambda_used(:, h) = lam;
    n_at_bound = n_at_bound + info.at_bound;

    theta(:, h + 1) = Bhat(2:1 + K, :)' * b1n;

    % FMAR quasi-Bayesian bands: NW sandwich at the posterior mean
    Uh = Yh - Zh * Bhat;
    Uh = Uh - mean(Uh, 1);              % FMAR demean the residuals
    ZtZinv = (Zh' * Zh) \ eye(m);
    L = h + 1;
    w = (L + 1 - (1:L)) / (L + 1);
    for i = 1:K
        Gsc = Zh .* Uh(:, i);           % (N_h x m) score rows
        S = Gsc' * Gsc;
        for l = 1:min(L, Nh - 1)
            Gl = Gsc(l+1:end, :)' * Gsc(1:end-l, :);
            S  = S + w(l) * (Gl + Gl');
        end
        V = ZtZinv * S * ZtZinv;
        sd_th = sqrt(max(b1n' * V(2:1 + K, 2:1 + K) * b1n, 0));
        lo(i, h + 1) = theta(i, h + 1) - zcrit * sd_th;
        hi(i, h + 1) = theta(i, h + 1) + zcrit * sd_th;
    end
end

blp.theta_mean  = theta;
blp.theta_med   = theta;                % conjugate posterior: mean = median
blp.lo          = lo;
blp.hi          = hi;
blp.lambda      = lambda_used;
blp.prior_theta = prior_theta;
blp.diag.lambda_at_bound = n_at_bound;
blp.diag.h1_mode = h1_mode;
blp.diag.h_start = h_start;

assert(all(isfinite(theta(:))), 'estimate_blp_fmar: non-finite estimates.');
end

% =====================================================================
function mode_str = fmar_h1_mode(cfg)
% Resolve cfg.fmar.h1_mode, defaulting to the FMAR convention 'bvar'
% so every earlier configuration struct keeps its exact behaviour.
mode_str = 'bvar';
if isfield(cfg, 'fmar') && isstruct(cfg.fmar) && ...
        isfield(cfg.fmar, 'h1_mode') && ~isempty(cfg.fmar.h1_mode)
    mode_str = cfg.fmar.h1_mode;
end
end
