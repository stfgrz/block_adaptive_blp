function blp = estimate_blp_blockadaptive(Y, cfg, var_est, lambda_mat)
% PURPOSE
% -------
% BLOCK-ADAPTIVE VAR-centred Bayesian Local Projection (the
% methodological object of the thesis).  Relative to a global BLP, each
% coefficient block g gets its own local "escape" scale tau_{g,h}:
%
%     beta_{g,h} - mu_{g,h}^VAR ~ N(0, [sigma2] lambda_h^2 tau_{g,h}^2 D_g),
%     tau_{g,h} ~ half-Cauchy(0, 1),
%
% where the bracketed sigma2 appears only in FMAR mode (see below).
% Interpretation: tau small  => block g stays close to the VAR prior;
%                 tau large  => block g is allowed to escape from it.
% The shrinkage is on DEVIATIONS FROM THE VAR CENTRE, not towards zero.
%
% TWO OPERATING MODES (cfg.mode)
% ------------------------------
% 'prototype' (original behaviour, unchanged):
%     raw data; prior centre from priors/var_implied_lp_prior.m (OLS
%     VAR); D from priors/build_block_prior.m (std-ratio scaling);
%     lambda from the grid marginal likelihood or fixed; sampler with
%     prior NOT scaling with sigma2; bands = posterior quantiles.
%
% 'fmar' (nested in the published FMAR baseline):
%     var_est must be the output of estimate_bvar_niw.  Everything is
%     inherited from the FMAR machinery so that fixing tau_g = 1
%     reproduces estimate_blp_fmar's posterior mean EXACTLY:
%       * data DETRENDED with the BVAR deterministic component;
%       * prior centre = companion power of the BVAR (constant at 0);
%       * relative scales d_j = 1/psi_v(h) with psi from
%         priors/fmar_prior_scale.m (NW long-run LP variances), and
%         d_intercept = Vc; the prior variance scales with sigma2
%         (opts.prior_scales_with_sigma2 = true), matching the NIW
%         prior variance Sigma_ii * lambda^2 / psi_v;
%       * per-equation sigma2 prior IG(a0, b0) with a0 = (d-K+1)/2 = 3/2
%         and b0 = psi_i(h)/2: this is the exact MARGINAL of the ith
%         diagonal of Sigma ~ IW(diag(psi), d = K+2), so the equation-
%         by-equation sampler mimics the FMAR system prior as closely
%         as a per-equation model can (the one approximation: the
%         sampler ignores cross-equation Sigma correlation, which does
%         not affect the tau = 1 posterior MEAN of beta -- that is
%         sigma2-free -- but does affect posterior spread);
%       * lambda_h from priors/select_lambda_fmar.m (or passed in via
%         lambda_mat, e.g. estimate_blp_fmar's lambda so both use the
%         SAME tightness -- recommended in the Monte Carlo);
%       * PRIMARY bands = the same FMAR quasi-Bayesian NW sandwich as
%         estimate_blp_fmar, centred at the block-adaptive posterior
%         mean, so coverage differences across the two estimators
%         reflect the point estimator only; posterior-quantile bands
%         are returned as SECONDARY output (.lo_post/.hi_post).
%       * h = 0 is the shared identification step (theta = b1n); h >= 1
%         are estimated as LPs.  NOTE: at h = 1 the FMAR baseline
%         reports the BVAR itself, while this estimator runs the h = 1
%         LP with the BVAR-centred prior; the two are close but not
%         identical -- documented design choice, tested at h >= 2.
%
% MODEL / EQUATIONS
% -----------------
% Estimated equation by equation and horizon by horizon with the Gibbs
% sampler in samplers/gibbs_block_horseshoe.m; ALL conditionals are
% documented there (primed set for FMAR mode).  Blocks: scheme
% 'per_variable' (block g = all lag coefficients of variable g).
%
% INPUTS
% ------
% Y          : (T x K) raw data.
% cfg        : configuration struct (cfg.mode selects the behaviour).
% var_est    : estimate_var output ('prototype') or estimate_bvar_niw
%              output ('fmar').
% lambda_mat : OPTIONAL (K x H) tightness per equation/horizon (e.g.
%              blp_global.lambda or blp_fmar.lambda).  If omitted or
%              empty, lambda is selected inside, mode-appropriately.
%
% OUTPUTS
% -------
% blp : struct with fields
%   .theta_mean/.theta_med (K x (H+1)) posterior mean/median IRFs
%   .lo/.hi                (K x (H+1)) PRIMARY bands (posterior
%                          quantiles in 'prototype'; FMAR sandwich in
%                          'fmar')
%   .lo_post/.hi_post      (K x (H+1)) posterior-quantile bands (only
%                          in 'fmar'; equals .lo/.hi in 'prototype')
%   .lambda                (K x H) tightness actually used
%   .prior_theta           (K x (H+1)) IRF implied by the prior centre
%   .tau_mean/.tau_med     (K x G x H) posterior block-scale summaries
%   .diag                  sampler diagnostics
%
% DIMENSIONS
% ----------
% G = K blocks of p coefficients each; m = 1 + K*p regressors.
%
% NOTES
% -----
% Posterior tau summaries are on the tau (not tau^2) scale.  RNG is
% controlled by the caller.

if nargin < 4
    lambda_mat = [];
end
if ~isfield(cfg, 'mode'), cfg.mode = 'prototype'; end
fmar = strcmp(cfg.mode, 'fmar');

[T, K] = size(Y);
p = cfg.p;  H = cfg.H;
alpha = 1 - cfg.ci_level;
zcrit = normal_quantile(1 - alpha / 2);
b1n = var_est.b1n;
m = 1 + K * p;

% --- mode-specific data, blocks, centres --------------------------------
if fmar
    dt = var_deterministic_trend(Y, var_est.B, p);
    x  = dt.x;                          % detrended, ((T-p) x K)
else
    x  = Y;
end
Tx = size(x, 1);
Zall = build_lp_regressors(x, p);

% block membership is the same in both modes; D differs.
bp = build_block_prior(Y, cfg);         % uses raw-data std ratios ('prototype')
G  = bp.G;

gopts.n_burn     = cfg.gibbs.n_burn;
gopts.n_keep     = cfg.gibbs.n_keep;
gopts.sample_tau = true;
gopts.tau2_min   = cfg.gibbs.tau2_min;
gopts.tau2_max   = cfg.gibbs.tau2_max;
gopts.prior_scales_with_sigma2 = fmar;
% Nesting switch used by tests: force tau to a fixed value.
if isfield(cfg.blp, 'fix_tau') && ~isempty(cfg.blp.fix_tau)
    gopts.sample_tau = false;
    gopts.tau_fixed  = cfg.blp.fix_tau;
end

theta_mean = zeros(K, H + 1);  theta_med = zeros(K, H + 1);
lo = nan(K, H + 1);  hi = nan(K, H + 1);
lo_post = nan(K, H + 1);  hi_post = nan(K, H + 1);
theta_mean(:, 1) = b1n;  theta_med(:, 1) = b1n;
prior_theta = zeros(K, H + 1);  prior_theta(:, 1) = b1n;
lambda_used = zeros(K, H);
tau_mean = zeros(K, G, H);  tau_med = zeros(K, G, H);
lag1 = zeros(K, H);  n_clip = 0;
draws = cell(K, H);

J = [eye(K), zeros(K, K * (p - 1))];

for h = 1:H
    Zh = Zall(1:end - h, :);
    Yh = x(p + h:Tx, :);
    Nh = size(Zh, 1);

    % ---- prior centre --------------------------------------------------
    if fmar
        Fh = var_est.F^h;
        Mu = [zeros(1, K); (J * Fh)'];  % detrended data: constant centred at 0
    else
        Mu = var_implied_lp_prior(var_est, h);
    end
    prior_theta(:, h + 1) = Mu(2:1 + K, :)' * b1n;

    % ---- horizon-level scale objects (FMAR) -----------------------------
    if fmar
        psi_h = fmar_prior_scale(x, p, h);
        % lambda common to all equations at this horizon:
        if isempty(lambda_mat)
            lam_h = select_lambda_fmar(Yh, Zh, Mu, psi_h, h, cfg);
        end
        % objects reused by the sandwich bands:
        ZtZinv = (Zh' * Zh) \ eye(m);
        Ltr = h + 1;                     % FMAR truncation
        wNW = (Ltr + 1 - (1:Ltr)) / (Ltr + 1);
    end

    for i = 1:K
        yh = Yh(:, i);
        prior.mu       = Mu(:, i);
        prior.block_id = bp.block_id;

        if fmar
            % d_j = 1/psi_v for a coefficient on any lag of variable v;
            % intercept d_1 = Vc.  Combined with the sigma2 scaling this
            % reproduces the NIW prior variance Sigma_ii * lambda^2/psi_v.
            d = zeros(m, 1);
            d(1) = cfg.fmar.Vc;
            for lag = 1:p
                d(1 + (lag-1)*K + (1:K)) = 1 ./ psi_h;
            end
            prior.d  = d;
            % sigma2_i prior: exact IW-diagonal marginal, IG(3/2, psi_i/2)
            prior.a0 = (K + 2 - K + 1) / 2;          % = 3/2 for d = K+2
            prior.b0 = psi_h(i) / 2;
            if ~isempty(lambda_mat)
                prior.lambda = lambda_mat(i, h);
            else
                prior.lambda = lam_h;
            end
        else
            prior.d  = bp.d_all(:, i);
            prior.a0 = cfg.gibbs.a0;
            prior.b0 = cfg.gibbs.b0;
            if ~isempty(lambda_mat)
                prior.lambda = lambda_mat(i, h);
            else
                switch cfg.blp.lambda_mode
                    case 'fixed'
                        prior.lambda = cfg.blp.lambda_fixed;
                    case 'grid'
                        prior.lambda = select_global_lambda(yh, Zh, ...
                            prior.mu, prior.d, prior.block_id, cfg);
                    otherwise
                        error('estimate_blp_blockadaptive: unknown lambda_mode.');
                end
            end
        end
        lambda_used(i, h) = prior.lambda;

        out = gibbs_block_horseshoe(yh, Zh, prior, gopts);

        th_draws = out.beta_draws(:, 2:1 + K) * b1n;
        theta_mean(i, h + 1) = mean(th_draws);
        theta_med(i, h + 1)  = empirical_quantile(th_draws, 0.5);
        q = empirical_quantile(th_draws, [alpha/2, 1 - alpha/2]);
        lo_post(i, h + 1) = q(1);  hi_post(i, h + 1) = q(2);

        if fmar
            % PRIMARY bands: FMAR quasi-Bayesian NW sandwich centred at
            % the block-adaptive posterior mean (same construction as
            % estimate_blp_fmar so coverage is comparable).
            beta_pm = out.beta_mean;                 % (m x 1)
            u = yh - Zh * beta_pm;
            u = u - mean(u);
            Gsc = Zh .* u;
            S = Gsc' * Gsc;
            for l = 1:min(Ltr, Nh - 1)
                Gl = Gsc(l+1:end, :)' * Gsc(1:end-l, :);
                S  = S + wNW(l) * (Gl + Gl');
            end
            V = ZtZinv * S * ZtZinv;
            sd_th = sqrt(max(b1n' * V(2:1 + K, 2:1 + K) * b1n, 0));
            lo(i, h + 1) = theta_mean(i, h + 1) - zcrit * sd_th;
            hi(i, h + 1) = theta_mean(i, h + 1) + zcrit * sd_th;
        else
            lo(i, h + 1) = lo_post(i, h + 1);
            hi(i, h + 1) = hi_post(i, h + 1);
        end

        tau_draws = sqrt(out.tau2_draws);            % (n_keep x G)
        tau_mean(i, :, h) = mean(tau_draws, 1);
        for g = 1:G
            tau_med(i, g, h) = empirical_quantile(tau_draws(:, g), 0.5);
        end
        lag1(i, h) = out.diag.lag1_acorr_beta;
        n_clip = n_clip + out.diag.n_tau_clip;
        if cfg.blp.return_draws
            draws{i, h} = th_draws;
        end
    end
end

blp.theta_mean  = theta_mean;
blp.theta_med   = theta_med;
blp.lo          = lo;
blp.hi          = hi;
blp.lo_post     = lo_post;
blp.hi_post     = hi_post;
blp.lambda      = lambda_used;
blp.prior_theta = prior_theta;
blp.tau_mean    = tau_mean;
blp.tau_med     = tau_med;
blp.diag.lag1_acorr = lag1;
blp.diag.n_tau_clip = n_clip;
if cfg.blp.return_draws
    blp.draws = draws;
end

assert(all(isfinite(theta_mean(:))), ...
    'estimate_blp_blockadaptive: non-finite estimates.');
end
