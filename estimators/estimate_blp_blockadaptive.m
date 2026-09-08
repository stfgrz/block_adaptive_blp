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
%         are estimated as LPs.  NOTE on h = 1: with the FMAR default
%         cfg.fmar.h1_mode = 'bvar' the global baseline reports the
%         BVAR itself at h = 1 while this estimator runs an h = 1 LP,
%         so the two are close but NOT identical there (and the lambda
%         inherited through lambda_mat(:,1) is then the BVAR's
%         Minnesota tightness, not an LP tightness).  Setting
%         cfg.fmar.h1_mode = 'lp' makes the baseline run an LP at h = 1
%         too, which restores exact tau = 1 nesting at EVERY h >= 1 and
%         is the harmonised setting for RMSE comparisons.  Either way
%         the preferred Monte Carlo metric integrates over h = 2..H.
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
%   .theta_mean/.theta_med (K x (H+1)) posterior mean/median IRFs.
%                          Which estimate of the posterior MEAN goes into
%                          .theta_mean is set by cfg.blp.point_estimate:
%                            'rao_blackwell' (DEFAULT) = average of the
%                              conditional means E[beta | tau, sigma2, y],
%                              which targets the same posterior mean with
%                              strictly less Monte Carlo error because the
%                              exact Gaussian draw noise is integrated out
%                              rather than simulated.  This matters for a
%                              FAIR RMSE comparison: the FMAR baseline is a
%                              closed form with NO sampler noise at all, so
%                              charging the adaptive estimators for avoidable
%                              simulation noise would understate them.
%                            'draw_mean' = plain average of the IRF draws
%                              (the original behaviour; keep it to
%                              reproduce results produced before this
%                              option existed).
%                          Both are always returned as .theta_rb and
%                          .theta_draw_mean.
%   .lo/.hi                (K x (H+1)) PRIMARY bands (posterior
%                          quantiles in 'prototype'; FMAR sandwich in
%                          'fmar')
%   .lo_post/.hi_post      (K x (H+1)) posterior-quantile bands (only
%                          in 'fmar'; equals .lo/.hi in 'prototype')
%   .lambda                (K x H) tightness actually used
%   .prior_theta           (K x (H+1)) IRF implied by the prior centre
%   .tau_mean/.tau_med     (K x G x H) posterior block-scale summaries
%   .tau_q                 (K x G x H x nQ) posterior quantiles of tau
%                          at the probabilities in .tau_probs
%   .p_tau_gt1             (K x G x H) posterior P(tau_{i,g,h} > 1)
%   .theta_cond            (K x (H+1)) IRF implied by the CONDITIONAL
%                          posterior mean at the last draw's (tau,
%                          sigma2).  With cfg.blp.fix_tau set this is
%                          the exact closed-form (Monte-Carlo-error-free)
%                          posterior mean -- used by the nesting tests.
%   .diag                  sampler diagnostics (.lag1_acorr, .ess_beta,
%                          .ess_logtau, .n_tau_clip)
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
if isfield(cfg.blp, 'tau_probs') && ~isempty(cfg.blp.tau_probs)
    gopts.tau_probs = cfg.blp.tau_probs;
else
    gopts.tau_probs = [0.05 0.25 0.50 0.75 0.95];
end
nQ = numel(gopts.tau_probs);
% Sensitivity switch for the cross-equation-covariance approximation:
% 'sample' (default) draws sigma2_i per equation as documented;
% 'fixed_niw' holds it at the diagonal of the FMAR SYSTEM posterior mode
% of Sigma, i.e. at a value that did use the cross-equation information.
% See montecarlo/run_sensitivity_approximations.m and
% docs/APPROXIMATIONS.md.
sigma2_mode = 'sample';
if isfield(cfg.blp, 'sigma2_mode') && ~isempty(cfg.blp.sigma2_mode)
    sigma2_mode = cfg.blp.sigma2_mode;
end
assert(any(strcmp(sigma2_mode, {'sample', 'fixed_niw'})), ...
    'estimate_blp_blockadaptive: unknown cfg.blp.sigma2_mode ''%s''.', sigma2_mode);
if strcmp(sigma2_mode, 'fixed_niw')
    assert(fmar, 'cfg.blp.sigma2_mode = ''fixed_niw'' needs cfg.mode = ''fmar''.');
end
% Nesting switch used by tests: force tau to a fixed value.
if isfield(cfg.blp, 'fix_tau') && ~isempty(cfg.blp.fix_tau)
    gopts.sample_tau = false;
    gopts.tau_fixed  = cfg.blp.fix_tau;
end

theta_mean = zeros(K, H + 1);  theta_med = zeros(K, H + 1);
theta_rb = zeros(K, H + 1);  theta_dm = zeros(K, H + 1);
point_mode = point_estimate_mode(cfg);
lo = nan(K, H + 1);  hi = nan(K, H + 1);
lo_post = nan(K, H + 1);  hi_post = nan(K, H + 1);
theta_mean(:, 1) = b1n;  theta_med(:, 1) = b1n;
theta_rb(:, 1) = b1n;  theta_dm(:, 1) = b1n;
prior_theta = zeros(K, H + 1);  prior_theta(:, 1) = b1n;
lambda_used = zeros(K, H);
tau_mean = zeros(K, G, H);  tau_med = zeros(K, G, H);
tau_q    = zeros(K, G, H, nQ);
p_tau_gt1 = zeros(K, G, H);
theta_cond = zeros(K, H + 1);  theta_cond(:, 1) = b1n;
lag1 = zeros(K, H);  n_clip = 0;
ess_beta = zeros(K, H);  ess_rb = zeros(K, H);  ess_logtau = zeros(K, G, H);
post_sd = zeros(K, H);  mcse_rb = zeros(K, H);
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
        if strcmp(sigma2_mode, 'fixed_niw')
            % System NIW posterior mode of Sigma at this horizon, using the
            % SAME lambda the estimator is about to use.
            if isempty(lambda_mat), lam_fix = lam_h; else, lam_fix = lambda_mat(1, h); end
            om = zeros(m, 1);  om(1) = cfg.fmar.Vc;
            for lag = 1:p
                om(1 + (lag-1)*K + (1:K)) = lam_fix^2 ./ psi_h;
            end
            [~, ~, Sig_sys] = niw_logml(Yh, Zh, Mu, om, psi_h, K + 2);
            sig2_fixed_h = diag(Sig_sys);
        end
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
            if strcmp(sigma2_mode, 'fixed_niw')
                gopts.sigma2_fixed = sig2_fixed_h(i);
            end
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
        theta_dm(i, h + 1) = mean(th_draws);
        theta_rb(i, h + 1) = out.beta_rb_mean(2:1 + K)' * b1n;
        if strcmp(point_mode, 'draw_mean')
            theta_mean(i, h + 1) = theta_dm(i, h + 1);
        else
            theta_mean(i, h + 1) = theta_rb(i, h + 1);
        end
        theta_med(i, h + 1)  = empirical_quantile(th_draws, 0.5);
        q = empirical_quantile(th_draws, [alpha/2, 1 - alpha/2]);
        lo_post(i, h + 1) = q(1);  hi_post(i, h + 1) = q(2);

        if fmar
            % PRIMARY bands: FMAR quasi-Bayesian NW sandwich centred at
            % the block-adaptive posterior mean (same construction as
            % estimate_blp_fmar so coverage is comparable).
            if strcmp(point_mode, 'draw_mean')
                beta_pm = out.beta_mean;             % (m x 1)
            else
                beta_pm = out.beta_rb_mean;          % same estimand, less noise
            end
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

        tau_mean(i, :, h)  = out.tau_mean';
        tau_med(i, :, h)   = out.tau_med';
        tau_q(i, :, h, :)  = reshape(out.tau_q, [1, G, 1, nQ]);
        p_tau_gt1(i, :, h) = out.p_tau_gt1';
        theta_cond(i, h + 1) = out.beta_cond_mean(2:1 + K)' * b1n;
        lag1(i, h) = out.diag.lag1_acorr_beta;
        ess_beta(i, h) = out.diag.ess_beta;
        rb_series      = out.beta_cond_draws(:, 2:1 + K) * b1n;
        ess_rb(i, h)   = ess_series(rb_series);
        post_sd(i, h)  = std(th_draws);
        mcse_rb(i, h)  = std(rb_series) / sqrt(ess_rb(i, h));
        ess_logtau(i, :, h) = out.diag.ess_logtau';
        n_clip = n_clip + out.diag.n_tau_clip;
        if cfg.blp.return_draws
            draws{i, h} = th_draws;
        end
    end
end

blp.theta_mean  = theta_mean;
blp.theta_med   = theta_med;
blp.theta_rb    = theta_rb;
blp.theta_draw_mean = theta_dm;
blp.point_estimate  = point_mode;
blp.sigma2_mode     = sigma2_mode;
blp.lo          = lo;
blp.hi          = hi;
blp.lo_post     = lo_post;
blp.hi_post     = hi_post;
blp.lambda      = lambda_used;
blp.prior_theta = prior_theta;
blp.tau_mean    = tau_mean;
blp.tau_med     = tau_med;
blp.tau_q       = tau_q;
blp.tau_probs   = gopts.tau_probs(:)';
blp.p_tau_gt1   = p_tau_gt1;
blp.theta_cond  = theta_cond;
blp.diag.lag1_acorr = lag1;
blp.diag.ess_beta   = ess_beta;   % ESS of the raw IRF draws
blp.diag.ess_rb     = ess_rb;     % ESS of the REPORTED (Rao-Blackwellised) mean
blp.diag.post_sd    = post_sd;   % posterior sd of the IRF at each (i, h)
blp.diag.mcse_rb    = mcse_rb;   % Monte Carlo se of the reported mean
blp.diag.mcse_ratio = mcse_rb ./ max(post_sd, realmin);
blp.diag.ess_logtau = ess_logtau;
blp.diag.n_tau_clip = n_clip;
if cfg.blp.return_draws
    blp.draws = draws;
end

assert(all(isfinite(theta_mean(:))), ...
    'estimate_blp_blockadaptive: non-finite estimates.');
end

% =====================================================================
function mode_str = point_estimate_mode(cfg)
% Resolve cfg.blp.point_estimate, defaulting to the Rao-Blackwellised
% posterior mean.  Configuration structs written before this option
% existed therefore switch to the lower-noise estimate of the SAME
% quantity; pass 'draw_mean' to reproduce the earlier numbers exactly.
mode_str = 'rao_blackwell';
if isfield(cfg, 'blp') && isstruct(cfg.blp) && ...
        isfield(cfg.blp, 'point_estimate') && ~isempty(cfg.blp.point_estimate)
    mode_str = cfg.blp.point_estimate;
end
assert(any(strcmp(mode_str, {'rao_blackwell', 'draw_mean'})), ...
    'estimate_blp_blockadaptive: unknown cfg.blp.point_estimate ''%s''.', mode_str);
end

function e = ess_series(v)
% Initial-positive-sequence effective sample size (Geyer 1992) of one
% scalar chain -- used for the mixing diagnostic of the reported IRF.
n = numel(v);
v = v - mean(v);
den = sum(v.^2);
if den <= 0, e = n; return; end
maxlag = min(n - 2, 200);
rho = zeros(1, maxlag);
for l = 1:maxlag
    rho(l) = sum(v(1:end-l) .* v(l+1:end)) / den;
end
ssum = 0;
for l = 1:2:maxlag - 1
    pair = rho(l) + rho(l + 1);
    if pair <= 0, break; end
    ssum = ssum + pair;
end
e = min(max(n / (1 + 2 * ssum), 1), n);
end
