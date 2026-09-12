function blp = estimate_blp_blockpooled(Y, cfg, bvar, lambda_mat)
% PURPOSE
% -------
% HORIZON-POOLED block-adaptive Bayesian Local Projection: the third
% estimator of the study, kept separate from
%   (1) estimate_blp_fmar.m            -- global FMAR baseline, tau = 1,
%   (2) estimate_blp_blockadaptive.m   -- tau_{i,g,h} INDEPENDENT by h.
%
% WHY IT EXISTS
% -------------
% In (2) every tau_{i,g,h} is identified by a single horizon-h
% regression, i.e. by p_g = p coefficient deviations from the VAR
% centre.  That is very little information: the resulting tau path is
% noisy, and the noise is passed straight into the posterior mean of
% beta_h, adding variance at horizons where the VAR prior is in fact
% fine.  Since the object being estimated -- "how badly does the VAR
% prior fail for block g?" -- varies smoothly with the horizon, pooling
% the tau path across h is the natural fix, and it nests both extremes:
%
%     log tau_{i,g,h} = phi * log tau_{i,g,h-1} + e,  e ~ N(0, kappa^2),
%
%   kappa -> 0    : one scale per block, fully pooled across horizons;
%   kappa -> inf  : back to the independent estimator (2);
%   kappa sampled : the data choose, per equation and block.
% With phi = 1 (default) this is a RANDOM WALK on log tau, i.e. a
% first-difference shrinkage prior; tau_{i,g,1} keeps the same
% half-Cauchy(0,1) anchor as (2), so the h = 1 marginal prior is
% unchanged.  All conditionals are documented in
% samplers/gibbs_block_pooled_horizons.m.
%
% EVERYTHING ELSE IS IDENTICAL TO (2)
% -----------------------------------
% Same detrending, same BVAR prior centre, same NW long-run scales
% psi(h), same lambda_h (passed in via lambda_mat so the three
% estimators differ ONLY through the treatment of tau), same conjugate
% sigma2-scaled prior, same FMAR quasi-Bayesian NW sandwich bands
% around the posterior mean.  Fixing tau = 1 therefore reproduces the
% FMAR closed form here too (tests/test_fmar_nesting_exact.m).
%
% INPUTS
% ------
% Y          : (T x K) raw data.
% cfg        : configuration struct; requires cfg.mode = 'fmar'.
%              Pooling settings live in cfg.blp.pool:
%                .phi        AR coefficient on log tau (default 1 = RW)
%                .a_kappa/.b_kappa  IG prior for kappa^2 (default 2, 0.1)
%                .kappa_fixed       optional: hold kappa there instead
%                                   of sampling it (sensitivity runs)
%                .mh_step    initial single-site Metropolis step (0.5)
%                .mh_step_level  initial level-move step (default 0.3)
%                .n_tau_sweeps   Metropolis sweeps of the tau block per
%                                Gibbs iteration (default 5)
% bvar       : output of estimate_bvar_niw.
% lambda_mat : OPTIONAL (K x H) tightness per equation/horizon (pass
%              blp_fmar.lambda to share the baseline's tightness).
%
% POINT ESTIMATE
% --------------
% As in estimate_blp_blockadaptive, cfg.blp.point_estimate selects
% 'rao_blackwell' (default: average the conditional means
% E[beta_h | tau, sigma2, y], which targets the same posterior mean with
% strictly less sampler noise) or 'draw_mean' (average the draws).  This
% matters more here than for the independent estimator: the pooled
% sampler's log-tau chain mixes slowly (see docs/APPROXIMATIONS.md), and
% Rao-Blackwellisation removes the part of that noise which would
% otherwise be charged to the estimator in the RMSE comparison.
%
% OUTPUTS
% -------
% blp : same interface as estimate_blp_blockadaptive, plus
%   .kappa_mean (K x G)  posterior mean smoothing sd of log tau
%   .diag.acc_rate (K x G x H)   Metropolis acceptance rates
%   .diag.ess_logtau (K x G x H) effective sample size of log tau
%
% DIMENSIONS
% ----------
% G = K blocks of p coefficients; m = 1 + K*p regressors; H horizons
% sampled JOINTLY per equation (one chain per equation, not per
% equation x horizon as in the independent estimator).
%
% NOTES
% -----
% Cost is comparable to (2): the same number of m x m Cholesky solves
% per sweep, plus G*H cheap scalar Metropolis steps.  RNG is controlled
% by the caller.

if nargin < 4, lambda_mat = []; end
if ~isfield(cfg, 'mode'), cfg.mode = 'prototype'; end
assert(strcmp(cfg.mode, 'fmar'), ...
    ['estimate_blp_blockpooled: implemented for cfg.mode = ''fmar'' only ' ...
     '(the pooling prior sits on top of the FMAR conjugate scaling).']);

[T, K] = size(Y);  %#ok<ASGLU>
p = cfg.p;  H = cfg.H;
alpha = 1 - cfg.ci_level;
zcrit = normal_quantile(1 - alpha / 2);
b1n = bvar.b1n;
m = 1 + K * p;

dt = var_deterministic_trend(Y, bvar.B, p);
x  = dt.x;
Tx = size(x, 1);
Zall = build_lp_regressors(x, p);

bp = build_block_prior(Y, cfg);
G  = bp.G;

pool = pool_options(cfg);
gopts.n_burn     = cfg.gibbs.n_burn;
gopts.n_keep     = cfg.gibbs.n_keep;
gopts.sample_tau = true;
gopts.phi        = pool.phi;
gopts.a_kappa    = pool.a_kappa;
gopts.b_kappa    = pool.b_kappa;
gopts.mh_step       = pool.mh_step;
gopts.mh_step_level = pool.mh_step_level;
gopts.n_tau_sweeps  = pool.n_tau_sweeps;
if ~isempty(pool.kappa_fixed), gopts.kappa_fixed = pool.kappa_fixed; end
gopts.x_min = 0.5 * log(cfg.gibbs.tau2_min);
gopts.x_max = 0.5 * log(cfg.gibbs.tau2_max);
if isfield(cfg.blp, 'fix_tau') && ~isempty(cfg.blp.fix_tau)
    gopts.sample_tau = false;
    gopts.tau_fixed  = cfg.blp.fix_tau;
end
% Blocks held at tau = 1 while the others adapt (cfg.blocks.fixed_tau;
% see default_config.m).  Empty = the original behaviour.
fixed_tau_blocks = [];
if isfield(cfg, 'blocks') && isfield(cfg.blocks, 'fixed_tau') && ~isempty(cfg.blocks.fixed_tau)
    fixed_tau_blocks = cfg.blocks.fixed_tau(:)';
    assert(all(fixed_tau_blocks >= 1 & fixed_tau_blocks <= G), ...
        'estimate_blp_blockpooled: cfg.blocks.fixed_tau must index blocks 1..%d.', G);
    gopts.fixed_blocks = fixed_tau_blocks;
end
if isfield(cfg.blp, 'tau_probs') && ~isempty(cfg.blp.tau_probs)
    tau_probs = cfg.blp.tau_probs;
else
    tau_probs = [0.05 0.25 0.50 0.75 0.95];
end
nQ = numel(tau_probs);

% Projection that turns a beta draw into the structural IRF value:
% theta_i(h) = b1n' * beta_h(2:1+K).
proj = zeros(m, 1);
proj(2:1 + K) = b1n;
gopts.proj = proj;

% --- horizon-level objects shared across equations ----------------------
% same optional prior-scale floor as estimate_blp_fmar (must be identical
% in the three estimators for the tau = 1 nesting)
psi_floor = isfield(cfg.fmar, 'psi_floor') && ~isempty(cfg.fmar.psi_floor) ...
            && logical(cfg.fmar.psi_floor);
J = [eye(K), zeros(K, K * (p - 1))];
Zh_all = cell(1, H);  Yh_all = cell(1, H);  Mu_all = cell(1, H);
psi_all = cell(1, H);  d_all = cell(1, H);  lam_all = zeros(K, H);
prior_theta = zeros(K, H + 1);  prior_theta(:, 1) = b1n;
for h = 1:H
    Zh_all{h} = Zall(1:end - h, :);
    Yh_all{h} = x(p + h:Tx, :);
    Fh = bvar.F^h;
    Mu_all{h} = [zeros(1, K); (J * Fh)'];
    prior_theta(:, h + 1) = Mu_all{h}(2:1 + K, :)' * b1n;
    psi_all{h} = fmar_prior_scale(x, p, h, psi_floor);
    d = zeros(m, 1);
    d(1) = cfg.fmar.Vc;
    for lag = 1:p
        d(1 + (lag-1)*K + (1:K)) = 1 ./ psi_all{h};
    end
    d_all{h} = d;
    if isempty(lambda_mat)
        lam_all(:, h) = select_lambda_fmar(Yh_all{h}, Zh_all{h}, Mu_all{h}, ...
                                           psi_all{h}, h, cfg);
    else
        lam_all(:, h) = lambda_mat(:, h);
    end
end

% --- containers ---------------------------------------------------------
theta_mean = zeros(K, H + 1);  theta_med = zeros(K, H + 1);
theta_cond = zeros(K, H + 1);
theta_rb = zeros(K, H + 1);  theta_dm = zeros(K, H + 1);
point_mode = 'rao_blackwell';
if isfield(cfg.blp, 'point_estimate') && ~isempty(cfg.blp.point_estimate)
    point_mode = cfg.blp.point_estimate;
end
assert(any(strcmp(point_mode, {'rao_blackwell', 'draw_mean'})), ...
    'estimate_blp_blockpooled: unknown cfg.blp.point_estimate ''%s''.', point_mode);
lo = nan(K, H + 1);  hi = nan(K, H + 1);
lo_post = nan(K, H + 1);  hi_post = nan(K, H + 1);
theta_mean(:, 1) = b1n;  theta_med(:, 1) = b1n;  theta_cond(:, 1) = b1n;
theta_rb(:, 1) = b1n;  theta_dm(:, 1) = b1n;
tau_mean = zeros(K, G, H);  tau_med = zeros(K, G, H);
tau_q = zeros(K, G, H, nQ);  p_tau_gt1 = zeros(K, G, H);
kappa_mean = zeros(K, G);
beta_all = zeros(m, K, H);              % full posterior-mean coefficient vectors
acc_rate = zeros(K, G, H);  acc_level = zeros(K, G);
ess_logtau = zeros(K, G, H);
lag1 = zeros(K, H);  ess_proj = zeros(K, H);  ess_rb = zeros(K, H);
post_sd = zeros(K, H);  mcse_rb = zeros(K, H);  n_clip = 0;

% --- one joint chain per equation ---------------------------------------
for i = 1:K
    hdata = struct('y', cell(1, H), 'Z', cell(1, H), 'mu', cell(1, H), ...
                   'd', cell(1, H), 'lambda', cell(1, H), ...
                   'a0', cell(1, H), 'b0', cell(1, H));
    for h = 1:H
        hdata(h).y      = Yh_all{h}(:, i);
        hdata(h).Z      = Zh_all{h};
        hdata(h).mu     = Mu_all{h}(:, i);
        hdata(h).d      = d_all{h};
        hdata(h).lambda = lam_all(i, h);
        hdata(h).a0     = (K + 2 - K + 1) / 2;      % = 3/2, as in (2)
        hdata(h).b0     = psi_all{h}(i) / 2;
    end

    prior.block_id = bp.block_id;
    out = gibbs_block_pooled_horizons(hdata, prior, gopts);

    for h = 1:H
        th = out.proj_draws(:, h);
        theta_dm(i, h + 1) = mean(th);
        theta_rb(i, h + 1) = proj' * out.beta_rb_mean(:, h);
        if strcmp(point_mode, 'draw_mean')
            theta_mean(i, h + 1) = theta_dm(i, h + 1);
        else
            theta_mean(i, h + 1) = theta_rb(i, h + 1);
        end
        theta_med(i, h + 1)  = empirical_quantile(th, 0.5);
        q = empirical_quantile(th, [alpha/2, 1 - alpha/2]);
        lo_post(i, h + 1) = q(1);  hi_post(i, h + 1) = q(2);
        theta_cond(i, h + 1) = proj' * out.beta_cond_mean(:, h);
        if strcmp(point_mode, 'draw_mean')
            beta_all(:, i, h) = out.beta_mean(:, h);
        else
            beta_all(:, i, h) = out.beta_rb_mean(:, h);
        end

        tdr = sqrt(squeeze(out.tau2_draws(:, :, h)));      % (n_keep x G)
        if G == 1, tdr = tdr(:); end
        tau_mean(i, :, h) = mean(tdr, 1);
        for g = 1:G
            tau_med(i, g, h)   = empirical_quantile(tdr(:, g), 0.5);
            tau_q(i, g, h, :)  = reshape(empirical_quantile(tdr(:, g), tau_probs), [1 1 1 nQ]);
            p_tau_gt1(i, g, h) = mean(tdr(:, g) > 1);
        end

        % PRIMARY bands: same FMAR NW sandwich as the other estimators.
        if strcmp(point_mode, 'draw_mean')
            beta_pm = out.beta_mean(:, h);
        else
            beta_pm = out.beta_rb_mean(:, h);
        end
        Zh = Zh_all{h};  yh = Yh_all{h}(:, i);
        Nh = size(Zh, 1);
        u = yh - Zh * beta_pm;  u = u - mean(u);
        Gsc = Zh .* u;
        S = Gsc' * Gsc;
        Ltr = h + 1;
        wNW = (Ltr + 1 - (1:Ltr)) / (Ltr + 1);
        for l = 1:min(Ltr, Nh - 1)
            Gl = Gsc(l+1:end, :)' * Gsc(1:end-l, :);
            S  = S + wNW(l) * (Gl + Gl');
        end
        ZtZinv = (Zh' * Zh) \ eye(m);
        V = ZtZinv * S * ZtZinv;
        sd_th = sqrt(max(b1n' * V(2:1 + K, 2:1 + K) * b1n, 0));
        lo(i, h + 1) = theta_mean(i, h + 1) - zcrit * sd_th;
        hi(i, h + 1) = theta_mean(i, h + 1) + zcrit * sd_th;
    end

    kappa_mean(i, :)   = mean(out.kappa_draws, 1);
    acc_rate(i, :, :)  = reshape(out.diag.acc_rate, [1, G, H]);
    acc_level(i, :)    = out.diag.acc_rate_level';
    ess_logtau(i, :, :) = reshape(out.diag.ess_logtau, [1, G, H]);
    lag1(i, :) = out.diag.lag1_acorr_beta;
    ess_proj(i, :) = out.diag.ess_proj;
    ess_rb(i, :)   = out.diag.ess_rb;
    for h = 1:H
        post_sd(i, h) = std(out.proj_draws(:, h));
        mcse_rb(i, h) = std(out.rb_draws(:, h)) / sqrt(out.diag.ess_rb(h));
    end
    n_clip = n_clip + out.diag.n_tau_clip;
end

blp.theta_mean  = theta_mean;
blp.theta_med   = theta_med;
blp.theta_cond  = theta_cond;
blp.theta_rb    = theta_rb;
blp.theta_draw_mean = theta_dm;
blp.point_estimate  = point_mode;
blp.lo          = lo;          blp.hi      = hi;
blp.lo_post     = lo_post;     blp.hi_post = hi_post;
blp.lambda      = lam_all;
blp.prior_theta = prior_theta;
blp.beta_mean   = beta_all;             % (m x K x H) reported posterior-mean coefficients
blp.tau_mean    = tau_mean;    blp.tau_med  = tau_med;
blp.tau_q       = tau_q;       blp.tau_probs = tau_probs(:)';
blp.p_tau_gt1   = p_tau_gt1;
blp.fixed_tau_blocks = fixed_tau_blocks;   % blocks held at tau = 1 ([] = none)
blp.kappa_mean  = kappa_mean;
blp.pool        = pool;
blp.diag.lag1_acorr = lag1;
blp.diag.acc_rate       = acc_rate;
blp.diag.acc_rate_level = acc_level;
blp.diag.ess_logtau = ess_logtau;
blp.diag.ess_proj   = ess_proj;   % ESS of the raw IRF draws
blp.diag.ess_rb     = ess_rb;     % ESS of the REPORTED (Rao-Blackwellised) mean
blp.diag.post_sd    = post_sd;   % posterior sd of the IRF at each (i, h)
blp.diag.mcse_rb    = mcse_rb;   % Monte Carlo se of the reported mean
blp.diag.mcse_ratio = mcse_rb ./ max(post_sd, realmin);
blp.diag.n_tau_clip = n_clip;

assert(all(isfinite(theta_mean(:))), ...
    'estimate_blp_blockpooled: non-finite estimates.');
end

% =====================================================================
function pool = pool_options(cfg)
% Resolve cfg.blp.pool with documented defaults, so configuration
% structs written before the pooled estimator existed remain valid.
pool = struct('phi', 1.0, 'a_kappa', 2.0, 'b_kappa', 0.1, ...
              'mh_step', 0.5, 'mh_step_level', 0.3, 'n_tau_sweeps', 5, ...
              'kappa_fixed', []);
if isfield(cfg, 'blp') && isstruct(cfg.blp) && isfield(cfg.blp, 'pool') ...
        && isstruct(cfg.blp.pool)
    f = fieldnames(pool);
    for k = 1:numel(f)
        if isfield(cfg.blp.pool, f{k}) && ~isempty(cfg.blp.pool.(f{k}))
            pool.(f{k}) = cfg.blp.pool.(f{k});
        end
    end
end
assert(pool.phi > 0 && pool.phi <= 1, ...
    'estimate_blp_blockpooled: cfg.blp.pool.phi must lie in (0, 1].');
end
