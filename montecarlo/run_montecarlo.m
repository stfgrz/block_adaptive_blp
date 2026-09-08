function mc = run_montecarlo(cfg, dgp_name)
% PURPOSE
% -------
% Monte Carlo experiment for ONE DGP: repeatedly simulate data, run the
% selected estimators, and store IRFs, interval bounds, posterior block
% scales and diagnostics for later summary.
%
% MODEL / EQUATIONS
% -----------------
% Estimators compared (all normalised to a unit impact of the recursive
% shock on cfg.shock_var).  Which ones are REPORTED is controlled by
% cfg.mc.estimators; the default reports all of them.
% cfg.mode = 'prototype' (original stack):
%   'lp'      'LP'         ordinary local projection      (estimate_lp)
%   'var'     'VAR'        iterated VAR(p) IRF            (estimate_var)
%   'global'  'BLP-glob'   global VAR-centred Bayesian LP (estimate_blp_global)
%   'block'   'BLP-block'  block-adaptive BLP             (estimate_blp_blockadaptive)
%   'pooled'  -- not available in prototype mode (skipped with a note)
% cfg.mode = 'fmar' (published FMAR machinery, see estimate_blp_fmar.m):
%   'lp'      'LP'         ordinary local projection      (estimate_lp)
%   'var'     'BVAR'       GLP/FMAR Bayesian VAR          (estimate_bvar_niw)
%   'global'  'BLP-FMAR'   FMAR global Bayesian LP        (estimate_blp_fmar)
%   'block'   'BLP-block'  block-adaptive, tau independent by horizon
%   'pooled'  'BLP-pooled' block-adaptive, log tau smoothed across
%                          horizons              (estimate_blp_blockpooled)
% In both modes the block-adaptive estimators REUSE the lambdas selected
% by the global BLP so that the three differ ONLY through the treatment
% of the local scales tau.
%
% REPRODUCIBILITY
% ---------------
% Replication r under DGP d uses rng(cfg.seed + 100000*d + r, 'twister')
% set ONCE at the top of the replication; everything downstream (data
% simulation, CI simulation, every sampler) consumes that single
% stream, so any replication can be re-run in isolation.
% The VAR/BVAR and the global BLP are ALWAYS computed (the block
% estimators need their prior centre and lambda), whether or not they
% are reported, so dropping estimators from cfg.mc.estimators never
% changes the random stream of the ones that remain -- with the single
% exception of 'pooled', which is drawn last.
% mc.meta records everything needed to re-run the experiment.
%
% CHUNKING (for running one experiment across several processes)
% --------------------------------------------------------------
% cfg.mc.rep_range = [r0 r1] runs ONLY replications r0..r1 of the same
% design.  Because replication r's seed is a pure function of r, a chunk
% is bit-identical to the corresponding slice of the full run, and
% montecarlo/merge_montecarlo.m glues chunks back together.  The stored
% struct records .rep_index so a partial result can never be mistaken
% for a complete one.
%
% INPUTS
% ------
% cfg      : configuration struct (cfg.mc.* control the experiment;
%            cfg.mc.gibbs_n_burn/keep override the sampler size).
% dgp_name : 'correct' | 'sparse' | 'intermediate' | 'dense'.
%
% OUTPUTS
% -------
% mc : struct with fields
%   .dgp_name, .description, .misspec_block, .dgp_params
%   .theta_true (K x (H+1))       true IRF (identical across reps:
%                                 DGP parameters are fixed)
%   .est_names  {1 x nE}
%   .theta      (nE x K x (H+1) x R) point estimates
%                                 (posterior means for the BLPs)
%   .lo, .hi    (nE x K x (H+1) x R) interval bounds
%   .tau_mean   (K x G x H x R)   posterior mean block scales (BLP-block;
%                                 kept at the top level for backward
%                                 compatibility with earlier results)
%   .tau        struct with per-estimator scale summaries:
%                 .block / .pooled, each with
%                 .mean/.med/.p_gt1 (K x G x H x R) and
%                 .q (K x G x H x nQ x R), the posterior quantiles of
%                 each tau at the probabilities in .probs
%   .lambda     (nE x K x H x R)  tightness actually used (NaN for LP/VAR)
%   .lo_post/.hi_post (nE x K x (H+1) x R) posterior-quantile bands for
%               the sampled estimators (NaN elsewhere).  The PRIMARY
%               bands in .lo/.hi are FMAR's quasi-Bayesian NW sandwich;
%               keeping both lets the sensitivity analysis report what
%               that choice costs in coverage.
%   .diag       struct arrays: VAR stability, sampler mixing, ESS,
%               Metropolis acceptance, clip counts
%   .meta       standardised reproducibility record (see mc_meta below)
%   .cfg        snapshot of the configuration actually used
%
% DIMENSIONS
% ----------
% R = cfg.mc.n_rep replications, nE = number of reported estimators.
%
% NOTES
% -----
% Plain for-loop by default.  If cfg.mc.use_parfor is true AND parpool
% is available the loop can be switched manually; the prototype keeps
% the serial loop as the verified reference implementation.

% --- Resolve DGP -------------------------------------------------------
switch dgp_name
    case 'correct',      sim_fun = @simulate_var_dgp;              d_id = 1;
    case 'sparse',       sim_fun = @simulate_sparse_misspec_dgp;   d_id = 2;
    case 'dense',        sim_fun = @simulate_dense_misspec_dgp;    d_id = 3;
    case 'intermediate', sim_fun = @simulate_intermediate_misspec_dgp; d_id = 4;
    otherwise, error('run_montecarlo: unknown DGP "%s".', dgp_name);
end

if ~isfield(cfg, 'mode'), cfg.mode = 'prototype'; end
fmar = strcmp(cfg.mode, 'fmar');

% --- Which estimators are reported -------------------------------------
if isfield(cfg.mc, 'estimators') && ~isempty(cfg.mc.estimators)
    est_keys = cfg.mc.estimators(:)';
else
    est_keys = {'lp', 'var', 'global', 'block'};
end
if fmar
    name_of = struct('lp', 'LP', 'var', 'BVAR', 'global', 'BLP-FMAR', ...
                     'block', 'BLP-block', 'pooled', 'BLP-pooled');
else
    name_of = struct('lp', 'LP', 'var', 'VAR', 'global', 'BLP-glob', ...
                     'block', 'BLP-block', 'pooled', 'BLP-pooled');
    if any(strcmp(est_keys, 'pooled'))
        fprintf(['run_montecarlo: ''pooled'' is a FMAR-mode estimator; ' ...
                 'skipped in prototype mode.\n']);
        est_keys = est_keys(~strcmp(est_keys, 'pooled'));
    end
end
for k = 1:numel(est_keys)
    assert(isfield(name_of, est_keys{k}), ...
        'run_montecarlo: unknown estimator key "%s".', est_keys{k});
end
est_names = cell(1, numel(est_keys));
for k = 1:numel(est_keys), est_names{k} = name_of.(est_keys{k}); end
nE = numel(est_keys);
want_pooled = any(strcmp(est_keys, 'pooled'));

% Cheaper sampler inside the loop (documented in default_config.m).
cfg.gibbs.n_burn = cfg.mc.gibbs_n_burn;
cfg.gibbs.n_keep = cfg.mc.gibbs_n_keep;
cfg.blp.return_draws = false;

R = cfg.mc.n_rep;
if isfield(cfg.mc, 'rep_range') && ~isempty(cfg.mc.rep_range)
    rep_index = cfg.mc.rep_range(1):cfg.mc.rep_range(2);
    assert(all(rep_index >= 1) && all(rep_index <= R), ...
        'run_montecarlo: cfg.mc.rep_range must lie inside 1..cfg.mc.n_rep.');
else
    rep_index = 1:R;
end
Rn = numel(rep_index);
K = cfg.K;  H = cfg.H;  G = K;          % scheme 'per_variable' => G = K

R = Rn;                                 % array sizes follow the chunk
theta = zeros(nE, K, H + 1, R);
lo    = nan(nE, K, H + 1, R);
hi    = nan(nE, K, H + 1, R);
lambda_store = nan(nE, K, H, R);
% Secondary (posterior-quantile) bands for the sampled estimators: the
% PRIMARY bands are FMAR's quasi-Bayesian NW sandwich, and comparing the
% two coverages is how the sensitivity analysis prices that choice
% (docs/APPROXIMATIONS.md, item 3).
lo_post = nan(nE, K, H + 1, R);
hi_post = nan(nE, K, H + 1, R);
tauB = zeros(K, G, H, R);  tauB_med = zeros(K, G, H, R);  tauB_p1 = zeros(K, G, H, R);
% Posterior QUANTILES of every tau_{i,g,h}, averaged over replications at
% the end: the mean and median alone do not say whether a large scale is
% a shifted posterior or a heavy right tail.
tauB_q = [];  tauP_q = [];
if want_pooled
    tauP = zeros(K, G, H, R);  tauP_med = zeros(K, G, H, R);  tauP_p1 = zeros(K, G, H, R);
    kappaP = zeros(K, G, R);
    accP = zeros(R, 1);  essP = zeros(R, 1);
end
diag_var_stable = false(R, 1);
diag_max_eig    = zeros(R, 1);
diag_lag1_glob  = zeros(R, 1);
diag_lag1_block = zeros(R, 1);
diag_n_clip     = zeros(R, 1);
diag_ess_block  = zeros(R, 1);
diag_ess_tau    = zeros(R, 1);
diag_all_finite = false(R, 1);
seeds           = zeros(R, 1);

if Rn == cfg.mc.n_rep
    fprintf('Monte Carlo: DGP = %s, mode = %s, R = %d, estimators = %s\n', ...
            dgp_name, cfg.mode, Rn, strjoin(est_names, ', '));
else
    fprintf('Monte Carlo: DGP = %s, mode = %s, reps %d-%d of %d, estimators = %s\n', ...
            dgp_name, cfg.mode, rep_index(1), rep_index(end), cfg.mc.n_rep, ...
            strjoin(est_names, ', '));
end
tstart = tic;
theta_true = [];

for r = 1:Rn
    seeds(r) = cfg.seed + 100000 * d_id + rep_index(r);
    rng(seeds(r), 'twister');

    dgp = sim_fun(cfg);
    if r == 1
        theta_true = dgp.theta_true;    % fixed parameters => fixed truth
        mc.description   = dgp.description;
        mc.misspec_block = dgp.misspec_block;
        if isfield(dgp, 'params'), mc.dgp_params = dgp.params;
        else,                      mc.dgp_params = struct(); end
    end

    if fmar
        var_est = estimate_bvar_niw(dgp.Y, cfg);
        lp_est  = estimate_lp(dgp.Y, cfg, var_est.b1n);
        blp_g   = estimate_blp_fmar(dgp.Y, cfg, var_est);
        blp_b   = estimate_blp_blockadaptive(dgp.Y, cfg, var_est, blp_g.lambda);
        if want_pooled
            blp_p = estimate_blp_blockpooled(dgp.Y, cfg, var_est, blp_g.lambda);
        end
    else
        var_est = estimate_var(dgp.Y, cfg);
        lp_est  = estimate_lp(dgp.Y, cfg, var_est.b1n);
        blp_g   = estimate_blp_global(dgp.Y, cfg, var_est);
        blp_b   = estimate_blp_blockadaptive(dgp.Y, cfg, var_est, blp_g.lambda);
    end

    for k = 1:nE
        switch est_keys{k}
            case 'lp'
                theta(k, :, :, r) = lp_est.theta;
                lo(k, :, :, r) = lp_est.lo;   hi(k, :, :, r) = lp_est.hi;
            case 'var'
                theta(k, :, :, r) = var_est.theta;
                lo(k, :, :, r) = var_est.theta_lo;  hi(k, :, :, r) = var_est.theta_hi;
            case 'global'
                theta(k, :, :, r) = blp_g.theta_mean;
                lo(k, :, :, r) = blp_g.lo;   hi(k, :, :, r) = blp_g.hi;
                lambda_store(k, :, :, r) = blp_g.lambda;
            case 'block'
                theta(k, :, :, r) = blp_b.theta_mean;
                lo(k, :, :, r) = blp_b.lo;   hi(k, :, :, r) = blp_b.hi;
                lo_post(k, :, :, r) = blp_b.lo_post;
                hi_post(k, :, :, r) = blp_b.hi_post;
                lambda_store(k, :, :, r) = blp_b.lambda;
            case 'pooled'
                theta(k, :, :, r) = blp_p.theta_mean;
                lo(k, :, :, r) = blp_p.lo;   hi(k, :, :, r) = blp_p.hi;
                lo_post(k, :, :, r) = blp_p.lo_post;
                hi_post(k, :, :, r) = blp_p.hi_post;
                lambda_store(k, :, :, r) = blp_p.lambda;
        end
    end

    tauB(:, :, :, r)     = blp_b.tau_mean;
    tauB_med(:, :, :, r) = blp_b.tau_med;
    if isfield(blp_b, 'p_tau_gt1'), tauB_p1(:, :, :, r) = blp_b.p_tau_gt1; end
    if isfield(blp_b, 'tau_q')
        if isempty(tauB_q), tauB_q = zeros([size(blp_b.tau_q), R]); end
        tauB_q(:, :, :, :, r) = blp_b.tau_q;
    end
    if want_pooled
        tauP(:, :, :, r)     = blp_p.tau_mean;
        tauP_med(:, :, :, r) = blp_p.tau_med;
        tauP_p1(:, :, :, r)  = blp_p.p_tau_gt1;
        if isempty(tauP_q), tauP_q = zeros([size(blp_p.tau_q), R]); end
        tauP_q(:, :, :, :, r) = blp_p.tau_q;
        kappaP(:, :, r)      = blp_p.kappa_mean;
        accP(r) = mean(blp_p.diag.acc_rate(:));
        essP(r) = mean(blp_p.diag.ess_logtau(:));
    end

    diag_var_stable(r) = var_est.is_stable;
    diag_max_eig(r)    = var_est.max_eig;
    if isfield(blp_g.diag, 'lag1_acorr')
        diag_lag1_glob(r) = mean(blp_g.diag.lag1_acorr(:));
    else
        diag_lag1_glob(r) = NaN;    % FMAR global BLP is closed-form (no chain)
    end
    diag_lag1_block(r) = mean(blp_b.diag.lag1_acorr(:));
    diag_n_clip(r)     = blp_b.diag.n_tau_clip;
    if isfield(blp_b.diag, 'ess_beta'),   diag_ess_block(r) = mean(blp_b.diag.ess_beta(:));   else, diag_ess_block(r) = NaN; end
    if isfield(blp_b.diag, 'ess_logtau'), diag_ess_tau(r)   = mean(blp_b.diag.ess_logtau(:)); else, diag_ess_tau(r)   = NaN; end
    th_r = theta(:, :, :, r);
    diag_all_finite(r) = all(isfinite(th_r(:)));

    if mod(r, 10) == 0 || r == Rn
        fprintf('  rep %3d / %3d  (elapsed %.1f s)\n', r, Rn, toc(tstart));
    end
end

mc.dgp_name   = dgp_name;
mc.theta_true = theta_true;
mc.est_names  = est_names;
mc.est_keys   = est_keys;
mc.theta = theta;  mc.lo = lo;  mc.hi = hi;
mc.lambda = lambda_store;
mc.lo_post = lo_post;  mc.hi_post = hi_post;
mc.tau_mean = tauB;                       % backward-compatible top level
mc.tau.block.mean = tauB;
mc.tau.block.med  = tauB_med;
mc.tau.block.p_gt1 = tauB_p1;
if ~isempty(tauB_q)
    % Stored PER REPLICATION (K x G x H x nQ x R), like every other tau
    % array, so that chunks of a parallel run merge by concatenation and
    % a merged result is identical to a serial one.  About 3.6 MB at
    % R = 500; tau_diagnostics averages over the last dimension.
    mc.tau.block.q      = tauB_q;
    mc.tau.block.probs  = blp_b.tau_probs;
end
if want_pooled
    mc.tau.pooled.mean  = tauP;
    mc.tau.pooled.med   = tauP_med;
    mc.tau.pooled.p_gt1 = tauP_p1;
    mc.tau.pooled.kappa = kappaP;
    if ~isempty(tauP_q)
        mc.tau.pooled.q     = tauP_q;
        mc.tau.pooled.probs = blp_p.tau_probs;
    end
    mc.diag.pooled_acc_rate  = accP;
    mc.diag.pooled_ess_logtau = essP;
end
mc.diag.var_stable = diag_var_stable;
mc.diag.max_eig    = diag_max_eig;
mc.diag.lag1_glob  = diag_lag1_glob;
mc.diag.lag1_block = diag_lag1_block;
mc.diag.n_tau_clip = diag_n_clip;
mc.diag.ess_beta_block  = diag_ess_block;
mc.diag.ess_logtau_block = diag_ess_tau;
mc.diag.all_finite = diag_all_finite;
mc.seeds = seeds;
mc.rep_index = rep_index;
mc.meta  = mc_meta(cfg, dgp_name, d_id, est_keys, est_names, mc.dgp_params, ...
                   mc.misspec_block, seeds, toc(tstart));
mc.meta.rep_index = rep_index;
mc.meta.R_done    = Rn;
mc.meta.is_complete = (Rn == cfg.mc.n_rep);
% Store a SAVE-SAFE snapshot: cfg.fmar.hyper_sd_rule is an anonymous
% function, which `save` cannot write (Octave) or writes lossily
% (MATLAB).  utils/cfg_to_savable.m turns handles into their source
% text; utils/cfg_from_saved.m converts them back on load.
mc.cfg = cfg_to_savable(cfg);

assert(all(diag_all_finite), 'run_montecarlo: non-finite estimates found.');
end

% =====================================================================
function meta = mc_meta(cfg, dgp_name, d_id, est_keys, est_names, ...
                        dgp_params, misspec_block, seeds, elapsed)
% Standardised reproducibility record attached to every Monte Carlo
% result.  Everything needed to re-run the experiment from a clean
% session lives here, so no downstream script has to guess at the
% settings behind a stored .mat file.
meta.dgp_name       = dgp_name;
meta.dgp_id         = d_id;
meta.dgp_params     = dgp_params;
meta.misspec_block  = misspec_block;
meta.mode           = cfg.mode;
meta.R              = cfg.mc.n_rep;
meta.T              = cfg.T;
meta.p              = cfg.p;
meta.K              = cfg.K;
meta.H              = cfg.H;
meta.burnin         = cfg.burnin;
meta.shock_var      = cfg.shock_var;
meta.ci_level       = cfg.ci_level;
meta.master_seed    = cfg.seed;
meta.seed_rule      = 'rng(cfg.seed + 100000*dgp_id + r, ''twister'')';
meta.seeds          = seeds(:)';
meta.est_keys       = est_keys;
meta.est_names      = est_names;
meta.blocks_scheme  = cfg.blocks.scheme;
meta.gibbs_n_burn   = cfg.gibbs.n_burn;
meta.gibbs_n_keep   = cfg.gibbs.n_keep;
meta.lp_vcov        = cfg.lp.vcov;
if strcmp(cfg.mode, 'fmar')
    meta.fmar = cfg.fmar;
    if ~isfield(meta.fmar, 'h1_mode'), meta.fmar.h1_mode = 'bvar'; end
    % A function handle does not survive `save` cleanly in every
    % version; store the rule's values at the horizons actually used.
    meta.fmar.hyper_sd_values = cfg.fmar.hyper_sd_rule(1:cfg.H);
    meta.fmar = rmfield(meta.fmar, 'hyper_sd_rule');
else
    meta.lambda_mode = cfg.blp.lambda_mode;
    meta.lambda_grid = cfg.blp.lambda_grid;
end
if isfield(cfg.blp, 'pool'), meta.pool = cfg.blp.pool; end
meta.elapsed_seconds = elapsed;
meta.created         = datestr(now, 'yyyy-mm-dd HH:MM:SS');  %#ok<TNOW1,DATST>
end
