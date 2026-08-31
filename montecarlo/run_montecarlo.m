function mc = run_montecarlo(cfg, dgp_name)
% PURPOSE
% -------
% Monte Carlo experiment for ONE DGP: repeatedly simulate data, run all
% four estimators, and store IRFs, interval bounds, posterior block
% scales and diagnostics for later summary.
%
% MODEL / EQUATIONS
% -----------------
% Estimators compared (all normalised to a unit impact of the recursive
% shock on cfg.shock_var):
% cfg.mode = 'prototype' (original stack):
%   1 'LP'        ordinary local projection             (estimate_lp)
%   2 'VAR'       iterated VAR(p) IRF                   (estimate_var)
%   3 'BLP-glob'  global VAR-centred Bayesian LP        (estimate_blp_global)
%   4 'BLP-block' block-adaptive VAR-centred Bayesian LP(estimate_blp_blockadaptive)
% cfg.mode = 'fmar' (published FMAR machinery, see estimate_blp_fmar.m):
%   1 'LP'        ordinary local projection             (estimate_lp)
%   2 'BVAR'      GLP/FMAR Bayesian VAR                 (estimate_bvar_niw)
%   3 'BLP-FMAR'  FMAR global Bayesian LP               (estimate_blp_fmar)
%   4 'BLP-block' block-adaptive BLP nested in FMAR     (estimate_blp_blockadaptive)
% In both modes the block-adaptive estimator REUSES the lambdas selected
% by the global BLP so that the two differ only through the local scales.
%
% REPRODUCIBILITY
% ---------------
% Replication r under DGP d uses rng(cfg.seed + 100000*d + r, 'twister')
% set ONCE at the top of the replication; everything downstream (data
% simulation, CI simulation, both Gibbs runs) consumes that single
% stream, so any replication can be re-run in isolation.
%
% INPUTS
% ------
% cfg      : configuration struct (cfg.mc.* control the experiment;
%            cfg.mc.gibbs_n_burn/keep override the sampler size).
% dgp_name : 'correct' | 'sparse' | 'dense'.
%
% OUTPUTS
% -------
% mc : struct with fields
%   .dgp_name, .description, .misspec_block
%   .theta_true (K x (H+1))       true IRF (identical across reps:
%                                 DGP parameters are fixed)
%   .est_names  {1 x 4}
%   .theta      (4 x K x (H+1) x R) point estimates
%                                 (posterior means for the BLPs)
%   .lo, .hi    (4 x K x (H+1) x R) interval bounds
%   .tau_mean   (K x G x H x R)   posterior mean block scales (BLP-block)
%   .diag       struct arrays: VAR stability, sampler mixing, clip counts
%   .cfg        snapshot of the configuration actually used
%
% DIMENSIONS
% ----------
% R = cfg.mc.n_rep replications.
%
% NOTES
% -----
% Plain for-loop by default.  If cfg.mc.use_parfor is true AND parpool
% is available the loop can be switched manually; the prototype keeps
% the serial loop as the verified reference implementation.

% --- Resolve DGP -------------------------------------------------------
switch dgp_name
    case 'correct',  sim_fun = @simulate_var_dgp;            d_id = 1;
    case 'sparse',   sim_fun = @simulate_sparse_misspec_dgp; d_id = 2;
    case 'dense',    sim_fun = @simulate_dense_misspec_dgp;  d_id = 3;
    otherwise, error('run_montecarlo: unknown DGP "%s".', dgp_name);
end

if ~isfield(cfg, 'mode'), cfg.mode = 'prototype'; end

% Cheaper sampler inside the loop (documented in default_config.m).
cfg.gibbs.n_burn = cfg.mc.gibbs_n_burn;
cfg.gibbs.n_keep = cfg.mc.gibbs_n_keep;
cfg.blp.return_draws = false;

R = cfg.mc.n_rep;
K = cfg.K;  H = cfg.H;  G = K;          % scheme 'per_variable' => G = K

theta = zeros(4, K, H + 1, R);
lo    = nan(4, K, H + 1, R);
hi    = nan(4, K, H + 1, R);
tauM  = zeros(K, G, H, R);
diag_var_stable = false(R, 1);
diag_max_eig    = zeros(R, 1);
diag_lag1_glob  = zeros(R, 1);
diag_lag1_block = zeros(R, 1);
diag_n_clip     = zeros(R, 1);
diag_all_finite = false(R, 1);

fprintf('Monte Carlo: DGP = %s, R = %d replications\n', dgp_name, R);
tstart = tic;
theta_true = [];

for r = 1:R
    rng(cfg.seed + 100000 * d_id + r, 'twister');

    dgp = sim_fun(cfg);
    if r == 1
        theta_true = dgp.theta_true;    % fixed parameters => fixed truth
        mc.description   = dgp.description;
        mc.misspec_block = dgp.misspec_block;
    end

    if strcmp(cfg.mode, 'fmar')
        % FMAR stack: BVAR comparator + FMAR global BLP + block-adaptive
        % BLP nested in FMAR (same lambda_h, same detrending, same
        % sandwich bands; see estimate_blp_blockadaptive.m).
        var_est = estimate_bvar_niw(dgp.Y, cfg);
        lp_est  = estimate_lp(dgp.Y, cfg, var_est.b1n);
        blp_g   = estimate_blp_fmar(dgp.Y, cfg, var_est);
        blp_b   = estimate_blp_blockadaptive(dgp.Y, cfg, var_est, blp_g.lambda);
    else
        var_est = estimate_var(dgp.Y, cfg);
        lp_est  = estimate_lp(dgp.Y, cfg, var_est.b1n);
        blp_g   = estimate_blp_global(dgp.Y, cfg, var_est);
        blp_b   = estimate_blp_blockadaptive(dgp.Y, cfg, var_est, blp_g.lambda);
    end

    theta(1, :, :, r) = lp_est.theta;   lo(1, :, :, r) = lp_est.lo;   hi(1, :, :, r) = lp_est.hi;
    theta(2, :, :, r) = var_est.theta;  lo(2, :, :, r) = var_est.theta_lo; hi(2, :, :, r) = var_est.theta_hi;
    theta(3, :, :, r) = blp_g.theta_mean; lo(3, :, :, r) = blp_g.lo;  hi(3, :, :, r) = blp_g.hi;
    theta(4, :, :, r) = blp_b.theta_mean; lo(4, :, :, r) = blp_b.lo;  hi(4, :, :, r) = blp_b.hi;
    tauM(:, :, :, r)  = blp_b.tau_mean;

    diag_var_stable(r) = var_est.is_stable;
    diag_max_eig(r)    = var_est.max_eig;
    if isfield(blp_g.diag, 'lag1_acorr')
        diag_lag1_glob(r) = mean(blp_g.diag.lag1_acorr(:));
    else
        diag_lag1_glob(r) = NaN;    % FMAR global BLP is closed-form (no chain)
    end
    diag_lag1_block(r) = mean(blp_b.diag.lag1_acorr(:));
    diag_n_clip(r)     = blp_b.diag.n_tau_clip;
    th_r = theta(:, :, :, r);
    diag_all_finite(r) = all(isfinite(th_r(:)));

    if mod(r, 10) == 0 || r == R
        fprintf('  rep %3d / %3d  (elapsed %.1f s)\n', r, R, toc(tstart));
    end
end

mc.dgp_name   = dgp_name;
mc.theta_true = theta_true;
if strcmp(cfg.mode, 'fmar')
    mc.est_names = {'LP', 'BVAR', 'BLP-FMAR', 'BLP-block'};
else
    mc.est_names = {'LP', 'VAR', 'BLP-glob', 'BLP-block'};
end
mc.theta = theta;  mc.lo = lo;  mc.hi = hi;
mc.tau_mean = tauM;
mc.diag.var_stable = diag_var_stable;
mc.diag.max_eig    = diag_max_eig;
mc.diag.lag1_glob  = diag_lag1_glob;
mc.diag.lag1_block = diag_lag1_block;
mc.diag.n_tau_clip = diag_n_clip;
mc.diag.all_finite = diag_all_finite;
mc.cfg = cfg;

assert(all(diag_all_finite), 'run_montecarlo: non-finite estimates found.');
end
