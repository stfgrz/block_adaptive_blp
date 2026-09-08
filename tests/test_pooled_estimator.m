function test_pooled_estimator()
% PURPOSE
% -------
% Check the HORIZON-POOLED block-adaptive estimator
% (estimators/estimate_blp_blockpooled.m) does what its prior says it
% does, independently of whether it wins on RMSE.
%
% WHAT IS CHECKED
% ---------------
% 1. INTERFACE.  Same field names and shapes as the independent
%    block-adaptive estimator, so every downstream consumer
%    (run_montecarlo, tau_diagnostics, the plots) is interchangeable.
%
% 2. THE TWO LIMITS OF THE POOLING PRIOR.  With
%       log tau_{g,h} = phi log tau_{g,h-1} + e,  e ~ N(0, kappa^2):
%    * kappa -> 0 must give an (almost) FLAT tau path: the roughness
%      mean_h |log tau_h - log tau_{h-1}| must be far below
%    * kappa large, which must give a path at least as rough as the
%      INDEPENDENT estimator's.
%    These two runs bracket the sampled-kappa default and are what makes
%    "pooling" a claim about the estimator rather than a label.
%
% 3. SMOOTHING AT THE DEFAULT.  With kappa sampled, the tau path must be
%    smoother than the independent estimator's on the same data and the
%    same lambda -- this is the whole point of the estimator and the
%    mechanism behind any variance reduction it delivers.
%
% 4. SAMPLER HEALTH.  Metropolis acceptance rates in a usable range and
%    a reported effective sample size for the IRF draws.  (The log-tau
%    ESS is deliberately NOT asserted to be large: it is genuinely low,
%    which is documented in docs/APPROXIMATIONS.md; what must be healthy
%    is the quantity the estimator reports, the IRF.)
%
% The tau = 1 nesting of this estimator is checked in
% tests/test_fmar_nesting_exact.m.

fprintf('test_pooled_estimator:\n');

cfg = default_config();
cfg.mode = 'fmar';
cfg.H = 10;  cfg.T = 180;
cfg.gibbs.n_burn = 250;  cfg.gibbs.n_keep = 500;

rng(9001, 'twister');
dgp  = simulate_sparse_misspec_dgp(cfg);
bvar = estimate_bvar_niw(dgp.Y, cfg);
blpf = estimate_blp_fmar(dgp.Y, cfg, bvar);

rng(9002, 'twister');
ind = estimate_blp_blockadaptive(dgp.Y, cfg, bvar, blpf.lambda);
rng(9003, 'twister');
poo = estimate_blp_blockpooled(dgp.Y, cfg, bvar, blpf.lambda);

% --- 1. interface --------------------------------------------------------
need = {'theta_mean', 'theta_med', 'lo', 'hi', 'lo_post', 'hi_post', ...
        'lambda', 'prior_theta', 'tau_mean', 'tau_med', 'tau_q', ...
        'p_tau_gt1', 'theta_cond'};
for k = 1:numel(need)
    assert(isfield(poo, need{k}), 'pooled estimator lacks field .%s', need{k});
    assert(isequal(size(poo.(need{k})), size(ind.(need{k}))), ...
        'pooled .%s has size %s, independent has %s', need{k}, ...
        mat2str(size(poo.(need{k}))), mat2str(size(ind.(need{k}))));
end
fprintf('  interface matches the independent estimator: OK\n');

% --- 2. the two limits ---------------------------------------------------
cfg_tight = cfg;  cfg_tight.blp.pool.kappa_fixed = 1e-3;   % kappa -> 0
cfg_loose = cfg;  cfg_loose.blp.pool.kappa_fixed = 5;      % kappa -> large
rng(9004, 'twister');
p_tight = estimate_blp_blockpooled(dgp.Y, cfg_tight, bvar, blpf.lambda);
rng(9005, 'twister');
p_loose = estimate_blp_blockpooled(dgp.Y, cfg_loose, bvar, blpf.lambda);

r_tight = path_roughness(p_tight.tau_mean);
r_loose = path_roughness(p_loose.tau_mean);
r_pool  = path_roughness(poo.tau_mean);
r_ind   = path_roughness(ind.tau_mean);

assert(r_tight < 0.05, ...
    'kappa -> 0 did not flatten the tau path (roughness %.4f)', r_tight);
assert(r_loose > 5 * r_tight, ...
    'kappa large did not loosen the tau path (%.4f vs %.4f)', r_loose, r_tight);
fprintf(['  pooling limits: roughness %.4f (kappa=1e-3) < %.4f (sampled) ' ...
         '< %.4f (kappa=5)\n'], r_tight, r_pool, r_loose);

% --- 3. smoothing at the default ----------------------------------------
assert(r_pool < r_ind, ...
    ['the pooled tau path (%.4f) is not smoother than the independent ' ...
     'one (%.4f) -- the pooling prior is not doing its job'], r_pool, r_ind);
fprintf('  default pooling is smoother than independent: %.4f < %.4f  OK\n', ...
        r_pool, r_ind);

% --- 4. sampler health ---------------------------------------------------
acc  = poo.diag.acc_rate(:);
accL = poo.diag.acc_rate_level(:);
assert(all(acc > 0.05 & acc < 0.95), ...
    'single-site acceptance out of range [%.2f, %.2f]', min(acc), max(acc));
assert(all(accL > 0.02 & accL < 0.98), ...
    'level-move acceptance out of range [%.2f, %.2f]', min(accL), max(accL));
% The health check is on the MONTE CARLO ERROR of the quantity actually
% REPORTED, relative to the posterior spread it is summarising.  Raw ESS
% is not comparable across the two series -- Rao-Blackwellisation removes
% the independent draw noise, which LOWERS the ESS of what remains while
% still lowering the error of the average -- so the diagnostic to assert
% is mcse / posterior sd.
mcse_ratio = poo.diag.mcse_ratio(:);
assert(all(mcse_ratio < 0.25), ...
    ['the reported IRF mean carries too much sampler noise: ' ...
     'max mcse/posterior-sd = %.3f'], max(mcse_ratio));
fprintf(['  acceptance: site %.2f, level %.2f; reported mean mcse/post-sd ' ...
         'max %.3f (ESS of RB series min %.0f, raw draws %.0f, log-tau ' ...
         '%.0f -- see docs/APPROXIMATIONS.md)\n'], ...
        mean(acc), mean(accL), max(mcse_ratio), min(poo.diag.ess_rb(:)), ...
        min(poo.diag.ess_proj(:)), min(poo.diag.ess_logtau(:)));

% Rao-Blackwell must not be worse than averaging the draws.
mcse_draw = poo.diag.post_sd ./ sqrt(max(poo.diag.ess_proj, 1));
assert(mean(poo.diag.mcse_rb(:)) <= mean(mcse_draw(:)), ...
    ['Rao-Blackwellisation increased the Monte Carlo error ' ...
     '(%.4f vs %.4f) -- that should be impossible'], ...
    mean(poo.diag.mcse_rb(:)), mean(mcse_draw(:)));
fprintf('  Rao-Blackwell mcse %.4f <= draw-average mcse %.4f: OK\n', ...
        mean(poo.diag.mcse_rb(:)), mean(mcse_draw(:)));

fprintf('PASS: test_pooled_estimator\n\n');
end

% =====================================================================
function r = path_roughness(tau)
% Mean absolute first difference of log tau across horizons, averaged
% over equations and blocks.  0 = perfectly flat path.
lt = log(max(tau, realmin));            % (K x G x H)
d  = abs(diff(lt, 1, 3));
r  = mean(d(:));
end
