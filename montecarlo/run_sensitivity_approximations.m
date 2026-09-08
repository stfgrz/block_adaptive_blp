function S = run_sensitivity_approximations(opts)
% PURPOSE
% -------
% Quantify what the three remaining METHODOLOGICAL APPROXIMATIONS of the
% block-adaptive estimator actually cost.  Each is documented in
% docs/APPROXIMATIONS.md; this function supplies the numbers that
% document sits on, so the claim "this approximation is harmless here"
% is an empirical statement about this design rather than an assertion.
%
% THE THREE APPROXIMATIONS
% ------------------------
% A. CROSS-EQUATION COVARIANCE IN THE ADAPTIVE SAMPLER.
%    The sampler runs equation by equation with a scalar sigma2_i and no
%    off-diagonal Sigma.  At tau = 1 this provably cannot matter for the
%    posterior MEAN of beta -- sigma2 cancels from conditional (1'), and
%    tests/test_fmar_nesting_exact.m verifies the resulting equality to
%    3e-16.  With tau SAMPLED it can matter indirectly, because the tau
%    conditional (3') divides by sigma2_i.  The test: re-run with
%    sigma2_i HELD at the diagonal of the system NIW posterior mode
%    (cfg.blp.sigma2_mode = 'fixed_niw'), which is a value that did use
%    the cross-equation information, and measure how far the IRF and the
%    tau surface move -- both in IRF units and relative to the posterior
%    standard deviation of the IRF itself.
%
% B. FIXED PRIOR CENTRE, FIXED lambda_h, FIXED IMPACT VECTOR.
%    The h >= 2 estimators condition on the BVAR posterior MEAN for the
%    prior centre and on a single lambda_h, and hold the impact vector
%    b1n at its point estimate (all three exactly as FMAR do).  Three
%    separate checks:
%      B1 IMPACT VECTOR.  theta_i(h) = b1n' beta_h(2:1+K) is LINEAR in
%         b1n, so the contribution of b1n uncertainty can be isolated
%         exactly: hold the coefficient posterior fixed and re-evaluate
%         theta at each NIW posterior draw of b1n.  Report the resulting
%         standard deviation against the reported band half-width.  This
%         is the cleanest of the three because no re-estimation is
%         needed and no approximation is introduced by the check itself.
%      B2 PRIOR CENTRE.  Re-run the whole block-adaptive estimator with
%         the prior centre (and b1n) taken from BVAR posterior DRAWS
%         instead of the posterior mean, for a handful of draws, and
%         report the spread of the resulting IRFs.
%      B3 lambda.  Compare the estimator when it INHERITS the global
%         BLP's lambda_h (the Monte Carlo default, which isolates the
%         effect of tau) against re-selecting lambda_h by marginal
%         likelihood inside the adaptive estimator.
%
% THE YARDSTICK (why a raw "the IRF moved by 0.01" means nothing)
% ---------------------------------------------------------------
% Every comparison here re-runs a sampler, so part of any difference is
% simply that a Markov chain was run twice.  Two things are therefore
% done.  First, each variant is run from the SAME seed as the baseline,
% so the streams are paired and only the modelling change differs.
% Second, a NOISE BASELINE is computed: the identical estimator run from
% a DIFFERENT seed.  A modelling change that moves the IRF by less than
% that baseline has not been shown to matter at all.  Every reported
% change is quoted next to it.
%
% C. QUASI-BAYESIAN HAC BANDS.
%    The primary intervals are FMAR's Newey-West sandwich around the
%    posterior mean, not posterior quantiles.  Both are computed and
%    stored by the Monte Carlo, so if a stored result is supplied
%    (opts.mc_file) this reports the two coverages and the two average
%    widths side by side, per horizon.  Without a stored result this
%    part is skipped with a message -- it is a repeated-sampling
%    statement and cannot be made on a single dataset.
%
% INPUTS (opts, all optional)
% ---------------------------
%   opts.n_data   number of simulated datasets for A/B (default 12)
%   opts.dgp      DGP name (default 'sparse')
%   opts.preset   mc_preset name for the configuration (default 'grid')
%   opts.n_draws  BVAR posterior draws used in B2 (default 20)
%   opts.mc_file  a saved Monte Carlo .mat (with `mc`) for part C
%   opts.out_csv  where to write the summary table
%   opts.seed     master seed (default 202)
%
% OUTPUTS
% -------
% S : struct with one field per part (.A, .B1, .B2, .B3, .C), each
%     carrying the numbers printed and written to the csv.
%
% NOTES
% -----
% Deterministic given opts.seed.  Runs at the exploratory preset by
% default: these are magnitude comparisons, not precision estimates.

if nargin < 1, opts = struct(); end
opts = dflt(opts, 'n_data', 12);
opts = dflt(opts, 'dgp', 'sparse');
opts = dflt(opts, 'preset', 'grid');
opts = dflt(opts, 'n_draws', 20);
opts = dflt(opts, 'seed', 202);
opts = dflt(opts, 'mc_file', '');
root = fileparts(fileparts(mfilename('fullpath')));
opts = dflt(opts, 'out_csv', fullfile(root, 'results', 'sensitivity_approximations.csv'));

cfg = mc_preset(opts.preset);
cfg.gibbs.n_burn = cfg.mc.gibbs_n_burn;
cfg.gibbs.n_keep = cfg.mc.gibbs_n_keep;
switch opts.dgp
    case 'correct',      sim = @simulate_var_dgp;
    case 'sparse',       sim = @simulate_sparse_misspec_dgp;
    case 'intermediate', sim = @simulate_intermediate_misspec_dgp;
    case 'dense',        sim = @simulate_dense_misspec_dgp;
    otherwise, error('run_sensitivity_approximations: unknown DGP %s', opts.dgp);
end

K = cfg.K;  H = cfg.H;
nD = opts.n_data;
dA_irf = zeros(nD, 1);  dA_rel = zeros(nD, 1);  dA_tau = zeros(nD, 1);
sd_b1n = zeros(nD, 1);  half_band = zeros(nD, 1);
sd_centre = zeros(nD, 1);
dB3_irf = zeros(nD, 1);  dB3_rel = zeros(nD, 1);
lam_inherit = zeros(nD, 1);  lam_reselect = zeros(nD, 1);
dN_irf = zeros(nD, 1);  dN_rel = zeros(nD, 1);  dN_tau = zeros(nD, 1);

fprintf('\n=== sensitivity of the three approximations (%s DGP, %d datasets) ===\n', ...
        opts.dgp, nD);
for r = 1:nD
    rng(opts.seed + r, 'twister');
    dgp = sim(cfg);
    Y = dgp.Y;

    bvar  = estimate_bvar_niw(Y, cfg);
    blp_f = estimate_blp_fmar(Y, cfg, bvar);

    rng(opts.seed + 1000 + r, 'twister');
    base = estimate_blp_blockadaptive(Y, cfg, bvar, blp_f.lambda);
    post_sd = base.diag.post_sd;                    % (K x H)

    % ---- A. cross-equation covariance --------------------------------
    cfgA = cfg;  cfgA.blp.sigma2_mode = 'fixed_niw';
    rng(opts.seed + 1000 + r, 'twister');           % same stream start
    altA = estimate_blp_blockadaptive(Y, cfgA, bvar, blp_f.lambda);
    dirf = abs(altA.theta_mean(:, 3:end) - base.theta_mean(:, 3:end));
    dA_irf(r) = max(dirf(:));
    dA_rel(r) = max(max(dirf ./ max(post_sd(:, 2:end), realmin)));
    dtau = abs(log(altA.tau_mean) - log(base.tau_mean));
    dA_tau(r) = max(dtau(:));

    % ---- B1. impact-vector uncertainty (EXACT, no re-estimation) -----
    % theta_i(h) = b1n' * beta_block(i,:,h)' is LINEAR in b1n, so
    % re-evaluating the SAME coefficient posterior mean at each NIW draw
    % of b1n isolates that channel with no approximation whatsoever.
    nd = size(bvar.b1n_draws, 2);
    th_b1 = zeros(K, H, nd);
    for j = 1:nd
        b1d = bvar.b1n_draws(:, j);
        for h = 1:H
            th_b1(:, h, j) = base.beta_block(:, :, h) * b1d;
        end
    end
    sd_b1n(r) = mean(mean(std(th_b1(:, 2:end, :), 0, 3)));
    half_band(r) = mean(mean((base.hi(:, 3:end) - base.lo(:, 3:end)) / 2));

    % ---- B2. prior-centre uncertainty --------------------------------
    % Only the prior CENTRE and the impact vector are varied here; the
    % deterministic trend is left at the posterior-mean B, because the
    % draws expose the companion matrix but not the constant, and
    % detrending uncertainty is a separate (smaller) channel.
    nd2 = min(opts.n_draws, size(bvar.F_draws, 3));
    th_c = zeros(K, H + 1, nd2);
    for j = 1:nd2
        bd = bvar;
        bd.F   = bvar.F_draws(:, :, j);
        bd.b1n = bvar.b1n_draws(:, j);
        rng(opts.seed + 2000 + 37 * r + j, 'twister');
        e = estimate_blp_blockadaptive(Y, cfg, bd, blp_f.lambda);
        th_c(:, :, j) = e.theta_mean;
    end
    sd_centre(r) = mean(mean(std(th_c(:, 3:end, :), 0, 3)));

    % ---- NOISE BASELINE: the same estimator, a different seed --------
    % Anything smaller than this has not been shown to be a consequence
    % of a modelling choice rather than of running a chain twice.
    rng(opts.seed + 5000 + r, 'twister');
    altN = estimate_blp_blockadaptive(Y, cfg, bvar, blp_f.lambda);
    dn = abs(altN.theta_mean(:, 3:end) - base.theta_mean(:, 3:end));
    dN_irf(r) = max(dn(:));
    dN_rel(r) = max(max(dn ./ max(post_sd(:, 2:end), realmin)));
    dnt = abs(log(altN.tau_mean) - log(base.tau_mean));
    dN_tau(r) = max(dnt(:));

    % ---- B3. inherited vs re-selected lambda -------------------------
    % SAME seed as the baseline, so the only difference is where
    % lambda_h came from.
    rng(opts.seed + 1000 + r, 'twister');
    altB3 = estimate_blp_blockadaptive(Y, cfg, bvar, []);   % re-select
    d3 = abs(altB3.theta_mean(:, 3:end) - base.theta_mean(:, 3:end));
    dB3_irf(r) = max(d3(:));
    dB3_rel(r) = max(max(d3 ./ max(post_sd(:, 2:end), realmin)));
    lam_inherit(r)  = mean(mean(base.lambda(:, 2:end)));
    lam_reselect(r) = mean(mean(altB3.lambda(:, 2:end)));

    if mod(r, 4) == 0 || r == nD
        fprintf('  dataset %d/%d done\n', r, nD);
    end
end

S.A  = struct('max_dIRF', mean(dA_irf), 'max_dIRF_over_postsd', mean(dA_rel), ...
              'max_dlogtau', mean(dA_tau));
S.B1 = struct('sd_from_b1n', mean(sd_b1n), 'reported_half_band', mean(half_band), ...
              'ratio', mean(sd_b1n) / mean(half_band));
S.B2 = struct('sd_from_prior_centre', mean(sd_centre), ...
              'reported_half_band', mean(half_band), ...
              'ratio', mean(sd_centre) / mean(half_band));
S.B3 = struct('max_dIRF', mean(dB3_irf), 'max_dIRF_over_postsd', mean(dB3_rel), ...
              'lambda_inherited', mean(lam_inherit), ...
              'lambda_reselected', mean(lam_reselect), ...
              'max_lambda_gap', max(abs(lam_inherit - lam_reselect)));
S.noise = struct('max_dIRF', mean(dN_irf), 'max_dIRF_over_postsd', mean(dN_rel), ...
                 'max_dlogtau', mean(dN_tau));

% ---- C. sandwich vs posterior-quantile bands -------------------------
S.C = struct('available', false);
if ~isempty(opts.mc_file) && exist(opts.mc_file, 'file') == 2
    L = load(opts.mc_file);
    s = summarize_montecarlo(L.mc);
    ib = find(strcmp(s.est_names, 'BLP-block'), 1);
    if ~isempty(ib) && any(~isnan(s.coverage_post(ib, :)))
        S.C.available = true;
        S.C.coverage_sandwich = mean(mean(s.coverage(ib, :, 2:end)));
        S.C.coverage_posterior = nanmean_(reshape(s.coverage_post(ib, :, 2:end), [], 1));
        S.C.width_sandwich = mean(mean(s.avg_len(ib, :, 2:end)));
        S.C.width_posterior = nanmean_(reshape(s.avg_len_post(ib, :, 2:end), [], 1));
        S.C.nominal = L.mc.meta.ci_level;
    end
else
    fprintf(['  [C] no stored Monte Carlo supplied (opts.mc_file): the ' ...
             'band comparison is a repeated-sampling statement and is ' ...
             'skipped.\n']);
end

print_and_write(S, opts);
end

% =====================================================================
function m = nanmean_(v)
v = v(~isnan(v));
if isempty(v), m = NaN; else, m = mean(v); end
end

function o = dflt(o, f, v)
if ~isfield(o, f) || isempty(o.(f)), o.(f) = v; end
end

function print_and_write(S, opts)
fprintf('\n--- NOISE BASELINE: the same estimator, a different seed ---\n');
fprintf('  (every number below must be read against these)\n');
fprintf('  max |IRF change| from re-running the chain alone             : %.4f\n', S.noise.max_dIRF);
fprintf('  the same, relative to the posterior sd of the IRF            : %.3f\n', S.noise.max_dIRF_over_postsd);
fprintf('  max |change in log tau|                                      : %.3f\n', S.noise.max_dlogtau);

fprintf('\n--- A. cross-equation covariance in the adaptive sampler ---\n');
fprintf('  max |IRF change| when sigma2 is fixed at the SYSTEM NIW value : %.4f  (noise %.4f)\n', ...
        S.A.max_dIRF, S.noise.max_dIRF);
fprintf('  the same, relative to the posterior sd of the IRF            : %.3f  (noise %.3f)\n', ...
        S.A.max_dIRF_over_postsd, S.noise.max_dIRF_over_postsd);
fprintf('  max |change in log tau|                                      : %.3f  (noise %.3f)\n', ...
        S.A.max_dlogtau, S.noise.max_dlogtau);
fprintf('  NOTE: these are MAXIMA over all (equation, horizon) cells of a\n');
fprintf('        replication, then averaged over replications -- a worst-cell\n');
fprintf('        measure, not a typical one.\n');
fprintf('\n--- B1. impact vector b1n held fixed ---\n');
fprintf('  sd of the IRF induced by b1n posterior uncertainty alone      : %.4f\n', S.B1.sd_from_b1n);
fprintf('  reported half band-width                                      : %.4f\n', S.B1.reported_half_band);
fprintf('  ratio (how much the bands understate, from this channel)      : %.3f\n', S.B1.ratio);
fprintf('\n--- B2. prior centre held at the BVAR posterior mean ---\n');
fprintf('  sd of the IRF across BVAR posterior draws of the centre       : %.4f\n', S.B2.sd_from_prior_centre);
fprintf('  ratio to the reported half band-width                         : %.3f\n', S.B2.ratio);
fprintf('\n--- B3. lambda inherited from the global BLP vs re-selected ---\n');
fprintf('  mean lambda inherited / re-selected                           : %.4f / %.4f\n', ...
        S.B3.lambda_inherited, S.B3.lambda_reselected);
fprintf('  largest gap between the two across datasets                   : %.2e\n', ...
        S.B3.max_lambda_gap);
fprintf('  max |IRF change| (same seed), and relative to the posterior sd : %.4f / %.3f\n', ...
        S.B3.max_dIRF, S.B3.max_dIRF_over_postsd);
if S.B3.max_lambda_gap < 1e-10
    fprintf(['  => the adaptive estimator''s own marginal-likelihood selection\n' ...
             '     returns EXACTLY the global BLP''s lambda (same objective, same\n' ...
             '     data, both evaluated at tau = 1), so inheriting it is not an\n' ...
             '     approximation at all.  What remains an approximation is that\n' ...
             '     lambda is selected at tau = 1 and never re-selected jointly\n' ...
             '     with tau; that would be a different, larger model.\n']);
end
if S.C.available
    fprintf('\n--- C. quasi-Bayesian sandwich vs posterior-quantile bands ---\n');
    fprintf('  nominal level                                                 : %.2f\n', S.C.nominal);
    fprintf('  coverage: sandwich %.3f, posterior quantiles %.3f\n', ...
            S.C.coverage_sandwich, S.C.coverage_posterior);
    fprintf('  average width: sandwich %.3f, posterior quantiles %.3f\n', ...
            S.C.width_sandwich, S.C.width_posterior);
end

d = fileparts(opts.out_csv);
if ~isempty(d) && exist(d, 'dir') ~= 7, mkdir(d); end
fid = fopen(opts.out_csv, 'w');
assert(fid > 0, 'run_sensitivity_approximations: cannot write %s', opts.out_csv);
fprintf(fid, 'approximation,quantity,value\n');
w = @(a, q, v) fprintf(fid, '%s,%s,%.10g\n', a, q, v);
w('A_cross_equation_covariance', 'max_dIRF', S.A.max_dIRF);
w('A_cross_equation_covariance', 'max_dIRF_over_posterior_sd', S.A.max_dIRF_over_postsd);
w('A_cross_equation_covariance', 'max_dlogtau', S.A.max_dlogtau);
w('B1_fixed_impact_vector', 'sd_from_b1n', S.B1.sd_from_b1n);
w('B1_fixed_impact_vector', 'reported_half_band', S.B1.reported_half_band);
w('B1_fixed_impact_vector', 'ratio', S.B1.ratio);
w('B2_fixed_prior_centre', 'sd_from_prior_centre', S.B2.sd_from_prior_centre);
w('B2_fixed_prior_centre', 'ratio', S.B2.ratio);
w('B3_fixed_lambda', 'lambda_inherited', S.B3.lambda_inherited);
w('B3_fixed_lambda', 'lambda_reselected', S.B3.lambda_reselected);
w('B3_fixed_lambda', 'max_dIRF', S.B3.max_dIRF);
w('B3_fixed_lambda', 'max_dIRF_over_posterior_sd', S.B3.max_dIRF_over_postsd);
w('B3_fixed_lambda', 'max_lambda_gap', S.B3.max_lambda_gap);
w('noise_baseline', 'max_dIRF', S.noise.max_dIRF);
w('noise_baseline', 'max_dIRF_over_posterior_sd', S.noise.max_dIRF_over_postsd);
w('noise_baseline', 'max_dlogtau', S.noise.max_dlogtau);
if S.C.available
    w('C_quasi_bayesian_bands', 'coverage_sandwich', S.C.coverage_sandwich);
    w('C_quasi_bayesian_bands', 'coverage_posterior', S.C.coverage_posterior);
    w('C_quasi_bayesian_bands', 'width_sandwich', S.C.width_sandwich);
    w('C_quasi_bayesian_bands', 'width_posterior', S.C.width_posterior);
    w('C_quasi_bayesian_bands', 'nominal', S.C.nominal);
end
fclose(fid);
fprintf('\nrun_sensitivity_approximations: wrote %s\n', opts.out_csv);
end
