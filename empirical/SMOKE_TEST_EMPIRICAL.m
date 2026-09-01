% SMOKE_TEST_EMPIRICAL  Five-minute check before committing a night to
% RUN_EMPIRICAL.
%
% PURPOSE
% -------
% RUN_EMPIRICAL is a multi-hour job whose first minute exercises almost
% none of what it depends on.  This script runs the SAME four estimators
% on the SAME real dataset, at a tiny p, H and chain length, and asserts
% every interface RUN_EMPIRICAL and run_null_calibration rely on:
%
%   * the dataset exists, has the documented shape, and carries isrw;
%   * estimate_bvar_niw accepts the 1 x K isrw vector and returns the
%     fields the drivers read (.B/.c/.A layout, .b1n, .theta, .lambda,
%     .max_eig, .F, .Psi, .theta_lo/.theta_hi);
%   * estimate_blp_fmar(Y, cfg, bvar) returns .lambda (K x H) and
%     .theta_mean;
%   * estimate_blp_blockadaptive(Y, cfg, bvar, lambda) returns
%     .tau_mean (K x G x H, G = K) and .lo/.hi (K x (H+1));
%   * estimate_lp_lagaug(Y, cfg, bvar) returns .theta (K x (H+1));
%   * simulate_fitted_bvar_dgp runs (its internal B vs c/A layout
%     assertion is what makes a signature mismatch fail loudly), which is
%     the one piece run_null_calibration adds on top.
%
% If this prints ALL CHECKS PASSED, RUN_EMPIRICAL will not fall over on
% an interface problem.  It says nothing about how long the real run
% takes, or about whether the results are interesting.
%
% USAGE
% -----
%   SMOKE_TEST_EMPIRICAL          % from anywhere; needs data/ea_dataset.mat
%
% NOTES
% -----
% * Deliberately no plotting: this must pass on a headless machine.
% * Runtime is a couple of minutes at p = 2, H = 6, 300 draws.

% --- path bootstrap (works from any directory, and when pasted) ----------
ea_this = mfilename('fullpath');
if isempty(ea_this), ea_this = fullfile(pwd, 'SMOKE_TEST_EMPIRICAL'); end
addpath(genpath(fileparts(fileparts(ea_this))));
P = ea_paths();

fprintf('\n=== SMOKE_TEST_EMPIRICAL ===\n');
t_all = tic;

% --- 0. dataset ----------------------------------------------------------
assert(exist(P.dataset, 'file') == 2, ...
       ['SMOKE_TEST_EMPIRICAL: %s not found.\n' ...
        'Build it first:\n' ...
        '    build_shock_series(); fetch_outcome_data(); assemble_dataset();'], ...
       P.dataset);
ds = load(P.dataset);
for f = {'Y', 'varnames', 'ym', 'isrw', 'shock_variant'}
    assert(isfield(ds, f{1}), 'dataset is missing field .%s', f{1});
end
[T, K] = size(ds.Y);
assert(K == 5, 'expected K = 5 variables, found %d', K);
assert(numel(ds.isrw) == K, 'isrw has %d entries, expected %d', numel(ds.isrw), K);
assert(numel(ds.varnames) == K, 'varnames has %d entries, expected %d', ...
       numel(ds.varnames), K);
assert(all(isfinite(ds.Y(:))), 'dataset contains non-finite values');
fprintf('[0] dataset OK: T = %d, K = %d, variant = %s, isrw = [%s]\n', ...
        T, K, ds.shock_variant, strtrim(sprintf('%d ', ds.isrw)));

% --- configuration: same shape as RUN_EMPIRICAL, tiny everywhere --------
cfg = default_config();
cfg.mode      = 'fmar';
cfg.shock_var = 1;
cfg.p         = 2;
cfg.H         = 6;
cfg.ci_level  = 0.90;
cfg.K         = K;
cfg.gibbs.n_burn = 100;
cfg.gibbs.n_keep = 200;
cfg.fmar.n_niw_draws = 100;
cfg.fmar.isrw = ds.isrw;          % the vectorised prior mean, the whole point
rng(20260107, 'twister');

% --- 1. BVAR --------------------------------------------------------------
t = tic;
bv = estimate_bvar_niw(ds.Y, cfg);
fprintf('[1] estimate_bvar_niw OK (%.1f s): lambda = %.3f, max |eig| = %.4f\n', ...
        toc(t), bv.lambda, bv.max_eig);
need = {'B', 'c', 'A', 'Sigma', 'F', 'Psi', 'b1n', 'theta', ...
        'theta_lo', 'theta_hi', 'lambda', 'max_eig'};
for f = need
    assert(isfield(bv, f{1}), 'bvar is missing field .%s', f{1});
end
assert(isequal(size(bv.B), [1 + K * cfg.p, K]), 'bvar.B has the wrong size');
assert(isequal(size(bv.A), [K, K, cfg.p]), 'bvar.A has the wrong size');
assert(isequal(size(bv.theta), [K, cfg.H + 1]), 'bvar.theta has the wrong size');
assert(numel(bv.b1n) == K, 'bvar.b1n has the wrong size');
assert(all(isfinite(bv.theta(:))), 'bvar.theta is not finite');
% b1n is the recursive impact vector normalised on the shock variable
assert(abs(bv.b1n(cfg.shock_var) - 1) < 1e-10, ...
       'bvar.b1n is not normalised to 1 on the shock variable');
% The 25 bp normalisation RUN_EMPIRICAL applies divides by this:
fprintf('    impact of the surprise on i1y = %.4f pp  =>  25bp scale k = %.2f\n', ...
        bv.theta(2, 1), 0.25 / bv.theta(2, 1));
if ~(bv.theta(2, 1) > 0)
    fprintf(2, '    WARNING: non-positive impact response; the 25bp normalisation\n');
    fprintf(2, '             in RUN_EMPIRICAL would flip the sign of every IRF.\n');
end

% --- 2. BLP-FMAR ----------------------------------------------------------
t = tic;
bf = estimate_blp_fmar(ds.Y, cfg, bv);
fprintf('[2] estimate_blp_fmar OK (%.1f s)\n', toc(t));
assert(isfield(bf, 'lambda') && isequal(size(bf.lambda), [K, cfg.H]), ...
       'blp_fmar.lambda must be K x H');
assert(isequal(size(bf.theta_mean), [K, cfg.H + 1]), 'blp_fmar.theta_mean size');
assert(all(isfinite(bf.theta_mean(:))), 'blp_fmar.theta_mean is not finite');

% --- 3. BLP-block ---------------------------------------------------------
t = tic;
bb = estimate_blp_blockadaptive(ds.Y, cfg, bv, bf.lambda);
fprintf('[3] estimate_blp_blockadaptive OK (%.1f s)\n', toc(t));
assert(isequal(size(bb.tau_mean), [K, K, cfg.H]), ...
       'blp_block.tau_mean must be K x G x H = %d x %d x %d, got %s', ...
       K, K, cfg.H, mat2str(size(bb.tau_mean)));
assert(isequal(size(bb.lo), [K, cfg.H + 1]) && isequal(size(bb.hi), [K, cfg.H + 1]), ...
       'blp_block .lo/.hi must be K x (H+1)');
assert(all(isfinite(bb.tau_mean(:))) && all(bb.tau_mean(:) > 0), ...
       'blp_block.tau_mean must be finite and positive');
assert(all(bb.lo(:, 2:end) <= bb.hi(:, 2:end)), 'blp_block bands are inverted');
fprintf('    tau_mean is %d x %d x %d, range [%.3f, %.3f]\n', ...
        K, K, cfg.H, min(bb.tau_mean(:)), max(bb.tau_mean(:)));

% --- 4. LP-LA -------------------------------------------------------------
t = tic;
lp = estimate_lp_lagaug(ds.Y, cfg, bv);
fprintf('[4] estimate_lp_lagaug OK (%.1f s)\n', toc(t));
assert(isequal(size(lp.theta), [K, cfg.H + 1]), 'lp_lagaug.theta size');
assert(all(isfinite(lp.theta(:))), 'lp_lagaug.theta is not finite');

% --- 5. the null-calibration extra: simulate from the fitted BVAR --------
t = tic;
dgp = simulate_fitted_bvar_dgp(bv, ds.Y, cfg, struct('method', 'resample'));
fprintf('[5] simulate_fitted_bvar_dgp OK (%.1f s)\n', toc(t));
assert(isequal(size(dgp.Y), size(ds.Y)), 'simulated dataset has the wrong size');
assert(isequal(size(dgp.theta), [K, cfg.H + 1]), 'dgp.theta size');
assert(all(isfinite(dgp.Y(:))), 'simulated data is not finite');

% --- 6. the statistic RUN_EMPIRICAL step C / the protocol actually use ---
h1 = 2;  h2 = min(12, cfg.H);
tau_bar = mean(bb.tau_mean(:, :, h1:h2), 3);
assert(isequal(size(tau_bar), [K, K]), 'early-h tau_bar must be K x K');
fprintf('[6] early-h tau_bar (h = %d..%d), rows = equations, cols = blocks:\n', h1, h2);
fprintf('        %s\n', sprintf('%8s', ds.varnames{:}));
for i = 1:K
    fprintf('    %-6s', ds.varnames{i});
    fprintf('%8.3f', tau_bar(i, :));
    fprintf('\n');
end

fprintf('\nALL CHECKS PASSED in %.1f s.\n', toc(t_all));
fprintf('RUN_EMPIRICAL should run clean; launch it with:\n');
fprintf('    RUN_EMPIRICAL\n\n');
