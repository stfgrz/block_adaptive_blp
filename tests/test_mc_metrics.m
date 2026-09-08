function test_mc_metrics()
% PURPOSE
% -------
% Check the Monte Carlo evaluation layer: the metric definitions, the
% two integrated-RMSE conventions, the standardised reproducibility
% record, and the CSV export.  These are the objects every number in the
% write-up is read off, so they get their own test rather than being
% trusted because they "look right" in a printout.
%
% WHAT IS CHECKED
% ---------------
% 1. METRIC DEFINITIONS recomputed by hand from mc.theta / mc.lo / mc.hi
%    for a randomly chosen (estimator, variable, horizon): bias,
%    variance, MSE, RMSE, coverage and average interval length must agree
%    to machine precision with summarize_montecarlo's output.
% 2. THE BIAS-VARIANCE IDENTITY  MSE = bias^2 + variance.
% 3. THE TWO INTEGRATED RMSE CONVENTIONS: irmse averages h = 1..H,
%    irmse_h2 averages h = 2..H, and they must differ exactly by the
%    h = 1 term.  This is the metric the "fair benchmark" fix turns on,
%    so it is pinned.
% 4. h = 1 HARMONISATION: with cfg.fmar.h1_mode = 'lp' the global
%    baseline must report an LP at h = 1 (so its h = 1 IRF differs from
%    the BVAR's), while irmse_h2 is UNCHANGED by the switch -- that
%    invariance is what makes irmse_h2 the metric to quote.
% 5. REPRODUCIBILITY: re-running run_montecarlo with the same cfg gives
%    bit-identical estimates, and mc.meta carries every field a reader
%    needs to re-run the experiment.
% 6. CSV EXPORT: the files are written, the metrics file has the right
%    number of data rows, and a spot-checked value round-trips.
%
% The Monte Carlo here is deliberately tiny (R = 4, H = 4): this test is
% about the accounting, not about the statistics.

fprintf('test_mc_metrics:\n');

cfg = default_config();
cfg.mode = 'fmar';
cfg.K = 3;  cfg.H = 4;  cfg.T = 140;
cfg.mc.n_rep = 4;
cfg.mc.gibbs_n_burn = 60;  cfg.mc.gibbs_n_keep = 120;
cfg.mc.estimators = {'lp', 'var', 'global', 'block', 'pooled'};
cfg.fmar.n_niw_draws = 60;

mc = run_montecarlo(cfg, 'sparse');
s  = summarize_montecarlo(mc);

[nE, K, Hp1, R] = size(mc.theta);
H = Hp1 - 1;
assert(nE == 5, 'expected 5 estimators, got %d', nE);

% --- 1. metric definitions ----------------------------------------------
rng(1, 'twister');
for trial = 1:6
    e = randi(nE);  i = randi(K);  h = randi(H);
    est = squeeze(mc.theta(e, i, h + 1, :));
    tru = mc.theta_true(i, h + 1);
    lo  = squeeze(mc.lo(e, i, h + 1, :));
    hi  = squeeze(mc.hi(e, i, h + 1, :));
    chk(s.bias(e, i, h),     mean(est - tru),                 'bias');
    chk(s.variance(e, i, h), mean((est - mean(est)).^2),      'variance');
    chk(s.mse(e, i, h),      mean((est - tru).^2),            'mse');
    chk(s.rmse(e, i, h),     sqrt(mean((est - tru).^2)),      'rmse');
    chk(s.coverage(e, i, h), mean(lo <= tru & tru <= hi),     'coverage');
    chk(s.avg_len(e, i, h),  mean(hi - lo),                   'avg_len');
end
fprintf('  metric definitions recomputed by hand: OK\n');

% --- 2. bias-variance identity ------------------------------------------
d = max(abs(s.mse(:) - (s.bias(:).^2 + s.variance(:))));
assert(d < 1e-12, 'MSE ~= bias^2 + variance (max gap %.3g)', d);
fprintf('  MSE = bias^2 + variance to %.1e: OK\n', d);

% --- 3. the two integrated conventions ----------------------------------
chk(max(max(abs(s.irmse - mean(s.rmse, 3)))), 0, 'irmse definition');
chk(max(max(abs(s.irmse_h2 - mean(s.rmse(:, :, 2:H), 3)))), 0, 'irmse_h2 definition');
implied = (H * s.irmse - s.rmse(:, :, 1)) / (H - 1);
chk(max(max(abs(implied - s.irmse_h2))), 0, 'irmse/irmse_h2 consistency');
fprintf('  irmse (h=1..%d) and irmse_h2 (h=2..%d) consistent: OK\n', H, H);

% --- 4. h = 1 harmonisation ----------------------------------------------
cfg_lp = cfg;  cfg_lp.fmar.h1_mode = 'lp';
mc_lp = run_montecarlo(cfg_lp, 'sparse');
s_lp  = summarize_montecarlo(mc_lp);
i_glob = find(strcmp(mc.est_names, 'BLP-FMAR'), 1);
i_bvar = find(strcmp(mc.est_names, 'BVAR'), 1);
d1_bvar = max(max(abs(squeeze(mc.theta(i_glob, :, 2, :)) - ...
                     squeeze(mc.theta(i_bvar, :, 2, :)))));
assert(d1_bvar < 1e-12, ...
    'with h1_mode=bvar the global BLP should EQUAL the BVAR at h=1 (gap %.3g)', d1_bvar);
d1_lp = max(max(abs(squeeze(mc_lp.theta(i_glob, :, 2, :)) - ...
                    squeeze(mc_lp.theta(i_bvar, :, 2, :)))));
assert(d1_lp > 1e-6, ...
    'with h1_mode=lp the global BLP should DIFFER from the BVAR at h=1');
% irmse_h2 must be EXACTLY invariant for the estimators whose h >= 2
% behaviour cannot depend on the h = 1 convention: LP, the BVAR and the
% global BLP are all computed before, or independently of, the h = 1
% choice.
det_names = {'LP', 'BVAR', 'BLP-FMAR'};
for k = 1:numel(det_names)
    e = find(strcmp(mc.est_names, det_names{k}), 1);
    dd = max(abs(s.irmse_h2(e, :) - s_lp.irmse_h2(e, :)));
    assert(dd < 1e-12, ...
        'irmse_h2 of %s must not depend on h1_mode (gap %.3g)', det_names{k}, dd);
end
% The two SAMPLED estimators shift a little, for two reasons that are
% both about the h = 1 regression and not about the h >= 2 estimand:
%   (a) they run an h = 1 Gibbs chain whose lambda changes with the
%       convention (under 'bvar' they inherit the BVAR's Minnesota
%       tightness, which is not an LP tightness at all), so the random
%       stream reaching h >= 2 differs;
%   (b) the POOLED estimator links tau across horizons on purpose, so a
%       different h = 1 fit genuinely, if weakly, moves later horizons.
% The shift must be small relative to the metric itself.
for nm = {'BLP-block', 'BLP-pooled'}
    e = find(strcmp(mc.est_names, nm{1}), 1);
    rel = max(abs(s.irmse_h2(e, :) - s_lp.irmse_h2(e, :)) ./ s.irmse_h2(e, :));
    assert(rel < 0.10, ...
        '%s irmse_h2 moved %.1f%% with h1_mode -- too much to be stream noise', ...
        nm{1}, 100 * rel);
end
d_legacy = max(max(abs(s.irmse - s_lp.irmse)));
assert(d_legacy > 1e-8, ...
    'the legacy irmse should move with h1_mode -- that is why it is not preferred');
fprintf(['  h1_mode: BVAR-equality at h=1 under ''bvar'', LP under ''lp''; ' ...
         'irmse_h2 exactly invariant for LP/BVAR/BLP-FMAR, legacy irmse ' ...
         'moves by %.3f: OK\n'], d_legacy);

% --- 5. reproducibility ---------------------------------------------------
mc2 = run_montecarlo(cfg, 'sparse');
assert(isequal(mc.theta, mc2.theta), 'run_montecarlo is not reproducible');
assert(isequal(mc.seeds, mc2.seeds), 'seed sequence changed between runs');
need = {'dgp_name', 'dgp_params', 'misspec_block', 'mode', 'R', 'T', 'p', ...
        'K', 'H', 'shock_var', 'ci_level', 'master_seed', 'seed_rule', ...
        'est_keys', 'est_names', 'blocks_scheme', 'gibbs_n_burn', ...
        'gibbs_n_keep', 'lp_vcov', 'created'};
for k = 1:numel(need)
    assert(isfield(mc.meta, need{k}), 'mc.meta lacks .%s', need{k});
end
assert(isfield(s, 'meta'), 'the summary does not carry .meta');
fprintf('  bit-identical re-run and complete mc.meta (%d fields): OK\n', ...
        numel(fieldnames(mc.meta)));

% --- 6. csv export --------------------------------------------------------
stem = fullfile(tempdir, sprintf('mc_metrics_test_%d', round(rand * 1e6)));
files = export_montecarlo_csv(s, stem);
assert(numel(files) >= 4, 'expected at least 4 exported files, got %d', numel(files));
for k = 1:numel(files)
    assert(exist(files{k}, 'file') == 2, 'missing export %s', files{k});
end
n_lines = count_lines([stem '_metrics.csv']);
assert(n_lines == nE * K * H + 1, ...
    'metrics csv has %d lines, expected %d', n_lines, nE * K * H + 1);
% spot-check one exported value against the summary
% The export prints %.10g, so the round-trip is exact to 10 significant
% digits -- far more than any table needs, but not to 1e-12.
v = read_csv_value([stem '_metrics.csv'], s.est_names{2}, 1, 2, 'rmse');
chk(v, s.rmse(2, 1, 2), 'csv rmse round-trip', 1e-9);
for k = 1:numel(files), delete(files{k}); end
fprintf('  csv export: %d files, %d metric rows, values round-trip: OK\n', ...
        numel(files), n_lines - 1);

fprintf('PASS: test_mc_metrics\n\n');
end

% =====================================================================
function chk(a, b, what, tol)
if nargin < 4 || isempty(tol), tol = 1e-12; end
assert(abs(a - b) < tol * max(1, abs(b)), ...
    'test_mc_metrics: %s mismatch (%.12g vs %.12g)', what, a, b);
end

function n = count_lines(f)
fid = fopen(f, 'r');  n = 0;
while ischar(fgetl(fid)), n = n + 1; end
fclose(fid);
end

function v = read_csv_value(f, est, ivar, ih, col)
% Pull one value out of the tidy metrics csv by (estimator, variable,
% horizon) and column name -- i.e. read it back the way a reader would.
fid = fopen(f, 'r');
hdr = strsplit(fgetl(fid), ',');
ic = find(strcmp(hdr, col), 1);
ie = find(strcmp(hdr, 'estimator'), 1);
iv = find(strcmp(hdr, 'variable'), 1);
ih_col = find(strcmp(hdr, 'horizon'), 1);
v = NaN;
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    f2 = strsplit(l, ',');
    if strcmp(f2{ie}, est) && str2double(f2{iv}) == ivar && ...
            str2double(f2{ih_col}) == ih
        v = str2double(f2{ic});  break
    end
end
fclose(fid);
assert(~isnan(v), 'read_csv_value: row not found');
end
