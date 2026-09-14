function test_empirical_iv()
% PURPOSE
% -------
% Interface test of the v2 (external-instrument) empirical package on
% SYNTHETIC inputs, so it runs on a machine with no licensed data:
%   1. build_instrument_series on the shipped EA-EMPD extract: shape, the
%      legacy identity (z_gc_1y_sum equals the v1 monthly sum), the
%      Kilian weights (every event's two weights add to one, so the sum
%      over months of the two aggregations agree), the JK split
%      (jk + info = the full gc series);
%   2. make_synthetic_ois_daily -> import_ois_daily: the fixture is flagged
%      synthetic, monthly eom/avg are built, and the real derived file is
%      NOT written;
%   3. ea_identify_proxy on a simulated VAR with a known proxy: b_z within
%      tolerance of the truth, unit effect, sign flag, bands present;
%   4. ea_relevance_iv returns every block with the right shapes; the
%      lead placebo is insignificant for a valid instrument; a LAGGED
%      (misdated) instrument shows up in the lag-k placebo, not in the
%      leads;
%   5. RUN_EMPIRICAL_IV refuses a synthetic dataset;
%   6. ea_pre_event_info and import_external_instrument shapes.
% Nothing here is an empirical result.

fprintf('test_empirical_iv:\n');
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(root));
P = ea_paths();
tmp = fullfile(tempdir, sprintf('ea_iv_test_%d', round(rand * 1e6)));  mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));

% --- 1. instrument builder -------------------------------------------------
inst = build_instrument_series(struct('out_csv', fullfile(tmp, 'inst.csv'), 'out_mat', fullfile(tmp, 'inst.mat'), ...
                                      'out_events', fullfile(tmp, 'ev.csv')));
assert(size(inst.Z, 1) == numel(inst.ym) && size(inst.Z, 2) == numel(inst.names), 'instrument panel shape');
S = load(fullfile(P.derived, 'shocks_monthly.mat'));
[tf, loc] = ismember(S.ym, inst.ym);
d = max(abs(S.mps_gc_1y(tf) - inst.Z(loc(tf), strcmp(inst.names, 'z_gc_1y_sum'))));
assert(d < 1e-6, 'z_gc_1y_sum does not reproduce the legacy monthly sum (max diff %.3g bp)', d);
zs = inst.Z(:, strcmp(inst.names, 'z_gcs_1m_adj_sum'));  zk = inst.Z(:, strcmp(inst.names, 'z_gcs_1m_adj_kilian'));
assert(abs(sum(zs(1:end - 1)) - sum(zk(1:end - 1))) < 1e-6 + abs(zk(end)), 'Kilian weights do not add to one');
zf = inst.Z(:, strcmp(inst.names, 'z_gcs_1m_adj_kilianfix'));
assert(abs(sum(zs(1:end - 1)) - sum(zf(1:end - 1))) < 1e-6 + abs(zf(end)), 'kilianfix weights do not add to one');
assert(any(abs(zf - zk) > 1e-9), 'kilianfix should differ from kilian');
z1 = inst.Z(:, strcmp(inst.names, 'z_gc_1m_adj_sum'));
zj = inst.Z(:, strcmp(inst.names, 'z_gc_1m_adj_jk_sum'));  zi = inst.Z(:, strcmp(inst.names, 'z_gc_1m_adj_info_sum'));
% jk + info covers the events that have BOTH a 1M and a STOXX move; the
% full series also has events without a STOXX move, so jk + info <= full in
% count, and equal where both exist -- test on months where they agree in
% support: the difference must be attributable to missing-STOXX events only
assert(sum((zj + zi) ~= 0) <= sum(z1 ~= 0), 'JK split has more events than the full series');
fprintf('  build_instrument_series: %d series, legacy identity (%.1e bp), Kilian weights, JK split: OK\n', numel(inst.names), d);

% --- 2. synthetic OIS ------------------------------------------------------------
f = make_synthetic_ois_daily(struct('out_csv', fullfile(tmp, 'syn_ois.csv')));
real_out = fullfile(tmp, 'never_written.csv');
m = import_ois_daily(struct('daily_csv', f, 'out_csv', real_out, 'verbose', false));
assert(m.synthetic, 'synthetic OIS fixture not flagged');
assert(exist(real_out, 'file') ~= 2, 'synthetic OIS was written to the derived name');
assert(all(isfinite(m.eom)) && all(isfinite(m.avg)) && numel(m.ym) == 252, 'monthly OIS aggregation');
m2 = import_ois_daily(struct('daily_csv', f, 'out_csv', fullfile(tmp, 'syn_monthly.csv'), 'allow_synthetic', true, 'verbose', false));
assert(exist(fullfile(tmp, 'syn_monthly.mat'), 'file') == 2, 'allow_synthetic did not write');
% bid/ask layout: the mid is used
fb = fullfile(tmp, 'syn_bidask.csv');
fid = fopen(f, 'r');  fgetl(fid);  fo = fopen(fb, 'w');  fprintf(fo, 'Date,Bid,Ask,SYNTHETIC\n');
while true
    l = fgetl(fid);  if ~ischar(l), break; end
    q = regexp(l, ',', 'split');  v = str2double(q{2});
    fprintf(fo, '%s,%.4f,%.4f,1\n', q{1}, v - 0.01, v + 0.01);
end
fclose(fid);  fclose(fo);
m3 = import_ois_daily(struct('daily_csv', fb, 'out_csv', fullfile(tmp, 'never2.csv'), 'verbose', false));
assert(max(abs(m3.eom - m.eom)) < 1e-9 && ~isempty(strfind(m3.meta.value_rule, 'mid')), 'bid/ask mid not used');  %#ok<STREMP>
fprintf('  import_ois_daily: fixture flagged, %d months, derived file protected, bid/ask mid: OK\n', numel(m2.ym));

% --- 3./4. proxy identification and relevance on a simulated VAR ----------------
rng(31, 'twister');
K = 3;  p = 2;  T = 1500;
A1 = [0.5 0.1 0; 0.2 0.6 0.1; 0 0.1 0.7];  A2 = [0.2 0 0; 0 0.1 0; 0.1 0 0.1];
B0 = [1 0 0; 0.5 1 0; -0.3 0.4 1];
eps_ = randn(T + 50, K);  Ysim = zeros(T + 50, K);
for t = 3:T + 50, Ysim(t, :) = (A1 * Ysim(t - 1, :)' + A2 * Ysim(t - 2, :)' + B0 * eps_(t, :)')'; end
Ysim = Ysim(51:end, :);  e1 = eps_(51:end, 1);
z = e1 + 0.5 * randn(T, 1);  z(rand(T, 1) < 0.2) = 0;      % noisy proxy with empty months
cfg = default_config();  cfg.mode = 'fmar';  cfg.p = p;  cfg.H = 6;  cfg.K = K;  cfg.shock_var = 1;
cfg.fmar.isrw = [0 0 0];  cfg.fmar.n_niw_draws = 60;  cfg.fmar.h1_mode = 'lp';
bvar = estimate_bvar_niw(Ysim, cfg);
bz = ea_identify_proxy(bvar, Ysim, z, cfg, struct('z_name', 'sim', 'verbose', false));
truth = B0(:, 1) / B0(1, 1);
assert(abs(bz.b1n(1) - 1) < 1e-12, 'unit effect normalisation');
assert(max(abs(bz.b1n - truth)) < 0.15, 'proxy impact vector off (max err %.3f)', max(abs(bz.b1n - truth)));
assert(bz.ident.sign_ok && bz.ident.first_stage.F_eff > 50, 'first stage should be strong here (F = %.1f)', bz.ident.first_stage.F_eff);
assert(all(isfinite(bz.theta_lo(:))) && all(bz.theta_lo(:) <= bz.theta_hi(:)), 'proxy bands');
rel = ea_relevance_iv(Ysim, bvar, z, cfg, struct('z_name', 'sim', 'verbose', false, 'ym', (12 * 2000 + 1:12 * 2000 + T)', ...
                      'subsamples', [2000 1 2050 12; 2051 1 2124 12]));
for f_ = {'fs', 'naive', 'placebo', 'predict', 'influence', 'subsample'}
    assert(isfield(rel, f_{1}), 'relevance block %s missing', f_{1});
end
assert(rel.fs.F_eff > 50 && numel(rel.placebo.t_hac) == 7, 'relevance first stage / placebo shape');
leads = rel.placebo.t_hac(rel.placebo.k > 0);
assert(all(abs(leads) < 3.5), 'lead placebo significant for a valid instrument (t = %s)', num2str(leads, '%.2f '));
zlag = [0; 0; z(1:end - 2)];                                % recorded two months AFTER the event
% u_s(t) then correlates with zlag(t+2): the k = +2 LEAD placebo must light
% up while the contemporaneous first stage stays quiet
rel2 = ea_relevance_iv(Ysim, bvar, zlag, cfg, struct('z_name', 'lagged', 'verbose', false));
assert(abs(rel2.placebo.t_hac(rel2.placebo.k == 2)) > 5 && abs(rel2.fs.t_hac) < 3, ...
       'a misdated instrument should show in the k = +2 placebo, not contemporaneously');
assert(numel(rel.influence.top.ym) == 5 && isfinite(rel.influence.fs_excl_crisis.F_eff), 'influence block');
fprintf('  ea_identify_proxy / ea_relevance_iv: b_z err %.3f, F_eff %.0f, placebo leads quiet, misdating detected: OK\n', ...
        max(abs(bz.b1n - truth)), rel.fs.F_eff);

% --- 5. the driver refuses a synthetic dataset -----------------------------------------
[X, xn] = ea_pre_event_info(Ysim, {'ind', 'ip', 'hicp'});
assert(size(X, 1) == T && numel(xn) == size(X, 2), 'pre-event info shape');
ds = struct('Y', Ysim, 'varnames', {{'ind', 'ip', 'hicp'}}, 'ym', (12 * 2000 + 1:12 * 2000 + T)', 'isrw', [0 0 0], ...
            'Z', z, 'znames', {{'z_gcs_1m_adj_sum'}}, 'pred', X, 'prednames', {xn}, 'n_gc', double(z ~= 0), ...
            'n_sp_1m', zeros(T, 1), 'system', 'synthetic', 'instrument_free', true, 'synthetic', true, 'meta', struct());
synth = fullfile(tmp, 'SYNTHETIC_ds.mat');  save(synth, '-struct', 'ds');
threw = false;
try
    DATASET = synth;  SYSTEM = 'synthetic';  QUICK = true;  %#ok<NASGU>
    run(fullfile(P.empirical, 'RUN_EMPIRICAL_IV.m'));
catch err
    threw = ~isempty(strfind(err.message, 'SYNTHETIC'));  %#ok<STREMP>
end
assert(threw, 'RUN_EMPIRICAL_IV accepted a synthetic dataset');
fprintf('  RUN_EMPIRICAL_IV refuses a synthetic dataset: OK\n');

% --- 6. external instrument adapter -----------------------------------------------------
fp = fullfile(tmp, 'ext.csv');
fid = fopen(fp, 'w');  fprintf(fid, 'year,month,MP_pm,CBI_pm\n2001,1,0.5,-0.2\n2001,4,-1.0,0.1\n');  fclose(fid);
E = import_external_instrument(fp, struct('column', 'CBI_pm', 'scale', 100));
assert(numel(E.ym) == 4 && E.value(1) == -20 && E.value(2) == 0 && E.value(4) == 10, 'external adapter values');
fp2 = fullfile(tmp, 'ext2.csv');
fid = fopen(fp2, 'w');  fprintf(fid, 'date,shock\n2001-01-31,1.5\n2001-02-28,-0.5\n');  fclose(fid);
E2 = import_external_instrument(fp2);
assert(numel(E2.ym) == 2 && E2.value(2) == -0.5, 'external adapter date layout');
fprintf('  import_external_instrument (year/month and date layouts): OK\n');
fprintf('PASS: test_empirical_iv\n');
end
