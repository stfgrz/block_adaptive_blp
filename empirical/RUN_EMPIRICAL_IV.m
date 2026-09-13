% RUN_EMPIRICAL_IV  Chapter 7 v2 driver: external-instrument identification
% on an instrument-free state vector (docs/CH7_REDESIGN.md).
%
% The legacy baseline (internal instrument) is RUN_EMPIRICAL.m and is not
% touched.  Everything here writes results/iv_<system>_p<p>_* files.
%
% PREREQUISITES
%   build_instrument_series();                 (public EA-EMPD extract, shipped)
%   assemble_dataset_v2(struct('system', S));  S in {lev4, lev4_yoy, ois4, ois6}
%   (ois4 / ois6 need the daily EUREON1M= file and ea_fetch_v2_series(); see
%   docs/MONDAY_DATA_CHECKLIST.md)
%
% STEPS (toggles; preset any of them in the workspace before running)
%   REL  relevance / timing / influence table for EVERY instrument in ds.Z,
%        full report for the headline instrument
%   A    K-variable system at p = P_BASE: BVAR, BLP-FMAR, BLP-block,
%        BLP-pooled, LP-LA, all identified with the HEADLINE instrument
%        (ea_identify_proxy), plus LP-IV with Anderson-Rubin sets; point
%        IRFs for the ALTERNATIVE instruments from the same posterior
%        coefficient blocks (the tau map does not depend on identification)
%   B    tau heatmap export and the reading protocol against the design's
%        null (null_calibration_iv_<system>_p<p>.mat; for lev4 at p = 12 the
%        legacy null of the same data, null_calibration_p12_levels.mat, is
%        accepted with a printed note)
%   C    LP-IV table for headline + alternatives (theta, HAC band, AR set)
%   D    cross-p comparison: tau maps at p in PLIST, read through each
%        design's own null (ea_cross_p_table); ratio / percentile / counts
%   E    out-of-sample block ablation on the protocol's escaping cells (or
%        the top cell) and a quiet control (run_block_ablation)
%
% HEADLINE INSTRUMENT RULE (pre-specified; HEADLINE_Z overrides)
%   surprise maturity = indicator maturity: indicator ois1m* / eonia* /
%   euribor1m -> '1m_adj'; i1y -> '1y'.  Event set 'gcs' (meetings +
%   speeches; AGKL).  Aggregation by the indicator's timing convention:
%   monthly average of same-day closes (ois1m_avg, eonia_avg) -> 'kilian';
%   monthly average of a pre-event fixing (i1y, euribor1m) -> 'kilianfix';
%   end-of-month (ois1m, ois1m_eom) -> 'sum'.  Alternatives: the 'gc'
%   version, the other aggregation, the JK version, the '_bs' version, the
%   path proxy, and every external series (z_ext_*).
%
% QUICK = true runs everything with tiny settings (interface check on the
% real data; results meaningless).  DATASET overrides the dataset file.

if ~exist('SYSTEM', 'var'),     SYSTEM = 'lev4'; end
if ~exist('P_BASE', 'var'),     P_BASE = 12; end
if ~exist('H_IRF', 'var'),      H_IRF = 48; end
if ~exist('PLIST', 'var'),      PLIST = [2 4 6 12]; end
if ~exist('DO_REL', 'var'),     DO_REL = true; end
if ~exist('DO_A', 'var'),       DO_A = true; end
if ~exist('DO_B', 'var'),       DO_B = true; end
if ~exist('DO_C', 'var'),       DO_C = true; end
if ~exist('DO_D', 'var'),       DO_D = true; end
if ~exist('DO_E', 'var'),       DO_E = true; end
if ~exist('DO_POOLED', 'var'),  DO_POOLED = true; end
if ~exist('HEADLINE_Z', 'var'), HEADLINE_Z = ''; end
if ~exist('ALT_Z', 'var'),      ALT_Z = {}; end
if ~exist('NULL_FILE', 'var'),  NULL_FILE = ''; end
if ~exist('N_ORIGINS', 'var'),  N_ORIGINS = 60; end
if ~exist('ORIGIN0_YM', 'var'), ORIGIN0_YM = [2010 1]; end
if ~exist('QUICK', 'var'),      QUICK = false; end
if ~exist('DATASET', 'var'),    DATASET = ''; end

ea_this = mfilename('fullpath');
if isempty(ea_this), ea_this = fullfile(pwd, 'RUN_EMPIRICAL_IV'); end
addpath(genpath(fileparts(fileparts(ea_this))));
P = ea_paths();
if exist(P.results, 'dir') ~= 7, mkdir(P.results); end

if isempty(DATASET), DATASET = fullfile(P.data, sprintf('ea_dataset_v2_%s.mat', SYSTEM)); end
assert(exist(DATASET, 'file') == 2, ...
       'RUN_EMPIRICAL_IV: %s not found. Run assemble_dataset_v2(struct(''system'', ''%s'')).', DATASET, SYSTEM);
ds = load(DATASET);
if isfield(ds, 'synthetic') && ds.synthetic
    error(['%s: this dataset is a SYNTHETIC FIXTURE (ds.synthetic = true), built for ' ...
           'interface testing only; no result computed from it is an empirical finding.'], mfilename);
end
assert(isfield(ds, 'Z') && isfield(ds, 'znames'), 'RUN_EMPIRICAL_IV: the dataset carries no instruments; rebuild with assemble_dataset_v2.');
Y = ds.Y;  [T, K] = size(Y);  vn = ds.varnames;  ym = ds.ym;

cfg = default_config();
cfg.mode = 'fmar';  cfg.fmar.h1_mode = 'lp';  cfg.fmar.psi_floor = true;
cfg.blocks.fixed_tau = [];                 % no instrument in the state vector: every block adapts
cfg.shock_var = 1;  cfg.p = P_BASE;  cfg.H = H_IRF;  cfg.K = K;  cfg.ci_level = 0.90;
cfg.fmar.isrw = ds.isrw;
if QUICK
    cfg.H = 8;  cfg.gibbs.n_burn = 60;  cfg.gibbs.n_keep = 120;  cfg.fmar.n_niw_draws = 50;
    PLIST = PLIST(PLIST <= 4);  if isempty(PLIST), PLIST = 2; end
    N_ORIGINS = 4;
end
h_early = [2, min(12, cfg.H)];
tag = sprintf('iv_%s_p%d', SYSTEM, cfg.p);
if QUICK, tag = [tag '_quick']; end
[dkey, dinfo] = ea_design_key(ds, cfg, struct('h_early', h_early));

% --- headline and alternative instruments ------------------------------------
ind = vn{1};
% aggregation by the indicator's timing convention (CH7_REDESIGN Sec. 3):
% end-of-month level -> sum; monthly average of same-day closes (OIS, EONIA)
% -> kilian; monthly average of a fixing set BEFORE the events (Euribor,
% 11:00 CET) -> kilianfix
if any(strcmp(ind, {'i1y', 'euribor1m'})),              agg = 'kilianfix'; agg_other = 'sum';
elseif any(strcmp(ind, {'ois1m_avg', 'eonia_avg'})),   agg = 'kilian';    agg_other = 'sum';
else,                                                  agg = 'sum';       agg_other = 'kilian'; end
if strcmp(ind, 'i1y'), sur = '1y'; else, sur = '1m_adj'; end
if isempty(HEADLINE_Z), HEADLINE_Z = sprintf('z_gcs_%s_%s', sur, agg); end
if isempty(ALT_Z)
    ALT_Z = {sprintf('z_gc_%s_%s', sur, agg), sprintf('z_gcs_%s_%s', sur, agg_other), ...
             sprintf('z_gc_%s_jk_%s', sur, agg), [HEADLINE_Z '_bs'], sprintf('z_gc_path1y_%s', agg)};
    ALT_Z = [ALT_Z, ds.znames(strncmp(ds.znames, 'z_ext_', 6))];
end
ALT_Z = ALT_Z(ismember(ALT_Z, ds.znames) & ~strcmp(ALT_Z, HEADLINE_Z));
zcol = @(nm) ds.Z(:, strcmp(ds.znames, nm));
assert(any(strcmp(ds.znames, HEADLINE_Z)), 'RUN_EMPIRICAL_IV: headline instrument %s not in the dataset', HEADLINE_Z);
z_head = zcol(HEADLINE_Z);
ev_of = @(nm) tern(strncmp(nm, 'z_gc_', 5), ds.n_gc > 0, ds.n_gc > 0 | ds.n_sp_1m > 0);

fprintf('RUN_EMPIRICAL_IV: %s\n  system %s: T = %d, K = %d [%s], p = %d, H = %d, window [%d %d], tag %s\n  design key %s\n  headline instrument %s (indicator %s: maturity %s, aggregation %s); %d alternatives\n', ...
        DATASET, SYSTEM, T, K, strjoin(vn, ', '), cfg.p, cfg.H, h_early, tag, dkey, HEADLINE_Z, ind, sur, agg, numel(ALT_Z));

% ==============================================================================
if DO_REL
fprintf('\n[REL] relevance / timing / influence for %d instruments ...\n', numel(ds.znames));
rng(20260901, 'twister');
bv0 = estimate_bvar_niw(Y, cfg);
out_csv = fullfile(P.results, sprintf('%s_relevance.csv', tag));
fid = fopen(out_csv, 'w');
fprintf(fid, ['instrument,b,t_ehw,t_hac,F_eff,R2,N,n_nonzero,sd_z_bp,naive_b,naive_t,naive_F,' ...
              'lead1_t,lead2_t,lead3_t,lag1_t,pred_F,pred_p,ar1,F_excl_crisis,F_winsor,F_2001_08,F_2009_11,F_2012_19,top_month,top_dfbeta,mop_verdict\n']);
REL = struct();
for j = 1:numel(ds.znames)
    nm = ds.znames{j};  z = zcol(nm);
    if all(~isfinite(z)) || nnz(z(isfinite(z))) < 20, continue; end
    r = ea_relevance_iv(Y, bv0, z, cfg, struct('ym', ym, 'z_name', nm, 'predictors', ds.pred, ...
                        'pred_names', {ds.prednames}, 'event_months', ev_of(nm), 'verbose', false));
    REL.(regexprep(nm, '[^A-Za-z0-9_]', '_')) = r;
    tk = @(k) r.placebo.t_hac(r.placebo.k == k);
    fsub = nan(1, 3);
    for q = 1:numel(r.subsample)
        w = r.subsample(q).window(1);
        if w == 2001, fsub(1) = r.subsample(q).fs.F_eff; elseif w == 2009, fsub(2) = r.subsample(q).fs.F_eff; elseif w == 2012, fsub(3) = r.subsample(q).fs.F_eff; end
    end
    pf = NaN;  pp = NaN;  if r.predict.available, pf = r.predict.F_hac;  pp = r.predict.p_value; end
    fw = NaN;  if isfield(r.influence, 'fs_winsorised'), fw = r.influence.fs_winsorised.F_eff; end
    fprintf(fid, '%s,%.6f,%.3f,%.3f,%.3f,%.4f,%d,%d,%.3f,%.6f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d-%02d,%.3f,%s\n', ...
            nm, r.fs.b, r.fs.t_ehw, r.fs.t_hac, r.fs.F_eff, r.fs.R2, r.fs.N, r.fs.n_nonzero, r.fs.sd_z, ...
            r.naive.b, r.naive.t_hac, r.naive.F_eff, tk(1), tk(2), tk(3), tk(-1), pf, pp, r.predict.ar1_all, ...
            r.influence.fs_excl_crisis.F_eff, fw, fsub, r.influence.top.year(1), r.influence.top.month(1), ...
            r.influence.top.dfbeta(1), r.fs.mop_verdict);
end
fclose(fid);
fprintf('  wrote %s\n', out_csv);
fprintf('  headline instrument, full report:\n');
relh = ea_relevance_iv(Y, bv0, z_head, cfg, struct('ym', ym, 'z_name', HEADLINE_Z, 'predictors', ds.pred, ...
                       'pred_names', {ds.prednames}, 'event_months', ev_of(HEADLINE_Z), 'verbose', true, ...
                       'out_csv', fullfile(P.results, sprintf('%s_relevance_headline.csv', tag))));
save(fullfile(P.results, sprintf('%s_relevance.mat', tag)), 'REL', 'relh', 'HEADLINE_Z', 'ALT_Z', 'dkey', 'dinfo', 'tag');
end

% ==============================================================================
step_a_file = fullfile(P.results, sprintf('%s.mat', tag));
if DO_A
fprintf('\n[A] %d-variable system at p = %d, H = %d, proxy-identified with %s ...\n', K, cfg.p, cfg.H, HEADLINE_Z);
t_a = tic;
rng(20260107, 'twister');
bvar = estimate_bvar_niw(Y, cfg);
fprintf('  BVAR: lambda %.3f, max |eig| %.4f (%.1f s)\n', bvar.lambda, bvar.max_eig, toc(t_a));
bz = ea_identify_proxy(bvar, Y, z_head, cfg, struct('z_name', HEADLINE_Z));
blpf = estimate_blp_fmar(Y, cfg, bz);
fprintf('  BLP-FMAR (%.1f s); lambda multimodal at %d of %d horizons\n', toc(t_a), blpf.diag.lambda_multimodal, cfg.H);
blpb = estimate_blp_blockadaptive(Y, cfg, bz, blpf.lambda);
fprintf('  BLP-block (%.1f s)\n', toc(t_a));
if DO_POOLED, blpp = estimate_blp_blockpooled(Y, cfg, bz, blpf.lambda);  fprintf('  BLP-pooled (%.1f s)\n', toc(t_a)); else, blpp = []; end
lpla = estimate_lp_lagaug(Y, cfg, bz);
lpiv = estimate_lp_iv(Y, z_head, cfg);
fprintf('  LP-IV first stage F_eff %.2f; AR set unbounded at %d of %d cells (h >= 1)\n', ...
        lpiv.first_stage.F_eff, sum(sum(lpiv.ar_unbounded(:, 2:end))), K * cfg.H);
% alternative instruments: BVAR (with bands) and point IRFs from the SAME posterior blocks
ALT = struct('name', {}, 'bz', {}, 'theta_fmar', {}, 'theta_block', {}, 'theta_pooled', {}, 'lpiv', {});
for j = 1:numel(ALT_Z)
    za = zcol(ALT_Z{j});
    if all(~isfinite(za)) || nnz(za(isfinite(za))) < 20, continue; end
    bza = ea_identify_proxy(bvar, Y, za, cfg, struct('z_name', ALT_Z{j}, 'verbose', false));
    e = struct('name', ALT_Z{j}, 'bz', slim_bvar(bza), ...
               'theta_fmar', reproject(blpf.beta_mean, bza.b1n, K, cfg.H, 'fmar'), ...
               'theta_block', reproject(blpb.beta_mean, bza.b1n, K, cfg.H, 'fmar'), ...
               'theta_pooled', [], 'lpiv', estimate_lp_iv(Y, za, cfg));
    if DO_POOLED, e.theta_pooled = reproject(blpp.beta_mean, bza.b1n, K, cfg.H, 'fmar'); end
    ALT(end + 1) = e;  %#ok<SAGROW>
    fprintf('  alt %-26s F_eff %6.2f  b_z:', ALT_Z{j}, bza.ident.first_stage.F_eff);  fprintf(' %7.2f', bza.b1n);  fprintf('\n');
end
ds_info = struct('varnames', {vn}, 'ym', ym, 'system', SYSTEM, 'T', T, 'meta', ds.meta);
bvar = slim_bvar(bvar);  bz = slim_bvar(bz);
save(step_a_file, 'bvar', 'bz', 'blpf', 'blpb', 'blpp', 'lpla', 'lpiv', 'ALT', 'cfg', 'h_early', 'ds_info', ...
     'tag', 'HEADLINE_Z', 'ALT_Z', 'dkey', 'dinfo', 'T');
% IRF table for the headline (25 bp scale = 0.25 x unit-effect IRF)
out_csv = fullfile(P.results, sprintf('%s_irf.csv', tag));
fid = fopen(out_csv, 'w');
fprintf(fid, 'instrument,estimator,variable,h,theta_25bp,lo_25bp,hi_25bp\n');
ests = {'BVAR', bz.theta, bz.theta_lo, bz.theta_hi; 'BLP-FMAR', blpf.theta_mean, blpf.lo, blpf.hi; ...
        'BLP-block', blpb.theta_mean, blpb.lo, blpb.hi; 'LP-LA', lpla.theta, lpla.lo, lpla.hi; ...
        'LP-IV', lpiv.theta, lpiv.lo, lpiv.hi; 'LP-IV-AR', lpiv.theta, lpiv.ar_lo, lpiv.ar_hi};
if DO_POOLED, ests(end + 1, :) = {'BLP-pooled', blpp.theta_mean, blpp.lo, blpp.hi}; end
for e = 1:size(ests, 1)
    for i = 1:K
        for h = 0:cfg.H
            fprintf(fid, '%s,%s,%s,%d,%.6f,%.6f,%.6f\n', HEADLINE_Z, ests{e, 1}, vn{i}, h, ...
                    0.25 * ests{e, 2}(i, h + 1), 0.25 * ests{e, 3}(i, h + 1), 0.25 * ests{e, 4}(i, h + 1));
        end
    end
end
fclose(fid);
fprintf('  step A done in %.1f min -> %s, %s\n', toc(t_a) / 60, step_a_file, out_csv);
try
    figure('visible', 'off', 'Position', [50 50 1100 800]);
    hs = 0:cfg.H;
    for i = 1:K
        subplot(ceil(K / 2), 2, i);  hold on
        plot(hs, 0.25 * bz.theta(i, :), 'r-');  plot(hs, 0.25 * blpf.theta_mean(i, :), 'b-');
        plot(hs, 0.25 * blpb.theta_mean(i, :), 'k-');
        if DO_POOLED, plot(hs, 0.25 * blpp.theta_mean(i, :), 'm-'); end
        plot(hs, 0.25 * lpiv.theta(i, :), 'g--');
        plot(hs, 0.25 * lpiv.ar_lo(i, :), 'g:');  plot(hs, 0.25 * lpiv.ar_hi(i, :), 'g:');
        plot(hs, 0.25 * blpb.lo(i, :), 'k:');  plot(hs, 0.25 * blpb.hi(i, :), 'k:');
        plot(hs, zeros(size(hs)), '-', 'Color', [0.6 0.6 0.6]);
        title(sprintf('%s: 25 bp %s (F_{eff} = %.1f)', vn{i}, strrep(HEADLINE_Z, '_', '\_'), lpiv.first_stage.F_eff));
        if i == 1
            if DO_POOLED, legend('BVAR', 'BLP-FMAR', 'BLP-block', 'BLP-pooled', 'LP-IV', 'AR set', 'Location', 'best');
            else, legend('BVAR', 'BLP-FMAR', 'BLP-block', 'LP-IV', 'AR set', 'Location', 'best'); end
        end
    end
    print(fullfile(P.results, sprintf('fig_%s_irf.png', tag)), '-dpng', '-r110');
catch err
    fprintf('  (plotting skipped: %s)\n', err.message);
end
end

% ==============================================================================
if DO_B
fprintf('\n[B] tau heatmap and reading protocol ...\n');
assert(exist(step_a_file, 'file') == 2, 'RUN_EMPIRICAL_IV step B needs %s (step A).', step_a_file);
L = load(step_a_file);
sets = {'block', L.blpb};
if isfield(L, 'blpp') && ~isempty(L.blpp), sets(end + 1, :) = {'pooled', L.blpp}; end
out_csv = fullfile(P.results, sprintf('%s_tau_heatmap.csv', tag));
fid = fopen(out_csv, 'w');
fprintf(fid, 'scale_estimator,equation,block,h,lambda_h,tau_mean,tau_median,p_tau_gt1\n');
for e = 1:size(sets, 1)
    est = sets{e, 2};
    for i = 1:K, for g = 1:K, for h = 1:size(est.tau_mean, 3)
        fprintf(fid, '%s,%s,%s,%d,%.6g,%.5f,%.5f,%.4f\n', sets{e, 1}, vn{i}, vn{g}, h, L.blpf.lambda(1, h), ...
                est.tau_mean(i, g, h), est.tau_med(i, g, h), est.p_tau_gt1(i, g, h));
    end, end, end
end
fclose(fid);
fprintf('  wrote %s\n', out_csv);
try
    for e = 1:size(sets, 1)
        est = sets{e, 2};
        figure('visible', 'off', 'Position', [50 50 1100 800]);
        for i = 1:K
            subplot(ceil(K / 2), 2, i);
            imagesc(log(squeeze(est.tau_mean(i, :, :))));  colorbar
            title(sprintf('%s eq %s: log tau (%s)', SYSTEM, vn{i}, sets{e, 1}), 'Interpreter', 'none');
            set(gca, 'YTick', 1:K, 'YTickLabel', vn);  xlabel('horizon h');
        end
        print(fullfile(P.results, sprintf('fig_%s_tau_%s.png', tag, sets{e, 1})), '-dpng', '-r110');
    end
catch err
    fprintf('  (plotting skipped: %s)\n', err.message);
end
nf = NULL_FILE;
if isempty(nf), nf = fullfile(P.results, sprintf('null_calibration_%s.mat', strrep(tag, '_quick', ''))); end
if exist(nf, 'file') ~= 2 && strcmp(SYSTEM, 'lev4') && cfg.p == 12
    alt = fullfile(P.results, 'null_calibration_p12_levels.mat');
    if exist(alt, 'file') == 2
        fprintf('  note: %s not found; using the legacy null of the same data, %s\n', nf, alt);  nf = alt;
    end
end
if exist(nf, 'file') == 2
    nl = load(nf);
    if isfield(nl, 'H') && nl.H < cfg.H
        fprintf('  note: null simulated at H = %d, run at H = %d: exact for the independent estimator on h <= %d, an approximation for the pooled one.\n', nl.H, cfg.H, h_early(2));
    end
    prot = ea_apply_protocol(sets, nl, vn, h_early, [], fullfile(P.results, sprintf('%s_tau_protocol.csv', tag)), ...
                             struct('design_key', dkey, 'force', false));
    save(fullfile(P.results, sprintf('%s_protocol.mat', tag)), 'prot', 'nf', 'dkey');
else
    fprintf(['  no null for this design yet (%s).  Run\n' ...
             '    run_null_calibration(struct(''p'', %d, ''null'', struct(''n_rep'', 200, ''dataset'', ''%s'', ''stem'', ''%s'')))\n' ...
             '  and rerun with only DO_B = true.\n'], nf, cfg.p, DATASET, strrep(tag, '_quick', ''));
end
end

% ==============================================================================
if DO_C
fprintf('\n[C] LP-IV table, headline + alternatives ...\n');
assert(exist(step_a_file, 'file') == 2, 'RUN_EMPIRICAL_IV step C needs %s (step A).', step_a_file);
L = load(step_a_file);
out_csv = fullfile(P.results, sprintf('%s_lpiv.csv', tag));
fid = fopen(out_csv, 'w');
fprintf(fid, 'instrument,F_eff,variable,h,theta,se,lo,hi,ar_lo,ar_hi,ar_unbounded,theta_la,lo_la,hi_la\n');
items = {HEADLINE_Z, L.lpiv};
for q = 1:numel(L.ALT), items(end + 1, :) = {L.ALT(q).name, L.ALT(q).lpiv}; end  %#ok<AGROW>
for q = 1:size(items, 1)
    lp = items{q, 2};
    for i = 1:K, for h = 0:cfg.H
        fprintf(fid, '%s,%.3f,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f\n', items{q, 1}, lp.first_stage.F_eff, vn{i}, h, ...
                lp.theta(i, h + 1), lp.se(i, h + 1), lp.lo(i, h + 1), lp.hi(i, h + 1), lp.ar_lo(i, h + 1), lp.ar_hi(i, h + 1), ...
                lp.ar_unbounded(i, h + 1), lp.theta_la(i, h + 1), lp.lo_la(i, h + 1), lp.hi_la(i, h + 1));
    end, end
    fprintf('  %-26s F_eff %6.2f  AR unbounded at %3d of %d cells\n', items{q, 1}, lp.first_stage.F_eff, ...
            sum(sum(lp.ar_unbounded(:, 2:end))), K * cfg.H);
end
fclose(fid);
fprintf('  wrote %s\n', out_csv);
end

% ==============================================================================
if DO_D
fprintf('\n[D] cross-p tau maps, p in {%s}, read through each design''s own null ...\n', num2str(PLIST));
spec = struct('label', {}, 'result_file', {}, 'null_file', {}, 'held_blocks', {});
for ip = 1:numel(PLIST)
    cfp = cfg;  cfp.p = PLIST(ip);  cfp.H = min(cfg.H, 12);   % the statistic lives on h <= 12
    if QUICK, cfp.H = cfg.H; end
    ptag = sprintf('iv_%s_p%d', SYSTEM, PLIST(ip));  if QUICK, ptag = [ptag '_quick']; end
    rf = fullfile(P.results, sprintf('%s_tau.mat', ptag));
    t_d = tic;
    rng(20260200 + PLIST(ip), 'twister');
    bv = estimate_bvar_niw(Y, cfp);
    bf = estimate_blp_fmar(Y, cfp, bv);
    bb = estimate_blp_blockadaptive(Y, cfp, bv, bf.lambda);
    bp = [];  if DO_POOLED, bp = estimate_blp_blockpooled(Y, cfp, bv, bf.lambda); end
    relp = ea_relevance_iv(Y, bv, z_head, cfp, struct('ym', ym, 'z_name', HEADLINE_Z, 'verbose', false));
    hep = [2, min(12, cfp.H)];
    S_ = struct('blpb', bb, 'blpp', bp, 'cfg', cfp, 'varnames', {vn}, 'h_early', hep, ...
                'lambda', bf.lambda(1, :), 'relp', relp, 'T', T, 'ptag', ptag, 'dkey', ...
                ea_design_key(ds, cfp, struct('h_early', hep)));
    save(rf, '-struct', 'S_');
    nfp = fullfile(P.results, sprintf('null_calibration_%s.mat', strrep(ptag, '_quick', '')));
    if exist(nfp, 'file') ~= 2 && strcmp(SYSTEM, 'lev4') && PLIST(ip) == 12 && exist(fullfile(P.results, 'null_calibration_p12_levels.mat'), 'file') == 2
        nfp = fullfile(P.results, 'null_calibration_p12_levels.mat');
    end
    fprintf('  p = %2d: BVAR lambda %.3f, headline F_eff %.2f, %.1f min; null %s\n', PLIST(ip), bv.lambda, relp.fs.F_eff, toc(t_d) / 60, tern(exist(nfp, 'file') == 2, 'found', 'MISSING'));
    if exist(nfp, 'file') == 2
        spec(end + 1) = struct('label', sprintf('%s p=%d', SYSTEM, PLIST(ip)), 'result_file', rf, 'null_file', nfp, 'held_blocks', []);  %#ok<SAGROW>
    end
end
if ~isempty(spec)
    ea_cross_p_table(spec, fullfile(P.results, sprintf('iv_%s_cross_p%s.csv', SYSTEM, tern(QUICK, '_quick', ''))));
else
    fprintf('  (no nulls found for these designs; run run_null_calibration with stem iv_%s_p<p> for each p)\n', SYSTEM);
end
end

% ==============================================================================
if DO_E
fprintf('\n[E] out-of-sample block ablation ...\n');
pf = fullfile(P.results, sprintf('%s_protocol.mat', tag));
cells = {};  ctrl = {};
if exist(pf, 'file') == 2
    Q = load(pf);  pr = Q.prot;
    [cells, ctrl] = cells_from_protocol(pr, vn);
end
if isempty(cells)
    assert(exist(step_a_file, 'file') == 2, 'RUN_EMPIRICAL_IV step E needs step A or B output.');
    L = load(step_a_file);
    tb = mean(L.blpb.tau_mean(:, :, h_early(1):h_early(2)), 3);
    [~, imax] = max(tb(:));  [i1, g1] = ind2sub(size(tb), imax);
    [~, imin] = min(tb(:));  [i0, g0] = ind2sub(size(tb), imin);
    cells = {vn{i1}, vn{g1}};  ctrl = {vn{i0}, vn{g0}};
    fprintf('  no protocol escapes available: ablating the top cell %s <- %s, control %s <- %s\n', vn{i1}, vn{g1}, vn{i0}, vn{g0});
end
o0 = find(ym == 12 * ORIGIN0_YM(1) + ORIGIN0_YM(2), 1);  if isempty(o0), o0 = round(0.6 * T); end
abl = run_block_ablation(struct('dataset', DATASET, 'cells', {cells}, 'control_cells', {ctrl}, 'p', cfg.p, ...
                               'H', min(12, cfg.H), 'h_eval', h_early, 'origin0', o0, 'n_origins', N_ORIGINS, ...
                               'out_stem', tag, 'quick', QUICK, ...
                               'cfg_over', struct('fmar', struct('isrw', ds.isrw, 'psi_floor', true, 'h1_mode', 'lp'), ...
                                                  'blocks', struct('fixed_tau', []))));  %#ok<NASGU>
end

fprintf('\nRUN_EMPIRICAL_IV finished (tag %s).\n', tag);

% -------------------------------------------------------------------------
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function b = slim_bvar(b)
if isfield(b, 'F_draws'), b = rmfield(b, 'F_draws'); end
end

function th = reproject(beta, b1n, K, H, kind)
% point IRF under another impact vector from stored posterior-mean blocks:
% 'fmar': beta is (m x K x H) with the y(t) block in rows 2..K+1;
% 'block': beta is (K x K x H) = beta_block(i, :, h).
th = zeros(K, H + 1);  th(:, 1) = b1n;
for h = 1:H
    if strcmp(kind, 'fmar'), G = beta(2:1 + K, :, h)'; else, G = beta(:, :, h); end
    th(:, h + 1) = G * b1n;
end
end

function [cells, ctrl] = cells_from_protocol(pr, vn)
% escaping cells (block estimator, q95) and the quietest free cell
cells = {};  ctrl = {};
if isfield(pr, 'block') && isfield(pr.block, 'escape_q95')
    E = pr.block.escape_q95;  R = pr.block.ratio_q95;
elseif isfield(pr, 'escape_q95')
    E = pr.escape_q95;  R = pr.ratio_q95;
else
    return
end
[ii, gg] = find(E);
for q = 1:numel(ii), cells(end + 1, :) = {vn{ii(q)}, vn{gg(q)}}; end  %#ok<AGROW>
R(~isfinite(R)) = Inf;
[~, imin] = min(R(:));  [i0, g0] = ind2sub(size(R), imin);
ctrl = {vn{i0}, vn{g0}};
end
