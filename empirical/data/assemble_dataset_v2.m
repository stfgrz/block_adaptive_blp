function ds = assemble_dataset_v2(opts)
% ASSEMBLE_DATASET_V2  Estimation-ready dataset for the v2 design: a
% K-variable state vector WITHOUT the surprise, plus the monthly external
% instruments carried alongside (ds.Z), plus pre-event information.
%
% SYSTEMS (opts.system; or give opts.vars explicitly)
% ---------------------------------------------------------------------------
%   lev4      i1y, ip, hicp, stoxx                    2000m1-2019m12  (today)
%   lev4_yoy  i1y, ip, hicp_yoy, stoxx                2000m1-2019m12  (today)
%   ois4      ois1m, ip, hicp_sa, stoxx               2000m1-2019m12  (needs
%             the daily EUREON1M= file and hicp_sa_ea.csv)
%   ois6      ois1m, ip, hicp_sa, loans, lend_spread, stoxx   2003m1-2019m12
%   The first variable is the POLICY INDICATOR (cfg.shock_var = 1 downstream).
%
% VARIABLE REGISTRY (code : source file : transform : prior centre isrw)
%   i1y        rate1y_ea.csv     level, percent (12M Euribor, monthly avg)   1
%   ois1m      = ois1m_eom       end-of-month 1M OIS level, percent          1
%   ois1m_eom / ois1m_avg / ois1m_first   from derived/ois1m_ea_monthly.mat
%                                (import_ois_daily)                          1
%   eonia_avg  eonia_ea.csv      level, percent                              1
%   euribor1m  euribor1m_ea.csv  level, percent                              1
%   ip         ip_ea.csv         100 log                                     1
%   hicp       hicp_ea.csv       100 log, NSA                                1
%   hicp_sa    hicp_sa_ea.csv    100 log, ECB SA (ICP.M.U2.Y.000000.3.INX)   1
%   hicp_yoy   hicp_ea.csv       100 (log P_t - log P_{t-12})                1*
%   stoxx      stoxx50_ea.csv    100 log                                     1
%   loans      loans_nfc_ea.csv  100 log (adjusted loans to NFCs, notional
%                                stocks index, SA; from 2003m1)              1
%   ccb        ccb_nfc_ea.csv    level, percent (composite cost of
%                                borrowing, NFCs; from 2003m1)               1
%   lend_spread ccb minus the policy indicator (opts.spread_ref: 'indicator'
%                                (default) | 'i1y' | 'eonia_avg')            0
%   * isrw for hicp_yoy is a documented choice (persistent YoY inflation);
%     override any entry with opts.isrw_override = struct('hicp_yoy', 0).
%
% INSTRUMENTS (ds.Z, basis points, aligned with ds.ym)
% ---------------------------------------------------------------------------
%   every column of derived/instruments_monthly.mat (build_instrument_series;
%   built here if missing), each also in a '_bs' version orthogonalised on
%   pre-event public information over event months (Bauer-Swanson logic,
%   ea_pre_event_info; zero in non-event months), and any external series
%   given in opts.external = {struct('path', .., 'column', .., 'name', ..,
%   'scale', ..)} (import_external_instrument).
%
% SAMPLE.  opts.y0m0 (default [2000 1]) and opts.y1m1 ([2019 12]).  If a
% chosen series starts later, the start is moved to the latest first
% observation with a warning (opts.auto_shift = true, default) or an error.
%
% OUTPUT / SAVED (opts.out_mat, default data/ea_dataset_v2_<system>.mat)
%   ds.Y, .varnames, .ym, .isrw, .Z, .znames, .zunits, .pred, .prednames,
%   .n_gc, .n_sp_1m (event counts), .system, .instrument_free = true,
%   .shock_variant = 'none (external instruments in ds.Z)', .synthetic,
%   .meta (opts, sources, built_at).
%
% NOTES
% -----
% * Inputs are validated with ea_check_series exactly as in v1.
% * A SYNTHETIC OIS file (import_ois_daily flags it) makes ds.synthetic true
%   and the drivers refuse the dataset; only the interface test uses it.

if nargin < 1, opts = struct(); end
if exist('ea_paths', 'file') ~= 2
    addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
P = ea_paths();
opts = sd(opts, 'system', 'lev4');
systems = struct('lev4', {{'i1y', 'ip', 'hicp', 'stoxx'}}, 'lev4_yoy', {{'i1y', 'ip', 'hicp_yoy', 'stoxx'}}, ...
                 'ois4', {{'ois1m', 'ip', 'hicp_sa', 'stoxx'}}, ...
                 'ois6', {{'ois1m', 'ip', 'hicp_sa', 'loans', 'lend_spread', 'stoxx'}});
if ~isfield(opts, 'vars') || isempty(opts.vars)
    assert(isfield(systems, opts.system), 'assemble_dataset_v2: unknown system %s (have %s)', ...
           opts.system, strjoin(fieldnames(systems)', ', '));
    opts.vars = systems.(opts.system);
end
opts = sd(opts, 'y0m0', [2000 1]);  opts = sd(opts, 'y1m1', [2019 12]);
opts = sd(opts, 'raw_dir', P.raw);  opts = sd(opts, 'auto_shift', true);
opts = sd(opts, 'spread_ref', 'indicator');
opts = sd(opts, 'instruments_mat', fullfile(P.derived, 'instruments_monthly.mat'));
opts = sd(opts, 'ois_monthly_mat', fullfile(P.derived, 'ois1m_ea_monthly.mat'));
opts = sd(opts, 'external', {});
opts = sd(opts, 'out_mat', fullfile(P.data, sprintf('ea_dataset_v2_%s.mat', opts.system)));
opts = sd(opts, 'skip_plausibility', false);
if ~isfield(opts, 'isrw_override'), opts.isrw_override = struct(); end
K = numel(opts.vars);
ym0 = 12 * opts.y0m0(1) + opts.y0m0(2);  ym1 = 12 * opts.y1m1(1) + opts.y1m1(2);

% --- registry --------------------------------------------------------------------
reg = struct();
reg.i1y       = r('rate1y_ea.csv',   'level',  1);
reg.eonia_avg = r('eonia_ea.csv',    'level',  1);
reg.euribor1m = r('euribor1m_ea.csv','level',  1);
reg.ip        = r('ip_ea.csv',       'log100', 1);
reg.hicp      = r('hicp_ea.csv',     'log100', 1);
reg.hicp_sa   = r('hicp_sa_ea.csv',  'log100', 1);
reg.hicp_yoy  = r('hicp_ea.csv',     'yoy',    1);
reg.stoxx     = r('stoxx50_ea.csv',  'log100', 1);
reg.loans     = r('loans_nfc_ea.csv','log100', 1);
reg.ccb       = r('ccb_nfc_ea.csv',  'level',  1);
reg.lend_spread = r('ccb_nfc_ea.csv', 'spread', 0);
for f = {'ois1m', 'ois1m_eom', 'ois1m_avg', 'ois1m_first'}
    reg.(f{1}) = r('ois1m_ea_daily.csv', 'ois', 1);
end

% --- load every raw series once (validated on its own coverage) ------------------
series = struct();  synthetic = false;  sources = struct();
for k = 1:K
    v = opts.vars{k};
    assert(isfield(reg, v), 'assemble_dataset_v2: unknown variable code %s', v);
    e = reg.(v);
    key = strrep(strrep(e.file, '.csv', ''), '.', '_');
    if isfield(series, key), continue; end
    switch e.kind
        case 'ois'
            m = load_ois(opts, P);
            synthetic = synthetic || (isfield(m, 'synthetic') && m.synthetic);
            series.(key) = struct('ym', m.ym, 'eom', m.eom, 'avg', m.avg, 'first', m.first);
            sources.(key) = sprintf('%s (import_ois_daily; %s)', m.meta.daily_csv, tern(synthetic, 'SYNTHETIC FIXTURE', 'real'));
        otherwise
            fpath = fullfile(opts.raw_dir, e.file);
            assert(exist(fpath, 'file') == 2, ...
                   ['assemble_dataset_v2: %s is missing.\n' ...
                    '  legacy four: fetch_outcome_data();  v2 extras: ea_fetch_v2_series();\n' ...
                    '  1M OIS: see empirical/docs/MONDAY_DATA_CHECKLIST.md'], fpath);
            [yms, vs] = read_sdmx_csv(fpath);
            series.(key) = struct('ym', yms, 'val', vs);
            sources.(key) = read_sidecar(fpath);
    end
end

% --- effective sample: latest start among the chosen series ----------------------
starts = zeros(K, 1);  needs_lag12 = false(K, 1);
for k = 1:K
    e = reg.(opts.vars{k});  key = strrep(strrep(e.file, '.csv', ''), '.', '_');
    yms = series.(key).ym;
    starts(k) = yms(1);
    if strcmp(e.kind, 'yoy'), needs_lag12(k) = true;  starts(k) = yms(1) + 12; end
end
ym0_eff = max([ym0; starts]);
if ym0_eff > ym0
    msg = sprintf('assemble_dataset_v2: %s starts %d-%02d; sample start moved from %d-%02d to %d-%02d', ...
                  opts.vars{find(starts == ym0_eff, 1)}, floor((ym0_eff - 1) / 12), ym0_eff - 12 * floor((ym0_eff - 1) / 12), ...
                  opts.y0m0(1), opts.y0m0(2), floor((ym0_eff - 1) / 12), ym0_eff - 12 * floor((ym0_eff - 1) / 12));
    assert(opts.auto_shift, '%s (set opts.auto_shift = true to accept)', msg);
    warning('assemble_dataset_v2:shift', '%s', msg);
    ym0 = ym0_eff;
end
ym = (ym0:ym1)';  T = numel(ym);

% --- build Y ---------------------------------------------------------------------------
Y = nan(T, K);  isrw = ones(1, K);  labels = cell(1, K);
for k = 1:K
    v = opts.vars{k};  e = reg.(v);  key = strrep(strrep(e.file, '.csv', ''), '.', '_');
    S = series.(key);
    switch e.kind
        case 'level'
            validate(v, S.ym, S.val, ym0, ym1, opts);  Y(:, k) = pick(S.ym, S.val, ym);
            labels{k} = 'percent';
        case 'log100'
            validate(v, S.ym, S.val, ym0, ym1, opts);  Y(:, k) = 100 * log(pick(S.ym, S.val, ym));
            labels{k} = '100 log';
        case 'yoy'
            validate(v, S.ym, S.val, ym0 - 12, ym1, opts);
            lp = 100 * log(S.val);
            yoy = nan(size(lp));  yoy(13:end) = lp(13:end) - lp(1:end - 12);
            Y(:, k) = pick(S.ym, yoy, ym);  labels{k} = '100 log 12-month change';
        case 'ois'
            fld = 'eom';
            if strcmp(v, 'ois1m_avg'), fld = 'avg'; elseif strcmp(v, 'ois1m_first'), fld = 'first'; end
            validate(v, S.ym, S.(fld), ym0, ym1, struct('skip_plausibility', true));
            Y(:, k) = pick(S.ym, S.(fld), ym);  labels{k} = sprintf('percent (1M OIS, %s)', fld);
        case 'spread'
            validate(v, S.ym, S.val, ym0, ym1, opts);
            ccb = pick(S.ym, S.val, ym);
            switch opts.spread_ref
                case 'indicator', ref = [];      % filled after the loop (needs column 1)
                otherwise
                    e2 = reg.(opts.spread_ref);  key2 = strrep(strrep(e2.file, '.csv', ''), '.', '_');
                    if ~isfield(series, key2)
                        [yms, vs] = read_sdmx_csv(fullfile(opts.raw_dir, e2.file));  series.(key2) = struct('ym', yms, 'val', vs);
                    end
                    ref = pick(series.(key2).ym, series.(key2).val, ym);
            end
            if isempty(ref), Y(:, k) = ccb; spread_pending = k; else, Y(:, k) = ccb - ref; end  %#ok<NASGU>
            labels{k} = sprintf('percentage points (ccb NFC minus %s)', opts.spread_ref);
    end
    isrw(k) = e.isrw;
    if isfield(opts.isrw_override, v), isrw(k) = opts.isrw_override.(v); end
end
if exist('spread_pending', 'var')
    Y(:, spread_pending) = Y(:, spread_pending) - Y(:, 1);
end
bad = find(any(isnan(Y), 2));
if ~isempty(bad)
    error('assemble_dataset_v2: %d months with missing data (first: %d-%02d).', ...
          numel(bad), floor((ym(bad(1)) - 1) / 12), ym(bad(1)) - 12 * floor((ym(bad(1)) - 1) / 12));
end

% --- instruments -----------------------------------------------------------------------------
if exist(opts.instruments_mat, 'file') ~= 2
    fprintf('assemble_dataset_v2: building the instrument panel ...\n');
    build_instrument_series();
end
I = load(opts.instruments_mat);
[tf, loc] = ismember(ym, I.ym);
assert(all(tf), 'assemble_dataset_v2: the instrument panel does not cover %d-%02d..%d-%02d', ...
       opts.y0m0(1), opts.y0m0(2), opts.y1m1(1), opts.y1m1(2));
Z = I.Z(loc, :);  znames = I.names;
n_gc = I.n_gc(loc);  n_sp = I.n_sp_1m(loc);
for j = 1:numel(opts.external)
    ex = opts.external{j};
    eo = struct();  for f = {'column', 'scale', 'name', 'source'}, if isfield(ex, f{1}), eo.(f{1}) = ex.(f{1}); end, end
    E = import_external_instrument(ex.path, eo);
    v = nan(T, 1);  [tf, loc] = ismember(ym, E.ym);  v(tf) = E.value(loc(tf));
    Z(:, end + 1) = v;  znames{end + 1} = ['z_ext_' regexprep(E.name, '[^A-Za-z0-9]', '_')];  %#ok<AGROW>
    sources.(['external_' num2str(j)]) = sprintf('%s: %s (%s)', znames{end}, ex.path, E.source);
end
% pre-event information and the orthogonalised ('_bs') versions
[X, xn] = ea_pre_event_info(Y, opts.vars);
nz0 = numel(znames);
for j = 1:nz0
    z = Z(:, j);
    if strncmp(znames{j}, 'z_gc_', 5), ev = n_gc > 0; else, ev = n_gc > 0 | n_sp > 0; end
    o = ev & isfinite(z) & all(isfinite(X), 2);
    zb = zeros(T, 1);  zb(~isfinite(z)) = NaN;
    if sum(o) > size(X, 2) + 20
        Xo = [ones(sum(o), 1), X(o, :)];  b = Xo \ z(o);
        zb(o) = z(o) - Xo * b;
    else
        zb(:) = NaN;
    end
    Z(:, end + 1) = zb;  znames{end + 1} = [znames{j} '_bs'];  %#ok<AGROW>
end

% --- pack -------------------------------------------------------------------------------------
ds = struct('Y', Y, 'varnames', {opts.vars}, 'ym', ym, 'isrw', isrw, 'units', {labels}, ...
            'Z', Z, 'znames', {znames}, 'zunits', 'basis points (external instruments; NaN = unavailable)', ...
            'pred', X, 'prednames', {xn}, 'n_gc', n_gc, 'n_sp_1m', n_sp, ...
            'system', opts.system, 'instrument_free', true, ...
            'shock_variant', 'none (external instruments in ds.Z)', 'synthetic', synthetic);
ds.meta = opts;  ds.meta.sources = sources;
ds.meta.built_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');  %#ok<TNOW1,DATST>
fprintf('assemble_dataset_v2: system %s, T = %d (%d-%02d..%d-%02d), K = %d [%s], isrw = [%s], %d instruments%s\n', ...
        opts.system, T, floor((ym(1) - 1) / 12), ym(1) - 12 * floor((ym(1) - 1) / 12), ...
        floor((ym(end) - 1) / 12), ym(end) - 12 * floor((ym(end) - 1) / 12), K, strjoin(opts.vars, ', '), ...
        num2str(isrw), numel(znames), tern(synthetic, '  [SYNTHETIC OIS: dataset refused by the drivers]', ''));
fprintf('  col means: ');  fprintf('%.2f ', mean(Y));  fprintf('\n  col stds : ');  fprintf('%.2f ', std(Y));  fprintf('\n');
outdir = fileparts(opts.out_mat);
if ~isempty(outdir) && exist(outdir, 'dir') ~= 7, mkdir(outdir); end
save(opts.out_mat, '-struct', 'ds');
fprintf('  saved %s\n', opts.out_mat);
end

% ===================================================================================================
function e = r(file, kind, isrw)
e = struct('file', file, 'kind', kind, 'isrw', isrw);
end

function m = load_ois(opts, P)
if exist(opts.ois_monthly_mat, 'file') == 2
    m = load(opts.ois_monthly_mat);
    return
end
daily = fullfile(opts.raw_dir, 'ois1m_ea_daily.csv');
assert(exist(daily, 'file') == 2, ...
       ['assemble_dataset_v2: no 1M OIS data: neither %s nor %s exists.\n' ...
        '  Place the daily EUREON1M= export as %s and run import_ois_daily(); ' ...
        'see empirical/docs/MONDAY_DATA_CHECKLIST.md'], opts.ois_monthly_mat, daily, daily);
m = import_ois_daily(struct('daily_csv', daily, 'out_csv', fullfile(P.derived, 'ois1m_ea_monthly.csv')));
end

function validate(name, yms, vs, ym0, ym1, opts)
chk = ea_check_series(name, yms, vs, ym0, ym1);
if isfield(opts, 'skip_plausibility') && opts.skip_plausibility
    keep = {};
    for k = 1:numel(chk.msgs)
        if isempty(strfind(chk.msgs{k}, 'LOOKS SYNTHETIC')), keep{end + 1} = chk.msgs{k}; end  %#ok<STREMP,AGROW>
    end
    chk.msgs = keep;  chk.ok = isempty(keep);
end
if ~chk.ok
    msg = sprintf('assemble_dataset_v2: %s failed validation:', name);
    for k = 1:numel(chk.msgs), msg = sprintf('%s\n  - %s', msg, chk.msgs{k}); end
    error('%s', msg);
end
end

function v = pick(yms, vs, ym)
[tf, loc] = ismember(ym, yms);
v = nan(numel(ym), 1);  v(tf) = vs(loc(tf));
end

function s = read_sidecar(fpath)
s = 'unrecorded';
sf = [fpath '.source.txt'];
if exist(sf, 'file') == 2
    fid = fopen(sf, 'r');  txt = '';
    while true
        l = fgetl(fid);  if ~ischar(l), break; end
        if isempty(txt), txt = l; else, txt = [txt ' | ' l]; end  %#ok<AGROW>
    end
    fclose(fid);  s = txt;
end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function s = sd(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
