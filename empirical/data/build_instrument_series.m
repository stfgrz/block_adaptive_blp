function inst = build_instrument_series(opts)
% BUILD_INSTRUMENT_SERIES  Monthly external instruments from the FULL EA-EMPD.
%
% PURPOSE
% -------
% v2 of the shock construction (docs/CH7_REDESIGN.md Sec. 2-3).  Reads the
% 54-column extract of the EA-EMPD workbook (raw/ea_empd_events_full.csv;
% Altavilla, Gurkaynak, Kind & Laeven 2025, ECB WP 3157) and produces, for
% every combination below, a MONTHLY series in BASIS POINTS, zero in months
% without a contributing event (AGKL p. 23):
%
%   event set   gc  = Governing Council monetary-event window (GC_ME)
%               sp  = Executive Board / President speeches (EB, P) passing
%                     the AGKL baseline filters: regular trading day and
%                     inside 9:00-18:00 CEST (their exclusion of windows
%                     containing a high-relevance data release cannot be
%                     replicated here: no release calendar)
%               gcs = gc + sp
%   surprise    1m_adj = 30/(30 - m) * dOIS_1M, m = Days_until_next_GC
%                        (AGKL eqs. 1-3; m = 0 for meetings so the factor is
%                        1; speeches with m >= 30 are dropped because the
%                        contract then ends before the meeting; factor capped
%                        at opts.max_adj, default 6, events above it dropped)
%               1m_raw = dOIS_1M unadjusted;  3m = dOIS_3M;  1y = dOIS_1Y
%   info        (gc only) jk   = Jarocinski-Karadi "poor man's" sign
%                                restriction: keep the surprise when it and
%                                the STOXX50E window move have OPPOSITE
%                                signs (an exact zero counts as policy);
%                          info = the complement (same sign)
%   stance      (gc only) path1y = residual of dOIS_1Y on the 1m_adj
%                          surprise across GC events: a ROUGH forward-
%                          guidance proxy (EA-MPD rotation logic: factors
%                          2-3 orthogonal to the 1M OIS; ABGMR 2019), NOT
%                          the EA-MPD FG factor
%   aggregation sum       = within-month sum (matches an END-OF-MONTH
%                           indicator: Kilian 2024)
%               kilian    = day-weighted with carry-over (matches a MONTHLY-
%                           AVERAGE indicator whose daily observation on the
%                           event day already reflects the event, e.g. an OIS
%                           close or EONIA): an event on business day d of a
%                           month with D business days puts (D-d+1)/D of the
%                           surprise in month t and (d-1)/D in month t+1
%               kilianfix = the same for an indicator FIXED BEFORE the event
%                           each day (Euribor, 11:00 CET, ahead of a 13:45
%                           release): (D-d)/D in month t, d/D in month t+1.
%               opts.day_basis = 'calendar' uses calendar days instead of
%               business days in both.
%
% Column names: z_<set>_<surprise>[_jk|_info]_<agg>, agg in {sum, kilian,
% kilianfix}, e.g. z_gcs_1m_adj_sum, z_gc_1m_adj_jk_kilian, z_gc_path1y_sum.  Also n_gc, n_sp_1m, n_sp_1y
% (event counts per month).
%
% INPUTS (opts, all optional)
% ---------------------------
%   csv_path   raw/ea_empd_events_full.csv
%   y0, m0 / y1, m1   output grid, default 1999m1 .. 2025m12
%   max_adj    6        cap on the day-count factor
%   day_basis  'business' | 'calendar'
%   out_csv / out_mat   derived/instruments_monthly.{csv,mat}
%   out_events          derived/instruments_events.csv (event-level table)
%
% OUTPUT
% ------
%   inst.ym, .year, .month (n x 1); inst.names (1 x nz cell); inst.Z (n x nz)
%   in basis points; inst.n_gc, .n_sp_1m, .n_sp_1y; inst.meta (options,
%   counts, the path regression slope, reference statistics).
%
% NOTES
% -----
% * Base MATLAB / Octave.  The csv has quoted text fields (speech titles
%   with commas): a small quote-aware line splitter handles them.
% * Sign convention: positive = unexpected tightening.
% * The legacy builder (build_shock_series.m, 12-column extract, 1Y sums)
%   is untouched; the legacy baseline keeps using it.

if nargin < 1, opts = struct(); end
if exist('ea_paths', 'file') ~= 2
    addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
P = ea_paths();
opts = sd(opts, 'csv_path', fullfile(P.raw, 'ea_empd_events_full.csv'));
opts = sd(opts, 'y0', 1999);  opts = sd(opts, 'm0', 1);
opts = sd(opts, 'y1', 2025);  opts = sd(opts, 'm1', 12);
opts = sd(opts, 'max_adj', 6);
opts = sd(opts, 'day_basis', 'business');
opts = sd(opts, 'out_csv', fullfile(P.derived, 'instruments_monthly.csv'));
opts = sd(opts, 'out_mat', fullfile(P.derived, 'instruments_monthly.mat'));
opts = sd(opts, 'out_events', fullfile(P.derived, 'instruments_events.csv'));
assert(any(strcmp(opts.day_basis, {'business', 'calendar'})), ...
       'build_instrument_series: day_basis must be ''business'' or ''calendar''.');

% --- read the full extract ---------------------------------------------------
[hdr, rows] = read_quoted_csv(opts.csv_path);
col = @(name) find(strcmp(hdr, name), 1);
c_dt = col('Date_time');  c_et = col('Event_type');
c_nr = col('Non_regular_trading_day');  c_oh = col('Outside_regular_trading_hours');
c_m  = col('Days_until_next_GC');
c_1m = col('OIS_1M');  c_3m = col('OIS_3M');  c_1y = col('OIS_1Y');  c_sx = col('STOXX50E');
assert(~any(cellfun(@isempty, {c_dt, c_et, c_nr, c_oh, c_m, c_1m, c_3m, c_1y, c_sx})), ...
       'build_instrument_series: %s lacks one of the required columns.', opts.csv_path);
n_ev = numel(rows);
dt = cell(n_ev, 1);  et = cell(n_ev, 1);
nonreg = nan(n_ev, 1);  outhrs = nan(n_ev, 1);  m_next = nan(n_ev, 1);
o1m = nan(n_ev, 1);  o3m = nan(n_ev, 1);  o1y = nan(n_ev, 1);  sx = nan(n_ev, 1);
for i = 1:n_ev
    r = rows{i};
    dt{i} = r{c_dt};  et{i} = r{c_et};
    nonreg(i) = num(r{c_nr});  outhrs(i) = num(r{c_oh});  m_next(i) = num(r{c_m});
    o1m(i) = num(r{c_1m});  o3m(i) = num(r{c_3m});  o1y(i) = num(r{c_1y});  sx(i) = num(r{c_sx});
end
ev_y = zeros(n_ev, 1);  ev_mo = zeros(n_ev, 1);  ev_d = zeros(n_ev, 1);
for i = 1:n_ev
    ev_y(i)  = str2double(dt{i}(1:4));
    ev_mo(i) = str2double(dt{i}(6:7));
    ev_d(i)  = str2double(dt{i}(9:10));
end
ev_ym = 12 * ev_y + ev_mo;

is_gc = strcmp(et, 'GC_ME');
is_sp = (strcmp(et, 'EB') | strcmp(et, 'P')) & nonreg == 0 & outhrs == 0;

% --- event-level surprises ----------------------------------------------------
% day-count adjustment (AGKL eq. 3), meetings have m = 0 -> factor 1
fac = nan(n_ev, 1);
ok_m = ~isnan(m_next) & m_next >= 0 & m_next < 30;
fac(ok_m) = 30 ./ (30 - m_next(ok_m));
fac(fac > opts.max_adj) = NaN;                 % capped: event dropped for 1m_adj
s_1m_adj = fac .* o1m;
s_1m_raw = o1m;  s_3m = o3m;  s_1y = o1y;

% JK classification on the GC events (1m_adj vs STOXX; and 1y vs STOXX)
jk_1m = is_gc & ~isnan(s_1m_adj) & ~isnan(sx) & (s_1m_adj .* sx <= 0);
in_1m = is_gc & ~isnan(s_1m_adj) & ~isnan(sx) & (s_1m_adj .* sx > 0);
jk_1y = is_gc & ~isnan(s_1y) & ~isnan(sx) & (s_1y .* sx <= 0);
in_1y = is_gc & ~isnan(s_1y) & ~isnan(sx) & (s_1y .* sx > 0);

% path proxy: residual of the 1Y move on the adjusted 1M move, GC events
sel = is_gc & ~isnan(s_1m_adj) & ~isnan(s_1y);
Xp = [ones(sum(sel), 1), s_1m_adj(sel)];
bp = Xp \ s_1y(sel);
s_path = nan(n_ev, 1);
s_path(sel) = s_1y(sel) - Xp * bp;

% --- monthly grid and Kilian weights -----------------------------------------
ym0 = 12 * opts.y0 + opts.m0;  ym1 = 12 * opts.y1 + opts.m1;
ym = (ym0:ym1)';  n = numel(ym);
w_t  = nan(n_ev, 1);                           % share in month t, event-day obs reflects the event
w_tf = nan(n_ev, 1);                           % share in month t, indicator fixed before the event
for i = 1:n_ev
    if strcmp(opts.day_basis, 'business')
        days = datenum(ev_y(i), ev_mo(i), 1):datenum(ev_y(i), ev_mo(i), eomday(ev_y(i), ev_mo(i)));
        wd = weekday(days);  isbd = wd >= 2 & wd <= 6;
        D = sum(isbd);
        d_b = sum(isbd(1:ev_d(i)));           % business days up to and incl. the event day
        event_on_bd = isbd(ev_d(i));
    else
        D = eomday(ev_y(i), ev_mo(i));  d_b = ev_d(i);  event_on_bd = true;
    end
    if event_on_bd, w_t(i) = (D - d_b + 1) / D; else, w_t(i) = (D - d_b) / D; end
    w_tf(i) = (D - d_b) / D;
end

% --- assemble the columns -------------------------------------------------------
sets = struct('gc', is_gc, 'sp', is_sp, 'gcs', is_gc | is_sp);
setn = fieldnames(sets);
surp = struct('x1m_adj', s_1m_adj, 'x1m_raw', s_1m_raw, 'x3m', s_3m, 'x1y', s_1y);
surn = fieldnames(surp);
names = {};  Z = zeros(n, 0);
for a = 1:numel(setn)
    for b = 1:numel(surn)
        mask = sets.(setn{a});  v = surp.(surn{b});
        lab = surn{b}(2:end);                  % strip the leading x
        [zs, zk] = aggregate(v, mask, ev_ym, w_t, ym0, n);
        [~, zf]  = aggregate(v, mask, ev_ym, w_tf, ym0, n);
        names{end + 1} = sprintf('z_%s_%s_sum', setn{a}, lab);        Z(:, end + 1) = zs;  %#ok<AGROW>
        names{end + 1} = sprintf('z_%s_%s_kilian', setn{a}, lab);     Z(:, end + 1) = zk;  %#ok<AGROW>
        names{end + 1} = sprintf('z_%s_%s_kilianfix', setn{a}, lab);  Z(:, end + 1) = zf;  %#ok<AGROW>
    end
end
info_sets = {'1m_adj', s_1m_adj, jk_1m, in_1m; '1y', s_1y, jk_1y, in_1y};
for b = 1:size(info_sets, 1)
    for c = 1:2
        if c == 1, mask = info_sets{b, 3};  tagc = 'jk'; else, mask = info_sets{b, 4};  tagc = 'info'; end
        [zs, zk] = aggregate(info_sets{b, 2}, mask, ev_ym, w_t, ym0, n);
        [~, zf]  = aggregate(info_sets{b, 2}, mask, ev_ym, w_tf, ym0, n);
        names{end + 1} = sprintf('z_gc_%s_%s_sum', info_sets{b, 1}, tagc);        Z(:, end + 1) = zs;  %#ok<AGROW>
        names{end + 1} = sprintf('z_gc_%s_%s_kilian', info_sets{b, 1}, tagc);     Z(:, end + 1) = zk;  %#ok<AGROW>
        names{end + 1} = sprintf('z_gc_%s_%s_kilianfix', info_sets{b, 1}, tagc);  Z(:, end + 1) = zf;  %#ok<AGROW>
    end
end
[zs, zk] = aggregate(s_path, is_gc, ev_ym, w_t, ym0, n);
[~, zf]  = aggregate(s_path, is_gc, ev_ym, w_tf, ym0, n);
names{end + 1} = 'z_gc_path1y_sum';        Z(:, end + 1) = zs;
names{end + 1} = 'z_gc_path1y_kilian';     Z(:, end + 1) = zk;
names{end + 1} = 'z_gc_path1y_kilianfix';  Z(:, end + 1) = zf;

n_gc = zeros(n, 1);  n_sp_1m = zeros(n, 1);  n_sp_1y = zeros(n, 1);
for i = 1:n_ev
    k = ev_ym(i) - ym0 + 1;
    if k < 1 || k > n, continue; end
    if is_gc(i), n_gc(k) = n_gc(k) + 1; end
    if is_sp(i) && ~isnan(s_1m_adj(i)), n_sp_1m(k) = n_sp_1m(k) + 1; end
    if is_sp(i) && ~isnan(s_1y(i)),     n_sp_1y(k) = n_sp_1y(k) + 1; end
end

inst = struct('ym', ym, 'year', floor((ym - 1) / 12), 'month', ym - 12 * floor((ym - 1) / 12), ...
              'names', {names}, 'Z', Z, 'n_gc', n_gc, 'n_sp_1m', n_sp_1m, 'n_sp_1y', n_sp_1y);
inst.units = 'basis points; zero in months without a contributing event';
inst.meta = opts;
inst.meta.n_events_read = n_ev;
inst.meta.path_slope_1y_on_1m_adj = bp(2);
inst.meta.built_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');  %#ok<TNOW1,DATST>

% --- reference statistics, window 2001m1-2019m12 --------------------------------
w = ym >= 12 * 2001 + 1 & ym <= 12 * 2019 + 12;
we = ev_ym >= 12 * 2001 + 1 & ev_ym <= 12 * 2019 + 12;
ref = struct();
ref.gc_events = sum(is_gc & we);
ref.gc_events_with_1m = sum(is_gc & we & ~isnan(s_1m_adj));
ref.speeches_kept = sum(is_sp & we);
ref.speeches_with_1m_adj = sum(is_sp & we & ~isnan(s_1m_adj));
ref.speeches_dropped_m_ge_30 = sum(is_sp & we & ~isnan(o1m) & m_next >= 30);
ref.speeches_dropped_cap = sum(is_sp & we & ~isnan(o1m) & ok_m & 30 ./ (30 - m_next) > opts.max_adj);
ref.jk_kept_1m = sum(jk_1m & we);  ref.jk_classifiable_1m = sum((jk_1m | in_1m) & we);
fprintf('build_instrument_series: %d events read from %s\n', n_ev, opts.csv_path);
fprintf('  window 2001m1-2019m12: GC_ME %d (with 1M %d); speeches kept %d (with adjusted 1M %d; dropped m>=30: %d, cap: %d); JK keeps %d of %d\n', ...
        ref.gc_events, ref.gc_events_with_1m, ref.speeches_kept, ref.speeches_with_1m_adj, ...
        ref.speeches_dropped_m_ge_30, ref.speeches_dropped_cap, ref.jk_kept_1m, ref.jk_classifiable_1m);
fprintf('  path proxy: dOIS_1Y = %.2f + %.2f * dOIS_1M(adj) over %d GC events\n', bp(1), bp(2), sum(sel));
fprintf('  %-26s %8s %8s %6s\n', 'series', 'sd (bp)', 'zeros', 'corr');
zref = Z(:, strcmp(names, 'z_gc_1y_sum'));
for j = 1:numel(names)
    v = Z(w, j);  cc = corrcoef(v, zref(w));
    ref.sd.(names{j}) = std(v);  ref.zeros.(names{j}) = sum(v == 0);
    if any(strcmp(names{j}, {'z_gc_1m_adj_sum', 'z_gc_1m_adj_kilian', 'z_gcs_1m_adj_sum', ...
                              'z_gcs_1m_adj_kilian', 'z_gc_1y_sum', 'z_gcs_1y_sum', 'z_gc_1m_adj_jk_sum', 'z_gc_path1y_sum'}))
        fprintf('  %-26s %8.2f %8d %6.2f\n', names{j}, std(v), sum(v == 0), cc(1, 2));
    end
end
inst.meta.ref = ref;

% --- write ------------------------------------------------------------------------
outdir = fileparts(opts.out_csv);
if ~isempty(outdir) && exist(outdir, 'dir') ~= 7, mkdir(outdir); end
fid = fopen(opts.out_csv, 'w');
assert(fid > 0, 'build_instrument_series: cannot write %s', opts.out_csv);
fprintf(fid, 'year,month,n_gc,n_sp_1m,n_sp_1y');
fprintf(fid, ',%s', names{:});  fprintf(fid, '\n');
for k = 1:n
    fprintf(fid, '%d,%d,%d,%d,%d', inst.year(k), inst.month(k), n_gc(k), n_sp_1m(k), n_sp_1y(k));
    fprintf(fid, ',%.6f', Z(k, :));  fprintf(fid, '\n');
end
fclose(fid);
save(opts.out_mat, '-struct', 'inst');
% event-level table (diagnostics; which events contribute what)
fid = fopen(opts.out_events, 'w');
if fid > 0
    fprintf(fid, 'date_time,event_type,in_gc,in_sp,days_until_next_gc,adj_factor,s_1m_adj,s_1m_raw,s_3m,s_1y,stoxx,jk_1m,w_month_t,w_month_t_fix\n');
    for i = 1:n_ev
        if ~(is_gc(i) || is_sp(i)), continue; end
        fprintf(fid, '%s,%s,%d,%d,%g,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.4f\n', dt{i}, et{i}, is_gc(i), is_sp(i), ...
                m_next(i), fac(i), s_1m_adj(i), s_1m_raw(i), s_3m(i), s_1y(i), sx(i), jk_1m(i), w_t(i), w_tf(i));
    end
    fclose(fid);
end
fprintf('  wrote %s (%d months x %d series), %s, %s\n', opts.out_csv, n, numel(names), opts.out_mat, opts.out_events);
end

% =============================================================================
function [zs, zk] = aggregate(v, mask, ev_ym, w_t, ym0, n)
% within-month sum, and the Kilian day-weighted version with carry-over
zs = zeros(n, 1);  zk = zeros(n, 1);
for i = 1:numel(v)
    if ~mask(i) || isnan(v(i)), continue; end
    k = ev_ym(i) - ym0 + 1;
    if k >= 1 && k <= n
        zs(k) = zs(k) + v(i);
        zk(k) = zk(k) + w_t(i) * v(i);
    end
    if k + 1 >= 1 && k + 1 <= n
        zk(k + 1) = zk(k + 1) + (1 - w_t(i)) * v(i);
    end
end
end

function x = num(s)
if isempty(s), x = NaN; else, x = str2double(s); end
end

function s = sd(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

function [hdr, rows] = read_quoted_csv(path)
% Quote-aware csv reader for the EA-EMPD extract: fields that contain a
% comma are double-quoted (speech titles), so plain strsplit is wrong.
fid = fopen(path, 'r');
assert(fid > 0, 'build_instrument_series: cannot open %s', path);
first = fgetl(fid);
hdr = split_line(first);
rows = cell(6000, 1);  k = 0;
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    if isempty(l), continue; end
    k = k + 1;
    if k > numel(rows), rows{2 * k} = []; end
    rows{k} = split_line(l);
end
fclose(fid);
rows = rows(1:k);
nh = numel(hdr);
for i = 1:k
    if numel(rows{i}) < nh, rows{i}(end + 1:nh) = {''}; end
end
end

function f = split_line(s)
if ~any(s == '"')
    f = regexp(s, ',', 'split');
    return
end
f = {};  cur = '';  inq = false;  i = 1;  n = numel(s);
while i <= n
    ch = s(i);
    if inq
        if ch == '"'
            if i < n && s(i + 1) == '"', cur(end + 1) = '"';  i = i + 1; else, inq = false; end  %#ok<AGROW>
        else
            cur(end + 1) = ch;  %#ok<AGROW>
        end
    else
        if ch == '"'
            inq = true;
        elseif ch == ','
            f{end + 1} = cur;  cur = '';  %#ok<AGROW>
        else
            cur(end + 1) = ch;  %#ok<AGROW>
        end
    end
    i = i + 1;
end
f{end + 1} = cur;
end
