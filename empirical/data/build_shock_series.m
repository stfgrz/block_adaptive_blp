function shocks = build_shock_series(opts)
% BUILD_SHOCK_SERIES  Monthly ECB monetary-policy surprise series from EA-EMPD.
%
% PURPOSE
% -------
% Reads data/raw/ea_empd_events.csv (extracted from the EA-EMPD workbook of
% Altavilla, Gurkaynak, Kind & Laeven, "Monetary Transmission with Frequent
% Policy Events") and aggregates event-level high-frequency surprises to a
% monthly panel, following the aggregation of AGKL (ECB WP 3157, Sec. 4.3):
% sum of within-month event surprises, zero for months with no event.
%
% SERIES PRODUCED (all in BASIS POINTS at this stage)
% ---------------------------------------------------
%   mps_gc_1y   : baseline. GC monetary-event-window (Event_type GC_ME)
%                 1-year OIS surprises, summed by month.
%   mps_gc_3m   : robustness. Same, 3-month OIS.
%   mps_gc_jk   : Jarocinski-Karadi "poor man's sign restriction": keep the
%                 GC_ME 1Y surprise only when the STOXX50E surprise has the
%                 OPPOSITE sign (pure policy news); zero otherwise.
%   mps_gc_info : the complement (same-sign events = "information" type).
%   mps_all_1y  : GC_ME plus Executive-Board/President speech surprises in
%                 the 1Y OIS.  Speech filters follow AGKL's baseline:
%                 regular trading days only (Non_regular_trading_day == 0)
%                 and inside 9:00-18:00 CEST (Outside_regular_trading_hours
%                 == 0; NOTE: the Notes sheet's description of this flag is
%                 inverted -- the column NAME is correct, verified against
%                 event timestamps).
%   mps_all_tgt : AGKL-style day-count-adjusted TARGET surprise composite
%                 for GC + speeches: 1M OIS scaled by 30/(30-m) when the
%                 next GC meeting is m < 30 days away, 2M OIS scaled by
%                 60/(60-m) when 30 <= m < 60 (AGKL eq. (3); "use the
%                 two-month OIS rate when there is no meeting in the next
%                 month").  OUR ADDITIONAL RULE (documented deviation):
%                 events whose adjustment factor exceeds opts.max_adj (= 6)
%                 are dropped from this series only, to avoid amplifying
%                 microstructure noise near contract boundaries.
%   n_gc, n_speech : monthly event counts.
%
% INPUTS (all optional, via opts struct)
% --------------------------------------
%   opts.csv_path : events csv (default 'data/raw/ea_empd_events.csv')
%   opts.y0, m0   : first output month (default 1999, 1)
%   opts.y1, m1   : last output month  (default 2025, 12)
%   opts.max_adj  : cap on the day-count adjustment factor (default 6)
%   opts.out_csv  : output csv (default 'data/derived/shocks_monthly.csv')
%   opts.out_mat  : output mat (default 'data/derived/shocks_monthly.mat')
%
% OUTPUT
% ------
%   shocks : struct with fields ym (n x 1, = 12*year + month), year, month,
%            one (n x 1) vector per series above, and .meta.
%
% SANITY REFERENCE (soft checks, computed from the workbook shipped with
% this module, window 2001m1-2019m12; a quarterly EA-EMPD update will move
% these slightly -- warnings, not errors):
%   GC_ME events 221; monthly mps_gc_1y std 4.29 bp, min -18.20, max 21.19,
%   23 zero months; JK keeps 114/221 events (51.6%), monthly std 3.52 bp;
%   speeches kept after filters 1596/2701; mps_all_1y std 5.12 bp, corr
%   with mps_gc_1y 0.856, 1 zero month; mps_all_tgt std 4.90 bp.
%
% NOTES
% -----
% * Base MATLAB / Octave only (textscan; no readtable, no toolboxes).
% * GC_ME events have Days_until_next_GC == 0, so the day-count adjustment
%   factor is exactly 1 for meetings: the adjustment only matters for
%   speeches, as in AGKL.
% * Sign convention: positive surprise = unexpected TIGHTENING (rate up).

if nargin < 1, opts = struct(); end
opts = set_default(opts, 'csv_path', fullfile('data', 'raw', 'ea_empd_events.csv'));
opts = set_default(opts, 'y0', 1999);  opts = set_default(opts, 'm0', 1);
opts = set_default(opts, 'y1', 2025);  opts = set_default(opts, 'm1', 12);
opts = set_default(opts, 'max_adj', 6);
opts = set_default(opts, 'out_csv', fullfile('data', 'derived', 'shocks_monthly.csv'));
opts = set_default(opts, 'out_mat', fullfile('data', 'derived', 'shocks_monthly.mat'));

% --- read the events csv -------------------------------------------------
fid = fopen(opts.csv_path, 'r');
assert(fid > 0, 'build_shock_series: cannot open %s', opts.csv_path);
header = fgetl(fid);  %#ok<NASGU>  % fixed column order, documented above
C = textscan(fid, '%s %s %f %f %f %f %f %f %f %f %f %f', ...
             'Delimiter', ',', 'EmptyValue', NaN, 'ReturnOnError', false);
fclose(fid);

dstr   = C{1};            % 'yyyy-mm-dd HH:MM:SS'
etype  = C{2};
dgc    = C{3};            % Days_until_next_GC
nonreg = C{4};            % Non_regular_trading_day (1 = weekend/holiday)
outhrs = C{5};            % Outside_regular_trading_hours (1 = outside 9-18)
ois1m  = C{6};   ois2m = C{7};   ois3m = C{8};   %#ok<NASGU>
ois6m  = C{9};   %#ok<NASGU>
ois1y  = C{10};  ois2y = C{11};  %#ok<NASGU>
stoxx  = C{12};

n_ev = numel(dstr);
ev_year  = zeros(n_ev, 1);  ev_month = zeros(n_ev, 1);
for i = 1:n_ev
    ev_year(i)  = str2double(dstr{i}(1:4));
    ev_month(i) = str2double(dstr{i}(6:7));
end
ev_ym = 12 * ev_year + ev_month;

is_gc = strcmp(etype, 'GC_ME');
is_sp = strcmp(etype, 'EB') | strcmp(etype, 'P');

% --- monthly grid ---------------------------------------------------------
ym0 = 12 * opts.y0 + opts.m0;
ym1 = 12 * opts.y1 + opts.m1;
ym  = (ym0:ym1)';
n   = numel(ym);
Z   = zeros(n, 1);
shocks = struct('ym', ym, 'year', floor((ym - 1) / 12), ...
                'month', ym - 12 * floor((ym - 1) / 12), ...
                'mps_gc_1y', Z, 'mps_gc_3m', Z, 'mps_gc_jk', Z, ...
                'mps_gc_info', Z, 'mps_all_1y', Z, 'mps_all_tgt', Z, ...
                'n_gc', Z, 'n_speech', Z);

% --- speech filter (AGKL baseline) ---------------------------------------
sp_ok = is_sp & nonreg == 0 & outhrs == 0;

% --- day-count-adjusted target surprise per event ------------------------
tgt = nan(n_ev, 1);
for i = 1:n_ev
    if ~(is_gc(i) || sp_ok(i)), continue; end
    m = dgc(i);
    if isnan(m), continue; end
    if m < 30
        f = 30 / (30 - m);  x = ois1m(i);
    elseif m < 60
        f = 60 / (60 - m);  x = ois2m(i);
    else
        continue                          % >1 meeting horizon: drop (rare)
    end
    if f <= opts.max_adj && ~isnan(x)
        tgt(i) = f * x;
    end
end

% --- aggregate ------------------------------------------------------------
ev_jk_kept  = false(n_ev, 1);      % event-level flags so the sanity report
ev_gc_class = false(n_ev, 1);      % can be restricted to the ref window
for i = 1:n_ev
    k = ev_ym(i) - ym0 + 1;
    if k < 1 || k > n, continue; end
    if is_gc(i)
        shocks.n_gc(k) = shocks.n_gc(k) + 1;
        if ~isnan(ois1y(i))
            shocks.mps_gc_1y(k)  = shocks.mps_gc_1y(k)  + ois1y(i);
            shocks.mps_all_1y(k) = shocks.mps_all_1y(k) + ois1y(i);
        end
        if ~isnan(ois3m(i))
            shocks.mps_gc_3m(k) = shocks.mps_gc_3m(k) + ois3m(i);
        end
        if ~isnan(ois1y(i)) && ~isnan(stoxx(i))
            ev_gc_class(i) = true;
            if ois1y(i) * stoxx(i) < 0        % opposite sign: policy news
                shocks.mps_gc_jk(k) = shocks.mps_gc_jk(k) + ois1y(i);
                ev_jk_kept(i) = true;
            elseif ois1y(i) * stoxx(i) > 0    % same sign: information type
                shocks.mps_gc_info(k) = shocks.mps_gc_info(k) + ois1y(i);
            else                              % exact zero: neutral -> policy
                shocks.mps_gc_jk(k) = shocks.mps_gc_jk(k) + ois1y(i);
                ev_jk_kept(i) = true;
            end
        end
    elseif sp_ok(i) && ~isnan(ois1y(i))
        % n_speech counts speeches that CONTRIBUTE to the 1Y series
        % (passed filters, non-missing 1Y surprise)
        shocks.n_speech(k) = shocks.n_speech(k) + 1;
        shocks.mps_all_1y(k) = shocks.mps_all_1y(k) + ois1y(i);
    end
    if ~isnan(tgt(i))
        shocks.mps_all_tgt(k) = shocks.mps_all_tgt(k) + tgt(i);
    end
end

% --- soft sanity report (window 2001m1-2019m12) ---------------------------
w  = shocks.ym >= 12 * 2001 + 1 & shocks.ym <= 12 * 2019 + 12;
we = ev_ym >= 12 * 2001 + 1 & ev_ym <= 12 * 2019 + 12;
fprintf('build_shock_series: %d events read; window 2001m1-2019m12:\n', n_ev);
fprintf('  GC_ME events %d (ref 221) | mps_gc_1y std %.2f bp (ref 4.29), zero months %d (ref ~23)\n', ...
        sum(shocks.n_gc(w)), std_(shocks.mps_gc_1y(w)), sum(shocks.mps_gc_1y(w) == 0));
fprintf('  JK kept %d of %d classifiable GC events (ref 114/221) | std %.2f bp (ref 3.52)\n', ...
        sum(ev_jk_kept & we), sum(ev_gc_class & we), std_(shocks.mps_gc_jk(w)));
fprintf('  speeches kept %d (ref 1596) | mps_all_1y std %.2f bp (ref 5.12), corr with GC %.3f (ref 0.856)\n', ...
        sum(shocks.n_speech(w)), std_(shocks.mps_all_1y(w)), ...
        corr_(shocks.mps_all_1y(w), shocks.mps_gc_1y(w)));
fprintf('  mps_all_tgt std %.2f bp (ref 4.90)\n', std_(shocks.mps_all_tgt(w)));

shocks.meta = opts;

% --- write ---------------------------------------------------------------
outdir = fileparts(opts.out_csv);
if ~isempty(outdir) && ~exist(outdir, 'dir'), mkdir(outdir); end
fid = fopen(opts.out_csv, 'w');
fprintf(fid, ['year,month,mps_gc_1y,mps_gc_3m,mps_gc_jk,mps_gc_info,' ...
              'mps_all_1y,mps_all_tgt,n_gc,n_speech\n']);
for k = 1:n
    fprintf(fid, '%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d\n', ...
            shocks.year(k), shocks.month(k), shocks.mps_gc_1y(k), ...
            shocks.mps_gc_3m(k), shocks.mps_gc_jk(k), shocks.mps_gc_info(k), ...
            shocks.mps_all_1y(k), shocks.mps_all_tgt(k), ...
            shocks.n_gc(k), shocks.n_speech(k));
end
fclose(fid);
save(opts.out_mat, '-struct', 'shocks');
fprintf('build_shock_series: wrote %s and %s\n', opts.out_csv, opts.out_mat);
end

% -------------------------------------------------------------------------
function s = set_default(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

function v = std_(x)
v = sqrt(sum((x - mean(x)).^2) / (numel(x) - 1));
end

function r = corr_(x, y)
x = x - mean(x);  y = y - mean(y);
r = (x' * y) / sqrt((x' * x) * (y' * y));
end
