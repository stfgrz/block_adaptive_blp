function m = import_ois_daily(opts)
% IMPORT_OIS_DAILY  Daily 1-month OIS (EONIA swap, Refinitiv EUREON1M=) ->
% monthly policy-indicator series under several timing conventions.
%
% PURPOSE
% -------
% The v2 design (docs/DESIGN.md Sec. 3-4) pairs the 1-month OIS
% SURPRISE with a 1-month OIS LEVEL as policy indicator, and the monthly
% aggregation of the surprise must match how the level is dated.  This
% function turns a daily export into the monthly candidates and validates it.
%
% INPUT FILE (opts.daily_csv, default raw/ois1m_ea_daily.csv)
% -----------------------------------------------------------------------
% Any csv with one date column and one value column (percent per annum):
%   * header names: the date column is the first column whose name contains
%     'date', 'time' or 'timestamp' (else column 1); the value is the first
%     column whose name contains 'close', 'last', 'mid', 'value', 'price' or
%     'EUREON'; failing that the MID of a 'Bid' and an 'Ask'/'Offer' column
%     (a lone side is used with a note); failing that column 2;
%   * date formats: yyyy-mm-dd, yyyy/mm/dd, dd/mm/yyyy, dd.mm.yyyy,
%     dd-mmm-yyyy (an optional time part is ignored);
%   * non-numeric or empty values are skipped; the file may be in either
%     time order.
% Excel serial dates are NOT handled: export dates as text.
%
% OUTPUT / SAVED (opts.out_csv, default derived/ois1m_ea_monthly.csv, and the
% .mat next to it)
% -----------------------------------------------------------------------
%   m.ym, .year, .month           monthly grid covering the data
%   m.eom      last available daily quote of the month  (END-OF-MONTH level)
%   m.avg      mean of the daily quotes in the month     (MONTHLY-AVERAGE level)
%   m.first    first quote of the month
%   m.chg_eom  eom(t) - eom(t-1)
%   m.n_days   number of daily quotes in the month
%   m.synthetic  true if the file is the synthetic fixture (header contains
%                'SYNTHETIC'); such a file is never written to the real
%                derived name unless opts.allow_synthetic = true
%   m.meta     source path, coverage, checks
%
% CHECKS (loud)
% -------------
%   coverage of [need_y0m0, need_y1m1] (default 1999m12 .. 2019m12 so that a
%   2000m1 sample has a lagged end-of-month value); no month with fewer than
%   opts.min_days quotes (default 10); no gap longer than opts.max_gap_days
%   calendar days (default 10); daily changes not implausibly large
%   (|d| > 150 bp flagged); values between -2 and 8 percent.  The Oct-2019
%   EONIA re-basing (EONIA = euro STR + 8.5 bp from 2 Oct 2019) is a level
%   convention and needs no treatment inside the 2000-2019 sample; EONIA was
%   discontinued on 3 Jan 2022.
%
% NOTES
% -----
% Base MATLAB / Octave.  A provenance sidecar <daily_csv>.source.txt is
% written by the caller (ea_write_provenance) when the file is placed.

if nargin < 1, opts = struct(); end
if exist('ea_paths', 'file') ~= 2
    addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
P = ea_paths();
opts = sd(opts, 'daily_csv', fullfile(P.raw, 'ois1m_ea_daily.csv'));
opts = sd(opts, 'out_csv', fullfile(P.derived, 'ois1m_ea_monthly.csv'));
opts = sd(opts, 'need_y0m0', [1999 12]);
opts = sd(opts, 'need_y1m1', [2019 12]);
opts = sd(opts, 'min_days', 10);
opts = sd(opts, 'max_gap_days', 10);
opts = sd(opts, 'allow_synthetic', false);
opts = sd(opts, 'verbose', true);

assert(exist(opts.daily_csv, 'file') == 2, ...
       ['import_ois_daily: %s not found.\nPlace the daily EUREON1M= export there ' ...
        '(see empirical/docs/DATA.md).'], opts.daily_csv);

% --- read -------------------------------------------------------------------------
fid = fopen(opts.daily_csv, 'r');
hdr_line = fgetl(fid);
hdr = strtrim(regexp(hdr_line, '[,;\t]', 'split'));
delim = ',';
if numel(hdr) == 1
    if any(hdr_line == ';'), delim = ';'; elseif any(hdr_line == sprintf('\t')), delim = sprintf('\t'); end
    hdr = strtrim(regexp(hdr_line, delim, 'split'));
end
synthetic = ~isempty(strfind(upper(hdr_line), 'SYNTHETIC'));  %#ok<STREMP>
lh = lower(hdr);
c_date = find(~cellfun(@isempty, regexp(lh, 'date|time')), 1);
if isempty(c_date), c_date = 1; end
% value column: an explicit mid/close/last first; otherwise the average of a
% BID and an ASK column (Refinitiv exports often carry both); otherwise
% column 2.  A lone bid or ask is used with a note.
c_val = find(~cellfun(@isempty, regexp(lh, 'close|last|mid|value|price|eureon')) & (1:numel(lh)) ~= c_date, 1);
c_bid = find(~cellfun(@isempty, regexp(lh, 'bid')), 1);
c_ask = find(~cellfun(@isempty, regexp(lh, 'ask|offer')), 1);
value_rule = '';
if ~isempty(c_val)
    value_rule = sprintf('column %s', hdr{c_val});
elseif ~isempty(c_bid) && ~isempty(c_ask)
    value_rule = sprintf('mid = (%s + %s) / 2', hdr{c_bid}, hdr{c_ask});
elseif ~isempty(c_bid) || ~isempty(c_ask)
    c_val = [c_bid c_ask];  c_val = c_val(1);
    value_rule = sprintf('column %s (only one side of the quote available)', hdr{c_val});
else
    c_val = min(2, numel(hdr));  if c_val == c_date, c_val = c_date + 1; end
    value_rule = sprintf('column %s (positional fallback)', hdr{c_val});
end
use_mid = isempty(c_val);
dn = zeros(0, 1);  val = zeros(0, 1);  n_skip = 0;
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    if isempty(strtrim(l)), continue; end
    f = strtrim(regexp(l, delim, 'split'));
    if use_mid
        if numel(f) < max([c_date, c_bid, c_ask]), n_skip = n_skip + 1; continue; end
        vb = str2double(strrep(f{c_bid}, '"', ''));  va = str2double(strrep(f{c_ask}, '"', ''));
        if isnan(vb) && isnan(va), v = NaN; elseif isnan(vb), v = va; elseif isnan(va), v = vb; else, v = (vb + va) / 2; end
    else
        if numel(f) < max(c_date, c_val), n_skip = n_skip + 1; continue; end
        v = str2double(strrep(f{c_val}, '"', ''));
    end
    d = parse_date(f{c_date});
    if isnan(d) || isnan(v), n_skip = n_skip + 1; continue; end
    dn(end + 1, 1) = d;  val(end + 1, 1) = v;  %#ok<AGROW>
end
fclose(fid);
assert(numel(dn) > 100, 'import_ois_daily: only %d usable daily rows in %s.', numel(dn), opts.daily_csv);
[dn, o] = sort(dn);  val = val(o);
[dn, iu] = unique(dn, 'last');  val = val(iu);      % keep the last quote of a duplicated day

% --- checks -------------------------------------------------------------------------
dv = datevec(dn);
ym_d = 12 * dv(:, 1) + dv(:, 2);
gaps = diff(dn);
chk = struct();
chk.n_days = numel(dn);
chk.first = datestr(dn(1), 'yyyy-mm-dd');  chk.last = datestr(dn(end), 'yyyy-mm-dd');
chk.max_gap_days = max(gaps);
chk.n_big_moves = sum(abs(diff(val)) > 1.5);
chk.min_val = min(val);  chk.max_val = max(val);
chk.n_skipped_rows = n_skip;

% --- monthly aggregation ---------------------------------------------------------------
ym = (ym_d(1):ym_d(end))';  n = numel(ym);
eom = nan(n, 1);  avg = nan(n, 1);  firstq = nan(n, 1);  nd = zeros(n, 1);
for k = 1:n
    idx = ym_d == ym(k);
    nd(k) = sum(idx);
    if nd(k) > 0
        vv = val(idx);
        eom(k) = vv(end);  firstq(k) = vv(1);  avg(k) = mean(vv);
    end
end
need0 = 12 * opts.need_y0m0(1) + opts.need_y0m0(2);
need1 = 12 * opts.need_y1m1(1) + opts.need_y1m1(2);
inwin = ym >= need0 & ym <= need1;
problems = {};
if ym(1) > need0 || ym(end) < need1
    problems{end + 1} = sprintf('coverage %s..%s does not span the needed %d-%02d..%d-%02d', ...
        chk.first, chk.last, opts.need_y0m0(1), opts.need_y0m0(2), opts.need_y1m1(1), opts.need_y1m1(2));
end
thin = find(inwin & nd < opts.min_days);
if ~isempty(thin)
    problems{end + 1} = sprintf('%d month(s) in the window with fewer than %d quotes (first: %d-%02d)', ...
        numel(thin), opts.min_days, floor((ym(thin(1)) - 1) / 12), ym(thin(1)) - 12 * floor((ym(thin(1)) - 1) / 12));
end
if chk.max_gap_days > opts.max_gap_days
    [~, ig] = max(gaps);
    problems{end + 1} = sprintf('gap of %d days after %s', chk.max_gap_days, datestr(dn(ig), 'yyyy-mm-dd'));
end
if chk.min_val < -2 || chk.max_val > 8
    problems{end + 1} = sprintf('values outside [-2, 8] percent (min %.3f, max %.3f)', chk.min_val, chk.max_val);
end
chk.problems = problems;
if opts.verbose
    fprintf('import_ois_daily: %s: value = %s; %d daily quotes %s..%s, %d skipped rows; max gap %d days; %d moves > 150 bp; range [%.3f, %.3f]%s\n', ...
            opts.daily_csv, value_rule, chk.n_days, chk.first, chk.last, n_skip, chk.max_gap_days, chk.n_big_moves, ...
            chk.min_val, chk.max_val, tern(synthetic, '  [SYNTHETIC FIXTURE]', ''));
    for i = 1:numel(problems), fprintf('  PROBLEM: %s\n', problems{i}); end
    if isempty(problems), fprintf('  coverage and plausibility checks passed\n'); end
end
assert(isempty(problems) || synthetic, ...
       'import_ois_daily: the daily file failed %d check(s); see above.', numel(problems));

m = struct('ym', ym, 'year', floor((ym - 1) / 12), 'month', ym - 12 * floor((ym - 1) / 12), ...
           'eom', eom, 'avg', avg, 'first', firstq, 'chg_eom', [NaN; diff(eom)], 'n_days', nd, ...
           'synthetic', synthetic);
m.meta = struct('daily_csv', opts.daily_csv, 'checks', chk, 'units', 'percent per annum', ...
                'built_at', datestr(now, 'yyyy-mm-dd HH:MM:SS'), 'date_col', hdr{c_date}, 'value_rule', value_rule);  %#ok<TNOW1,DATST>

% --- write ---------------------------------------------------------------------------------
if synthetic && ~opts.allow_synthetic
    if opts.verbose, fprintf('  synthetic fixture: not written to %s (pass allow_synthetic to override)\n', opts.out_csv); end
    return
end
outdir = fileparts(opts.out_csv);
if ~isempty(outdir) && exist(outdir, 'dir') ~= 7, mkdir(outdir); end
fid = fopen(opts.out_csv, 'w');
assert(fid > 0, 'import_ois_daily: cannot write %s', opts.out_csv);
fprintf(fid, 'year,month,eom,avg,first,chg_eom,n_days%s\n', tern(synthetic, ',SYNTHETIC', ''));
for k = 1:n
    fprintf(fid, '%d,%d,%.6f,%.6f,%.6f,%.6f,%d%s\n', m.year(k), m.month(k), eom(k), avg(k), firstq(k), ...
            m.chg_eom(k), nd(k), tern(synthetic, ',1', ''));
end
fclose(fid);
[pth, nm] = fileparts(opts.out_csv);
save(fullfile(pth, [nm '.mat']), '-struct', 'm');
if opts.verbose, fprintf('  wrote %s and %s.mat (%d months)\n', opts.out_csv, fullfile(pth, nm), n); end
end

% ==========================================================================================
function d = parse_date(s)
s = strtrim(strrep(s, '"', ''));
s = regexprep(s, '[T ].*$', '');                         % drop a time part
d = NaN;
try
    if ~isempty(regexp(s, '^\d{4}-\d{2}-\d{2}$', 'once'))
        d = datenum(s, 'yyyy-mm-dd');
    elseif ~isempty(regexp(s, '^\d{4}/\d{2}/\d{2}$', 'once'))
        d = datenum(s, 'yyyy/mm/dd');
    elseif ~isempty(regexp(s, '^\d{2}/\d{2}/\d{4}$', 'once'))
        d = datenum(s, 'dd/mm/yyyy');
    elseif ~isempty(regexp(s, '^\d{2}\.\d{2}\.\d{4}$', 'once'))
        d = datenum(s, 'dd.mm.yyyy');
    elseif ~isempty(regexp(s, '^\d{2}-[A-Za-z]{3}-\d{4}$', 'once'))
        d = datenum(s, 'dd-mmm-yyyy');
    elseif ~isempty(regexp(s, '^\d{1,2}/\d{1,2}/\d{4}$', 'once'))
        d = datenum(s, 'dd/mm/yyyy');
    end
catch
    d = NaN;
end
end

function s = sd(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
