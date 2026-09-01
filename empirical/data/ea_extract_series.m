function n = ea_extract_series(in_path, out_path, filt)
% EA_EXTRACT_SERIES  Reduce a multi-series SDMX-CSV export to one series.
%
% PURPOSE
% -------
% The Eurostat Data Browser makes it very easy to export far more than you
% asked for: leaving the geo or the activity filter untouched yields one
% file containing dozens of countries and several NACE aggregates stacked
% on top of each other.  Such a file parses fine and has the right column
% names, but every TIME_PERIOD appears many times, so it is not a series.
% ea_check_series rejects it (duplicated TIME_PERIOD); this function is
% the fix, and avoids re-downloading.
%
% USAGE
% -----
%   % 1. See what is in the file (no output file written):
%   ea_extract_series('data/raw/ip_ea.csv');
%
%   % 2. Extract one series, overwriting the file in place:
%   ea_extract_series('data/raw/ip_ea.csv', 'data/raw/ip_ea.csv', ...
%                     struct('geo', 'EA20', 'nace_r2', 'B-D'));
%
% INPUTS
% ------
%   in_path  : the downloaded csv.
%   out_path : where to write the cleaned series.  May equal in_path (the
%              file is read fully before anything is written).  Omit, or
%              pass [], to run in INSPECT MODE: the distinct values of
%              every dimension column are printed and nothing is written.
%   filt     : struct whose FIELD NAMES are column names in the file and
%              whose values are the required entries, e.g.
%              struct('geo', 'EA20', 'nace_r2', 'B-D').  Only the columns
%              you name are constrained.  Case-sensitive, matched on the
%              CODE columns (geo = 'EA20', not 'Euro area').
%
% OUTPUT
% ------
%   n : number of observations written (or found, in inspect mode).
%   The output file has the three columns read_sdmx_csv needs:
%   freq,TIME_PERIOD,OBS_VALUE.
%
% NOTES
% -----
% * Base MATLAB / Octave only.
% * Rows with an empty OBS_VALUE are dropped (Eurostat writes those for
%   periods a country does not report).
% * If the filter still leaves more than one series, the error names the
%   columns that are still varying and their distinct values, so the next
%   attempt is a copy-paste away.

if nargin < 2, out_path = []; end
if nargin < 3, filt = struct(); end

fid = fopen(in_path, 'r');
assert(fid > 0, 'ea_extract_series: cannot open %s', in_path);
header = fgetl(fid);
assert(ischar(header), 'ea_extract_series: %s is empty.', in_path);
cols = split_csv_line(header);
ncol = numel(cols);
C = textscan(fid, repmat('%q', 1, ncol), 'Delimiter', ',', ...
             'ReturnOnError', false);
fclose(fid);

i_time = find(strcmpi(strtrim(cols), 'TIME_PERIOD'), 1);
i_val  = find(strcmpi(strtrim(cols), 'OBS_VALUE'), 1);
assert(~isempty(i_time) && ~isempty(i_val), ...
       'ea_extract_series: %s has no TIME_PERIOD/OBS_VALUE columns.', in_path);

nrow = numel(C{i_time});
keep = true(nrow, 1);

% --- apply the filter -----------------------------------------------------
fn = fieldnames(filt);
for k = 1:numel(fn)
    j = find(strcmp(strtrim(cols), fn{k}), 1);
    assert(~isempty(j), ...
           'ea_extract_series: %s has no column named "%s". Columns are: %s', ...
           in_path, fn{k}, strjoin(cols, ', '));
    keep = keep & strcmp(C{j}, filt.(fn{k}));
end
% drop rows without an observation
keep = keep & ~cellfun(@isempty, strtrim_cell(C{i_val}));

assert(any(keep), ...
       ['ea_extract_series: the filter matched no rows. Run ' ...
        'ea_extract_series(''%s'') with no other arguments to list the ' ...
        'available values.'], in_path);

% --- identify the columns that still vary --------------------------------
% Dimension columns are the ones that are neither TIME_PERIOD/OBS_VALUE nor
% a free-text label; a column still varying after the filter means the
% selection is not yet a single series.
varying = {};
for j = 1:ncol
    if j == i_time || j == i_val, continue; end
    u = unique(C{j}(keep));
    if numel(u) > 1
        varying{end + 1} = struct('name', cols{j}, 'vals', {u});  %#ok<AGROW>
    end
end

% --- inspect mode ---------------------------------------------------------
if isempty(out_path)
    fprintf('ea_extract_series: %s\n  %d rows, %d columns.\n', ...
            in_path, nrow, ncol);
    if isempty(varying)
        fprintf('  This is already a single series (%d observations).\n', sum(keep));
    else
        fprintf('  Columns still varying (constrain these with the filter):\n');
        for k = 1:numel(varying)
            v = varying{k}.vals;
            show = v(1:min(numel(v), 12));
            fprintf('    %-12s %d values: %s%s\n', varying{k}.name, numel(v), ...
                    strjoin(show', ', '), tern(numel(v) > 12, ' ...', ''));
        end
    end
    n = sum(keep);
    return
end

% --- must be a single series by now ---------------------------------------
tp = C{i_time}(keep);
if numel(unique(tp)) ~= numel(tp)
    msg = sprintf(['ea_extract_series: the filter still leaves %d rows for ' ...
                   '%d distinct periods, i.e. more than one series.\n' ...
                   'Constrain these columns too:\n'], numel(tp), numel(unique(tp)));
    for k = 1:numel(varying)
        v = varying{k}.vals;
        show = v(1:min(numel(v), 12));
        msg = sprintf('%s    %-12s : %s%s\n', msg, varying{k}.name, ...
                      strjoin(show', ', '), tern(numel(v) > 12, ' ...', ''));
    end
    error('%s', msg);
end

% --- write ----------------------------------------------------------------
val = str2double(strtrim_cell(C{i_val}(keep)));
ym  = zeros(numel(tp), 1);
ok  = true(numel(tp), 1);
for i = 1:numel(tp)
    t = strtrim(tp{i});
    if numel(t) < 7 || t(5) ~= '-' || isnan(val(i))
        ok(i) = false;  continue
    end
    ym(i) = 12 * str2double(t(1:4)) + str2double(t(6:7));
end
tp = tp(ok);  val = val(ok);  ym = ym(ok);
assert(~isempty(ym), 'ea_extract_series: no monthly observations survived.');
[~, ord] = sort(ym);
tp = tp(ord);  val = val(ord);

outdir = fileparts(out_path);
if ~isempty(outdir) && exist(outdir, 'dir') ~= 7
    [okdir, msg] = mkdir(outdir);
    assert(okdir == 1, 'ea_extract_series: cannot create %s (%s)', outdir, msg);
end
fo = fopen(out_path, 'w');
assert(fo > 0, 'ea_extract_series: cannot write %s', out_path);
fprintf(fo, 'freq,TIME_PERIOD,OBS_VALUE\n');
for i = 1:numel(tp)
    fprintf(fo, 'M,%s,%.6f\n', strtrim(tp{i}), val(i));
end
fclose(fo);

n = numel(tp);
fprintf('ea_extract_series: wrote %d observations (%s .. %s) to %s\n', ...
        n, strtrim(tp{1}), strtrim(tp{end}), out_path);
end

% -------------------------------------------------------------------------
function c = strtrim_cell(c)
for i = 1:numel(c), c{i} = strtrim(c{i}); end
end

function s = tern(cond, a, b)
if cond, s = a; else, s = b; end
end

function f = split_csv_line(line)
% Split a csv line on commas, respecting double quotes.
f = {};  cur = '';  inq = false;
for i = 1:numel(line)
    ch = line(i);
    if ch == '"'
        inq = ~inq;
    elseif ch == ',' && ~inq
        f{end + 1} = cur;  %#ok<AGROW>
        cur = '';
    else
        cur(end + 1) = ch;  %#ok<AGROW>
    end
end
f{end + 1} = cur;
end
