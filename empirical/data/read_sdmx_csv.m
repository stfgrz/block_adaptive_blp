function [ym, val] = read_sdmx_csv(path)
% READ_SDMX_CSV  Minimal reader for SDMX-CSV files from the ECB Data Portal
% ('format=csvdata') or the Eurostat dissemination API ('format=SDMX-CSV').
%
% Both formats are plain comma-separated files whose header row contains,
% among other dimension columns, TIME_PERIOD (monthly: 'YYYY-MM') and
% OBS_VALUE.  This reader locates those two columns by name, so it is
% robust to the (different) dimension layouts of the two portals.
%
% INPUTS
% ------
%   path : csv file path.
%
% OUTPUTS
% -------
%   ym  : (n x 1) month index, 12*year + month, sorted ascending.
%   val : (n x 1) observation values.
%
% NOTES
% -----
% * Base MATLAB / Octave only.  Handles quoted fields (Eurostat quotes
%   some labels); does NOT handle embedded newlines inside quotes (not
%   used by either portal).
% * Rows with non-monthly TIME_PERIOD or non-numeric OBS_VALUE are
%   skipped silently (e.g. annual aggregates, empty observations).

fid = fopen(path, 'r');
assert(fid > 0, 'read_sdmx_csv: cannot open %s', path);

header = fgetl(fid);
cols = split_csv_line(header);
i_time = find(strcmpi(strtrim(cols), 'TIME_PERIOD'), 1);
i_val  = find(strcmpi(strtrim(cols), 'OBS_VALUE'), 1);
assert(~isempty(i_time) && ~isempty(i_val), ...
       'read_sdmx_csv: %s has no TIME_PERIOD/OBS_VALUE columns.', path);

ym  = zeros(0, 1);
val = zeros(0, 1);
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    if isempty(line), continue; end
    f = split_csv_line(line);
    if numel(f) < max(i_time, i_val), continue; end
    t = strtrim(f{i_time});
    v = str2double(strtrim(f{i_val}));
    if numel(t) < 7 || t(5) ~= '-' || isnan(v), continue; end
    y = str2double(t(1:4));  m = str2double(t(6:7));
    if isnan(y) || isnan(m), continue; end
    ym(end + 1, 1)  = 12 * y + m;  %#ok<AGROW>
    val(end + 1, 1) = v;           %#ok<AGROW>
end
fclose(fid);

[ym, order] = sort(ym);
val = val(order);
assert(~isempty(ym), 'read_sdmx_csv: no monthly observations found in %s', path);
end

% -------------------------------------------------------------------------
function f = split_csv_line(line)
% Split a csv line on commas, respecting double quotes.
f = {};  cur = '';  inq = false;
for i = 1:numel(line)
    c = line(i);
    if c == '"'
        inq = ~inq;
    elseif c == ',' && ~inq
        f{end + 1} = cur;  %#ok<AGROW>
        cur = '';
    else
        cur(end + 1) = c;  %#ok<AGROW>
    end
end
f{end + 1} = cur;
end
