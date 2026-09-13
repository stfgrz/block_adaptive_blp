function ext = import_external_instrument(path, opts)
% IMPORT_EXTERNAL_INSTRUMENT  Read an externally supplied MONTHLY shock or
% instrument series (e.g. Jarocinski's updated ECB shocks, an
% informationally robust series obtained from its authors).
%
% INPUT
% -----
%   path : csv.  Accepted layouts (header names matched case-insensitively):
%          * columns 'year' and 'month' (any order), or
%          * one column named 'date', 'ym', 'period', 'time' or 'month'
%            holding 'yyyy-mm', 'yyyy-mm-dd', 'yyyymm' or 'yyyyMmm';
%          plus the value column named by opts.column (default: the first
%          numeric column that is not a date part).  Jarocinski's
%          shocks_ecb_mpd_me_m.csv works with opts.column = 'MP_pm' (or
%          'MP_median', 'CBI_pm', 'CBI_median').
%   opts : .column, .scale (multiply values, e.g. 100 if the file is in
%          percentage points and basis points are wanted; default 1),
%          .name (label; default = column), .source (free text recorded in
%          the output), .fill_zero (default true: months inside the file's
%          span with no row are set to 0, as event-based series intend;
%          false leaves NaN).
%
% OUTPUT
% ------
%   ext.ym, .year, .month, .value (n x 1), .name, .source, .units_note.
%
% The series is an EXTERNAL INSTRUMENT: it is attached to the dataset by
% assemble_dataset_v2 (opts.external) and used only by ea_identify_proxy /
% ea_relevance_iv / estimate_lp_iv.  Nothing here validates its
% construction; cite its authors and their licence.

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'column'), opts.column = ''; end
if ~isfield(opts, 'scale') || isempty(opts.scale), opts.scale = 1; end
if ~isfield(opts, 'source'), opts.source = path; end
if ~isfield(opts, 'fill_zero'), opts.fill_zero = true; end
fid = fopen(path, 'r');
assert(fid > 0, 'import_external_instrument: cannot open %s', path);
hl = fgetl(fid);
delim = ',';  if ~any(hl == ',') && any(hl == ';'), delim = ';'; end
hdr = strtrim(regexp(hl, delim, 'split'));
hdr = strrep(hdr, '"', '');
lh = lower(hdr);
raw = {};
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    if isempty(strtrim(l)), continue; end
    raw{end + 1} = strrep(strtrim(regexp(l, delim, 'split')), '"', '');  %#ok<AGROW>
end
fclose(fid);
n = numel(raw);
cy = find(strcmp(lh, 'year'), 1);  cm = find(strcmp(lh, 'month'), 1);
ym = nan(n, 1);
if ~isempty(cy) && ~isempty(cm) && cy ~= cm
    for i = 1:n
        r = raw{i};
        if numel(r) >= max(cy, cm), ym(i) = 12 * str2double(r{cy}) + str2double(r{cm}); end
    end
    datecols = [cy cm];
else
    cd = find(ismember(lh, {'date', 'ym', 'period', 'time', 'month', 'dates'}), 1);
    assert(~isempty(cd), 'import_external_instrument: no year/month or date column in %s', path);
    for i = 1:n
        r = raw{i};
        if numel(r) < cd, continue; end
        s = r{cd};
        tok = regexp(s, '^(\d{4})[-/Mm]?(\d{1,2})', 'tokens', 'once');
        if ~isempty(tok), ym(i) = 12 * str2double(tok{1}) + str2double(tok{2}); end
    end
    datecols = cd;
end
if isempty(opts.column)
    cand = setdiff(1:numel(hdr), datecols);
    cv = [];
    for c = cand
        v = str2double(raw{1}{min(c, numel(raw{1}))});
        if ~isnan(v), cv = c; break; end
    end
    assert(~isempty(cv), 'import_external_instrument: no numeric value column found in %s', path);
else
    cv = find(strcmpi(hdr, opts.column), 1);
    assert(~isempty(cv), 'import_external_instrument: column %s not in %s (have: %s)', ...
           opts.column, path, strjoin(hdr, ', '));
end
val = nan(n, 1);
for i = 1:n
    r = raw{i};
    if numel(r) >= cv, val(i) = str2double(r{cv}); end
end
keep = isfinite(ym);
ym = ym(keep);  val = val(keep) * opts.scale;
[ym, o] = sort(ym);  val = val(o);
grid = (ym(1):ym(end))';
v = nan(numel(grid), 1);
if opts.fill_zero, v(:) = 0; end
[tf, loc] = ismember(ym, grid);
v(loc(tf)) = val(tf);
ext = struct('ym', grid, 'year', floor((grid - 1) / 12), 'month', grid - 12 * floor((grid - 1) / 12), ...
             'value', v, 'name', hdr{cv}, 'source', opts.source, ...
             'units_note', sprintf('as in the file, times %g', opts.scale));
if isfield(opts, 'name') && ~isempty(opts.name), ext.name = opts.name; end
fprintf('import_external_instrument: %s, column %s, %d months %d-%02d..%d-%02d, %d non-zero\n', ...
        path, hdr{cv}, numel(grid), ext.year(1), ext.month(1), ext.year(end), ext.month(end), sum(v ~= 0 & isfinite(v)));
end
