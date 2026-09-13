function out_csv = make_synthetic_ois_daily(opts)
% MAKE_SYNTHETIC_OIS_DAILY  A clearly-labelled SYNTHETIC daily 1M-OIS file,
% for INTERFACE TESTING of import_ois_daily / assemble_dataset_v2 only.
%
% The numbers are simulated (a step path plus AR(1) noise on business days
% 1999-01-04 .. 2019-12-31); they are not market data and no result computed
% from them is an empirical finding.  Safeguards: the value column is named
% SYNTHETIC_EUREON1M= so import_ois_daily stamps m.synthetic = true and
% refuses to write the real derived file; the default output path is
% empirical/data/SYNTHETIC_ois1m_daily.csv (gitignored pattern SYNTHETIC_*).
%
% opts.out_csv (default above), opts.seed (default 7). Returns the path.

if nargin < 1, opts = struct(); end
if exist('ea_paths', 'file') ~= 2
    addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
P = ea_paths();
if ~isfield(opts, 'out_csv') || isempty(opts.out_csv)
    opts.out_csv = fullfile(P.data, 'SYNTHETIC_ois1m_daily.csv');
end
if ~isfield(opts, 'seed'), opts.seed = 7; end
s = rng;  rng(opts.seed, 'twister');
d0 = datenum(1999, 1, 4);  d1 = datenum(2019, 12, 31);
days = (d0:d1)';  wd = weekday(days);  days = days(wd >= 2 & wd <= 6);
n = numel(days);
% policy-rate step path: a level that changes by +-25 bp at random dates
level = 3.0;  path = zeros(n, 1);
for i = 1:n
    if rand < 1 / 45, level = level + 0.25 * sign(randn); end
    level = min(max(level, -0.5), 5);
    path(i) = level;
end
e = zeros(n, 1);
for i = 2:n, e(i) = 0.85 * e(i - 1) + 0.02 * randn; end
v = path + e;
fid = fopen(opts.out_csv, 'w');
assert(fid > 0, 'make_synthetic_ois_daily: cannot write %s', opts.out_csv);
fprintf(fid, 'Date,SYNTHETIC_EUREON1M=\n');
for i = 1:n, fprintf(fid, '%s,%.4f\n', datestr(days(i), 'yyyy-mm-dd'), v(i)); end
fclose(fid);
rng(s);
out_csv = opts.out_csv;
end
