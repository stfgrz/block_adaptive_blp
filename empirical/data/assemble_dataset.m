function ds = assemble_dataset(opts)
% ASSEMBLE_DATASET  Build the estimation-ready euro-area dataset.
%
% SYSTEM (K = 5, ordering = identification; see docs/CH7_DESIGN.md Sec. 2)
% -----------------------------------------------------------------------
%   1  mps    monthly monetary-policy surprise, PERCENTAGE POINTS (bp/100),
%             ordered FIRST: internal-instrument recursive identification
%             (Stock-Watson 2018; Plagborg-Moller-Wolf 2021; the same
%             design as AGKL Sec. 4.3).  cfg.shock_var = 1.
%   2  i1y    1-year nominal rate level, percent.
%   3  ip     100 * log industrial production index.
%   4  hicp   100 * log HICP index.
%   5  stoxx  100 * log EURO STOXX 50 index.
%
% SAMPLE
% ------
% Y spans opts.y0m0 .. opts.y1m1, DEFAULT 2000m1-2019m12 (T = 240): with
% cfg.p = 12 the first 12 rows are lag initialisation, so the effective
% estimation sample is 2001m1-2019m12, matching AGKL's January-2001 start.
% Keep Y fixed across the p in {2, 6, 12} dose-response runs so that only
% the prior's flexibility changes, not the data.
%
% INPUTS (opts, all optional)
% ---------------------------
%   opts.shock_variant : column of shocks_monthly to use as variable 1.
%                        default 'mps_gc_1y'.  Alternatives: 'mps_gc_3m',
%                        'mps_gc_jk', 'mps_gc_info', 'mps_all_1y',
%                        'mps_all_tgt'.
%   opts.y0m0 / y1m1   : [year month] window bounds, default [2000 1] /
%                        [2019 12].
%   opts.raw_dir       : default <repo>/empirical/data/raw.
%   opts.out_mat       : default <repo>/empirical/data/ea_dataset.mat
%                        (the variant name is appended automatically for
%                        non-baseline variants).
%   opts.skip_plausibility : default false.  The four outcome series are
%                        re-validated here (coverage + the synthetic-data
%                        fingerprints of ea_check_series) so a fabricated
%                        or short input cannot reach the estimators even
%                        if fetch_outcome_data was bypassed.
% All default paths are ABSOLUTE, so this runs from any directory.
%
% OUTPUT / SAVED FIELDS
% ---------------------
%   ds.Y (T x 5), ds.varnames, ds.ym (T x 1), ds.isrw (1 x 5) = [0 1 1 1 1]
%   (white-noise prior mean for the surprise, random-walk for the levels --
%   requires the vectorised isrw patch in estimate_bvar_niw, see
%   README_EMPIRICAL.md), ds.shock_variant, ds.meta.
%
% DIAGNOSTICS PRINTED
% -------------------
% Missing-value check; head/tail; and the first-stage relevance regression
% Delta i1y_t = a + b mps_t + e_t with heteroskedasticity-robust t.  AGKL
% report a first-stage F of about 4 for meetings-only surprises and about
% 13 when speeches are included (their Fig. 9, 1M OIS on both sides) --
% expect the same ordering here.

if nargin < 1, opts = struct(); end

% Self-sufficient on the path, whatever the working directory is.
if exist('ea_paths', 'file') ~= 2
    addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
P = ea_paths();

opts = set_default(opts, 'shock_variant', 'mps_gc_1y');
opts = set_default(opts, 'y0m0', [2000 1]);
opts = set_default(opts, 'y1m1', [2019 12]);
opts = set_default(opts, 'raw_dir', P.raw);
opts = set_default(opts, 'shocks_mat', fullfile(P.derived, 'shocks_monthly.mat'));
opts = set_default(opts, 'out_mat', P.dataset);
if ~isfield(opts, 'skip_plausibility'), opts.skip_plausibility = false; end

% --- load shocks ----------------------------------------------------------
assert(exist(opts.shocks_mat, 'file') > 0, ...
       'assemble_dataset: run build_shock_series first (%s missing).', opts.shocks_mat);
S = load(opts.shocks_mat);
assert(isfield(S, opts.shock_variant), ...
       'assemble_dataset: unknown shock variant %s', opts.shock_variant);

% --- align -----------------------------------------------------------------
ym0 = 12 * opts.y0m0(1) + opts.y0m0(2);
ym1 = 12 * opts.y1m1(1) + opts.y1m1(2);
ym  = (ym0:ym1)';
T   = numel(ym);

% --- load and VALIDATE outcomes --------------------------------------------
% Second line of defence after fetch_outcome_data: a series that does not
% cover the window, or that carries a synthetic-data fingerprint, must
% never reach the estimators.  See ea_check_series.m for the checks and
% empirical/data/raw/placeholder_rejected/README.md for why they exist.
[ym_ip, v_ip] = load_outcome(opts, 'ip_ea.csv',      ym0, ym1);
[ym_pi, v_pi] = load_outcome(opts, 'hicp_ea.csv',    ym0, ym1);
[ym_r,  v_r ] = load_outcome(opts, 'rate1y_ea.csv',  ym0, ym1);
[ym_sx, v_sx] = load_outcome(opts, 'stoxx50_ea.csv', ym0, ym1);

Y = nan(T, 5);
Y(:, 1) = pick(S.ym, S.(opts.shock_variant), ym) / 100;   % bp -> pp
Y(:, 2) = pick(ym_r,  v_r,  ym);                          % percent
Y(:, 3) = 100 * log(pick(ym_ip, v_ip, ym));
Y(:, 4) = 100 * log(pick(ym_pi, v_pi, ym));
Y(:, 5) = 100 * log(pick(ym_sx, v_sx, ym));

bad = find(any(isnan(Y), 2));
if ~isempty(bad)
    error('assemble_dataset: %d months with missing outcome data (first: %d-%02d). Check the raw csvs.', ...
          numel(bad), floor((ym(bad(1))-1)/12), ym(bad(1)) - 12*floor((ym(bad(1))-1)/12));
end

ds.Y = Y;
ds.varnames = {'mps', 'i1y', 'ip', 'hicp', 'stoxx'};
ds.ym = ym;
ds.isrw = [0 1 1 1 1];
ds.shock_variant = opts.shock_variant;
ds.meta = opts;

% --- diagnostics ------------------------------------------------------------
fprintf('assemble_dataset: T = %d months (%d-%02d to %d-%02d), K = 5, variant = %s\n', ...
        T, floor((ym(1)-1)/12), ym(1)-12*floor((ym(1)-1)/12), ...
        floor((ym(end)-1)/12), ym(end)-12*floor((ym(end)-1)/12), opts.shock_variant);
fprintf('  col means: '); fprintf('%.2f ', mean(Y)); fprintf('\n');
fprintf('  col stds : '); fprintf('%.2f ', std(Y));  fprintf('\n');

% First-stage relevance: Delta i1y_t = a + b * mps_t + e_t, both in pp,
% over the FULL assembled window (t = 2..T, i.e. 2000m2-2019m12 by
% default -- not the p = 12 effective estimation sample, which starts a
% year later; the relevance of the instrument is a property of the series,
% not of the lag order).
di = diff(Y(:, 2));  s = Y(2:end, 1);
X  = [ones(numel(s), 1), s];
b  = X \ di;
u  = di - X * b;
XtXi = inv(X' * X);                                     %#ok<MINV>
V  = XtXi * (X' * (X .* (u.^2 * ones(1, 2)))) * XtXi;   % EHW
tb = b(2) / sqrt(V(2, 2));
fprintf('  first stage: Delta i1y = %.3f + %.3f mps, robust t = %.2f, F = %.2f\n', ...
        b(1), b(2), tb, tb^2);
fprintf('  (AGKL benchmark: F ~ 4 meetings-only, ~ 13 with speeches; a weak-ish\n');
fprintf('   meetings-only first stage is expected and is itself a Ch. 7 talking point.)\n');

% --- save --------------------------------------------------------------------
out = opts.out_mat;
if ~strcmp(opts.shock_variant, 'mps_gc_1y')
    [pth, nm, ext] = fileparts(out);
    out = fullfile(pth, sprintf('%s_%s%s', nm, opts.shock_variant, ext));
end
outdir = fileparts(out);
if ~isempty(outdir) && exist(outdir, 'dir') ~= 7
    [okdir, msg] = mkdir(outdir);
    assert(okdir == 1, 'assemble_dataset: cannot create %s (%s)', outdir, msg);
end
save(out, '-struct', 'ds');
fprintf('assemble_dataset: saved %s\n', out);
end

% -------------------------------------------------------------------------
function s = set_default(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

function [ym_s, v_s] = load_outcome(opts, fname, ym0, ym1)
% Read one outcome csv and refuse it if it is unusable.
fpath = fullfile(opts.raw_dir, fname);
assert(exist(fpath, 'file') == 2, ...
       ['assemble_dataset: %s is missing.\n' ...
        'Run fetch_outcome_data() first; if the automatic download fails it ' ...
        'prints the exact manual download route for this series.'], fpath);
[ym_s, v_s] = read_sdmx_csv(fpath);
chk = ea_check_series(fname, ym_s, v_s, ym0, ym1);
if ~chk.ok && ~opts.skip_plausibility
    msg = sprintf('assemble_dataset: %s failed validation:', fname);
    for k = 1:numel(chk.msgs), msg = sprintf('%s\n  - %s', msg, chk.msgs{k}); end
    error('%s\n%s', msg, ...
          ['Fix the input (see fetch_outcome_data) rather than the check. ' ...
           'Pass opts.skip_plausibility = true only if you have inspected ' ...
           'the series and decided the fingerprint is a false alarm.']);
end
end

function v = pick(ym_src, val_src, ym_want)
% Align a source series onto the wanted month grid (NaN where absent).
v = nan(numel(ym_want), 1);
[tf, loc] = ismember(ym_want, ym_src);
v(tf) = val_src(loc(tf));
end
