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
%   opts.raw_dir       : default 'data/raw'.
%   opts.out_mat       : default 'data/ea_dataset.mat' (variant name is
%                        appended automatically for non-baseline variants).
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
opts = set_default(opts, 'shock_variant', 'mps_gc_1y');
opts = set_default(opts, 'y0m0', [2000 1]);
opts = set_default(opts, 'y1m1', [2019 12]);
opts = set_default(opts, 'raw_dir', fullfile('data', 'raw'));
opts = set_default(opts, 'shocks_mat', fullfile('data', 'derived', 'shocks_monthly.mat'));
opts = set_default(opts, 'out_mat', fullfile('data', 'ea_dataset.mat'));

% --- load shocks ----------------------------------------------------------
assert(exist(opts.shocks_mat, 'file') > 0, ...
       'assemble_dataset: run build_shock_series first (%s missing).', opts.shocks_mat);
S = load(opts.shocks_mat);
assert(isfield(S, opts.shock_variant), ...
       'assemble_dataset: unknown shock variant %s', opts.shock_variant);

% --- load outcomes ---------------------------------------------------------
[ym_ip,  v_ip]  = read_sdmx_csv(fullfile(opts.raw_dir, 'ip_ea.csv'));
[ym_pi,  v_pi]  = read_sdmx_csv(fullfile(opts.raw_dir, 'hicp_ea.csv'));
[ym_r,   v_r]   = read_sdmx_csv(fullfile(opts.raw_dir, 'rate1y_ea.csv'));
[ym_sx,  v_sx]  = read_sdmx_csv(fullfile(opts.raw_dir, 'stoxx50_ea.csv'));

% --- align -----------------------------------------------------------------
ym0 = 12 * opts.y0m0(1) + opts.y0m0(2);
ym1 = 12 * opts.y1m1(1) + opts.y1m1(2);
ym  = (ym0:ym1)';
T   = numel(ym);

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

% first-stage relevance: Delta i1y on mps (both in pp), effective 2001m1+
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
save(out, '-struct', 'ds');
fprintf('assemble_dataset: saved %s\n', out);
end

% -------------------------------------------------------------------------
function s = set_default(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

function v = pick(ym_src, val_src, ym_want)
% Align a source series onto the wanted month grid (NaN where absent).
v = nan(numel(ym_want), 1);
[tf, loc] = ismember(ym_want, ym_src);
v(tf) = val_src(loc(tf));
end
