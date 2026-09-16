function res = run_block_ablation(opts)
% PURPOSE
% -------
% Out-of-sample BLOCK-ABLATION validation of the escapes found by the
% block-adaptive BLP (empirical/docs/DESIGN.md, Section 8).  Without
% ground truth, an escape of cell (i, g) -- equation i letting the lag
% block of variable g leave the BVAR prior -- is credible if releasing
% that block buys predictive accuracy.  So, in a pseudo out-of-sample
% expanding-window design, the horizon-h direct forecast of y_{i,T0+h} is
% formed under
%   'free'    : the block-adaptive posterior of equation i, and
%   'ablated' : the same posterior with cell (i, g) FORCED back to tau = 1
%               via cfg.blocks.fixed_tau_mask,
% and the two forecast-error series are compared.  Prediction: forcing an
% escaping cell back worsens the MSFE; forcing a quiet CONTROL cell does
% not.
%
% MODEL / EQUATIONS
% -----------------
% At each origin T0 the whole pipeline is re-estimated on Y(1:T0, :):
%   bvar0 = estimate_bvar_niw,   blpf0 = estimate_blp_fmar (gives lambda_h),
%   blp   = estimate_blp_blockadaptive(..., blpf0.lambda) for equation i
%           only (cfg.blp.equations = i), free and ablated.
% The LPs are run on the BVAR-detrended data x_t = y_t - trend_t
% (utils/var_deterministic_trend.m), so the direct forecast is
%     yhat_{i,T0+h} = z_{T0}' beta_mean(:, i, h) + trend_{i,T0+h},
%     z_{T0} = [1; x_{T0}; x_{T0-1}; ...; x_{T0-p+1}]   (build_lp_regressors
%                                                        layout, last row),
% where trend is continued beyond T0 with the fitted VAR recursion
%     trend(t) = c + sum_j A_j trend(t-j),   t = T0+1..T0+h.
% References: the FMAR closed form (blpf0.beta_mean, same construction)
% and the BVAR iterated forecast of the levels from bvar0.c, bvar0.A.
%
% Per cell and per h in opts.h_eval (plus 'pooled', the per-origin average
% of the loss differential over h_eval), with e_r = ablated (restricted)
% and e_u = free (unrestricted) errors over n origins:
%   ratio   = MSFE_ablated / MSFE_free;
%   DM      = mean(d) / sqrt(HAC(d)/n) * HLN,  d_t = e_r^2 - e_u^2,
%             Newey-West Bartlett HAC with truncation h (max h_eval when
%             pooled) and the Harvey-Leybourne-Newbold small-sample
%             factor sqrt((n + 1 - 2h + h(h-1)/n)/n);
%   CW      = mean(f) / sqrt(HAC(f)/n),  Clark-West (2007)
%             f_t = e_r^2 - (e_u^2 - (yhat_r - yhat_u)^2),
%             the right correction for NESTED predictors (ablated is
%             nested in free).
% Positive DM / CW = the ablation forecasts worse.
%
% INPUTS (opts, all optional)
% ---------------------------
%   .dataset        path to a .mat with Y (T x K), varnames, ym, isrw
%                   (default ea_paths().dataset, the ois4 headline system),
%                   or that struct itself
%   .cells          (n x 2) cell array {equation_name, block_name} to
%                   ablate; names from varnames (block g = variable g
%                   under the per-variable block scheme), or indices
%   .control_cells  same format, quiet cells that should NOT matter
%   .p, .H          lag order and LP horizon (default 12, 12)
%   .h_eval         horizons entering the loss (default [2 12])
%   .origin0        first forecast origin (default round(0.6*T))
%   .origin_step    step between origins (default 1)
%   .n_origins      cap on the number of origins ([] = all up to T-max h)
%   .fixed_tau      cfg.blocks.fixed_tau for every run (default [])
%   .cfg_over       struct merged (recursively) onto the configuration
%   .gibbs_n_burn/.gibbs_n_keep   chain lengths (default 300 / 700)
%   .seed           base seed (default 7101); free and ablated runs of one
%                   (origin, equation) share the same seed
%   .out_stem       results are results/empirical/ablation_<stem>.mat/.csv
%   .out_dir        output folder (default ea_paths().results)
%   .quick          true = tiny settings for tests (p = 2, H = 4,
%                   h_eval = [2 4], 3 origins, chains 40 + 80, 40 NIW
%                   draws); explicit opts still win
%   .allow_synthetic  must be true to run on a dataset stamped
%                   ds.synthetic = true (interface tests only)
%   .verbose        print the table (default true)
%
% OUTPUTS
% -------
% res : struct
%   .table   struct array, one row per (cell, h) and one pooled row per
%            cell (h = 0): cell, role ('escape'/'control'), h, n,
%            msfe_free, msfe_ablated, ratio, dm_stat, cw_stat
%   .cells   struct array of the resolved cells (i, g, names, role)
%   .origins forecast origins; .h_eval; .varnames
%   .err     .free/.ablated (nO x H x nCells) forecast errors,
%            .fmar/.bvar (nO x H x K) reference errors
%   .yhat    the corresponding forecasts; .yact (nO x H x K) actuals
%   .tau_free/.p_gt1_free (nO x H x nCells) posterior tau summaries of
%            the ablated cell in the FREE run at each origin
%   .ref     MSFE of free / FMAR / BVAR per involved equation and h
%   .cfg, .opts, .files
%
% NOTES
% -----
% The trend continuation is a local helper (var_continue) so that
% utils/var_deterministic_trend.m stays untouched.  Numbers computed on
% the synthetic fixture are not empirical findings.

if nargin < 1, opts = struct(); end
if exist('ea_paths', 'file') ~= 2
    addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
P = ea_paths();

% --- options: user > quick > defaults ---------------------------------------
if isfield(opts, 'quick') && ~isempty(opts.quick) && opts.quick
    opts = setdef(opts, 'p', 2);              opts = setdef(opts, 'H', 4);
    opts = setdef(opts, 'h_eval', [2 4]);     opts = setdef(opts, 'n_origins', 3);
    opts = setdef(opts, 'gibbs_n_burn', 40);  opts = setdef(opts, 'gibbs_n_keep', 80);
    opts = setdef(opts, 'n_niw_draws', 40);
end
opts = setdef(opts, 'quick', false);
opts = setdef(opts, 'dataset', P.dataset);
opts = setdef(opts, 'cells', {});
opts = setdef(opts, 'control_cells', {});
opts = setdef(opts, 'p', 12);                 opts = setdef(opts, 'H', 12);
opts = setdef(opts, 'h_eval', [2 12]);
opts = setdef(opts, 'origin0', []);           opts = setdef(opts, 'origin_step', 1);
opts = setdef(opts, 'n_origins', []);
opts = setdef(opts, 'fixed_tau', []);         opts = setdef(opts, 'cfg_over', struct());
opts = setdef(opts, 'gibbs_n_burn', 300);     opts = setdef(opts, 'gibbs_n_keep', 700);
opts = setdef(opts, 'n_niw_draws', []);
opts = setdef(opts, 'seed', 7101);
opts = setdef(opts, 'out_stem', 'default');   opts = setdef(opts, 'out_dir', P.results);
opts = setdef(opts, 'allow_synthetic', false);
opts = setdef(opts, 'verbose', true);
assert(all(opts.h_eval >= 1 & opts.h_eval <= opts.H), ...
    'run_block_ablation: h_eval must lie in 1..H.');

% --- data ---------------------------------------------------------------------
if isstruct(opts.dataset)
    ds = opts.dataset;  ds_path = '<struct passed in>';
else
    ds_path = opts.dataset;
    assert(exist(ds_path, 'file') == 2, 'run_block_ablation: dataset not found: %s', ds_path);
    ds = load(ds_path);
end
if isfield(ds, 'synthetic') && ds.synthetic && ~opts.allow_synthetic
    error(['run_block_ablation: this dataset is a SYNTHETIC FIXTURE ' ...
           '(ds.synthetic = true); set opts.allow_synthetic = true only for ' ...
           'interface tests.  Nothing computed from it is an empirical result.']);
end
Y = ds.Y;  [T, K] = size(Y);
varnames = ds.varnames;
p = opts.p;  H = opts.H;  G = K;             % per-variable blocks: G = K
m = 1 + K * p;

% --- configuration -----------------------------------------------------------
cfg = default_config();
cfg.mode = 'fmar';
cfg.fmar.h1_mode   = 'lp';
cfg.fmar.psi_floor = true;
cfg.fmar.isrw      = ds.isrw;
cfg.blocks.fixed_tau = opts.fixed_tau;
cfg.K = K;  cfg.p = p;  cfg.H = H;  cfg.T = T;
cfg.gibbs.n_burn = opts.gibbs_n_burn;
cfg.gibbs.n_keep = opts.gibbs_n_keep;
if ~isempty(opts.n_niw_draws), cfg.fmar.n_niw_draws = opts.n_niw_draws; end
cfg = merge_struct(cfg, opts.cfg_over);
assert(strcmp(cfg.blocks.scheme, 'per_variable'), ...
    'run_block_ablation: written for the per-variable block scheme.');

% --- cells ---------------------------------------------------------------------
cells = [resolve_cells(opts.cells, varnames, 'escape'), ...
         resolve_cells(opts.control_cells, varnames, 'control')];
nC = numel(cells);
assert(nC >= 1, 'run_block_ablation: give at least one cell in opts.cells / opts.control_cells.');
eqs = unique([cells.i]);  nE = numel(eqs);
if ~isempty(opts.fixed_tau)
    held = ismember([cells.g], opts.fixed_tau(:)');
    if any(held)
        warning('run_block_ablation:heldCell', ...
            'cell(s) %s lie in a block already held by opts.fixed_tau; ablated == free there.', ...
            strjoin({cells(held).label}, ', '));
    end
end

% --- origins ---------------------------------------------------------------------
hmax = max(opts.h_eval);
if isempty(opts.origin0), opts.origin0 = round(0.6 * T); end
origins = opts.origin0 : opts.origin_step : (T - hmax);
if ~isempty(opts.n_origins), origins = origins(1:min(end, opts.n_origins)); end
nO = numel(origins);
assert(nO >= 1, 'run_block_ablation: no admissible origin (origin0 + max(h_eval) must be <= T).');
assert(origins(1) - p - H > m + 5, ...
    'run_block_ablation: origin0 = %d leaves too few observations for p = %d, H = %d.', ...
    origins(1), p, H);

% --- containers ----------------------------------------------------------------
yact      = nan(nO, H, K);
yhat_fmar = nan(nO, H, K);   yhat_bvar = nan(nO, H, K);
yhat_free = nan(nO, H, nC);  yhat_abl  = nan(nO, H, nC);
tau_free  = nan(nO, H, nC);  pgt1_free = nan(nO, H, nC);
yhat_free_eq = nan(nO, H, nE);

if opts.verbose
    fprintf('run_block_ablation: T = %d, K = %d, p = %d, H = %d, %d origins (%d..%d), %d cells, %d equation(s)\n', ...
        T, K, p, H, nO, origins(1), origins(end), nC, nE);
end
t_start = tic;
for o = 1:nO
    T0 = origins(o);
    Y0 = Y(1:T0, :);
    cfg0 = cfg;  cfg0.T = T0;

    % (1) BVAR and FMAR baseline on the window
    rng(opts.seed + o, 'twister');
    bvar0 = estimate_bvar_niw(Y0, cfg0);
    blpf0 = estimate_blp_fmar(Y0, cfg0, bvar0);

    % (2) regressor at the origin and deterministic continuation
    dt = var_deterministic_trend(Y0, bvar0.B, p);
    Z0 = build_lp_regressors(dt.x, p);
    z  = Z0(end, :)';                                 % t = T0
    trend_f = var_continue(dt.trend, bvar0.c, bvar0.A, H);   % (H x K), T0+1..T0+H
    ybvar   = var_continue(Y0,       bvar0.c, bvar0.A, H);   % iterated BVAR forecast

    for h = 1:H
        if T0 + h <= T, yact(o, h, :) = Y(T0 + h, :); end
        yhat_bvar(o, h, :) = ybvar(h, :);
        for i = 1:K
            yhat_fmar(o, h, i) = z' * blpf0.beta_mean(:, i, h) + trend_f(h, i);
        end
    end

    % (3) block-adaptive, one equation at a time: free, then each ablation
    for e = 1:nE
        i = eqs(e);
        cfg_e = cfg0;  cfg_e.blp.equations = i;
        seed_e = opts.seed + 1000 * o + i;
        cfg_e.blocks.fixed_tau_mask = [];
        rng(seed_e, 'twister');
        bf = estimate_blp_blockadaptive(Y0, cfg_e, bvar0, blpf0.lambda);
        for h = 1:H
            yhat_free_eq(o, h, e) = z' * bf.beta_mean(:, i, h) + trend_f(h, i);
        end
        for c = find([cells.i] == i)
            g = cells(c).g;
            mask = false(K, G);  mask(i, g) = true;
            cfg_e.blocks.fixed_tau_mask = mask;
            rng(seed_e, 'twister');                   % same stream as the free run
            ba = estimate_blp_blockadaptive(Y0, cfg_e, bvar0, blpf0.lambda);
            for h = 1:H
                yhat_free(o, h, c) = yhat_free_eq(o, h, e);
                yhat_abl(o, h, c)  = z' * ba.beta_mean(:, i, h) + trend_f(h, i);
                tau_free(o, h, c)  = bf.tau_mean(i, g, h);
                pgt1_free(o, h, c) = bf.p_tau_gt1(i, g, h);
            end
        end
    end
    if opts.verbose
        fprintf('  origin %3d/%3d (T0 = %d) done, %.0f s elapsed\n', o, nO, T0, toc(t_start));
    end
end

% --- errors --------------------------------------------------------------------
e_fmar = yact - yhat_fmar;   e_bvar = yact - yhat_bvar;
e_free = nan(nO, H, nC);     e_abl  = nan(nO, H, nC);
for c = 1:nC
    e_free(:, :, c) = yact(:, :, cells(c).i) - yhat_free(:, :, c);
    e_abl(:, :, c)  = yact(:, :, cells(c).i) - yhat_abl(:, :, c);
end

% --- statistics -----------------------------------------------------------------
rows = struct('cell', {}, 'role', {}, 'h', {}, 'n', {}, 'msfe_free', {}, ...
              'msfe_ablated', {}, 'ratio', {}, 'dm_stat', {}, 'cw_stat', {});
for c = 1:nC
    D = nan(nO, numel(opts.h_eval));  Fcw = nan(nO, numel(opts.h_eval));
    for k = 1:numel(opts.h_eval)
        h = opts.h_eval(k);
        [r, d, f] = cell_stats(e_abl(:, h, c), e_free(:, h, c), ...
                               yhat_abl(:, h, c), yhat_free(:, h, c), h);
        D(:, k) = d;  Fcw(:, k) = f;
        rows(end + 1) = struct('cell', cells(c).label, 'role', cells(c).role, 'h', h, ...
            'n', r.n, 'msfe_free', r.msfe_u, 'msfe_ablated', r.msfe_r, ...
            'ratio', r.ratio, 'dm_stat', r.dm, 'cw_stat', r.cw); %#ok<AGROW>
    end
    % pooled over h_eval: per-origin average of the loss differentials
    ok = all(isfinite(D), 2);
    dp = mean(D(ok, :), 2);  fp = mean(Fcw(ok, :), 2);
    e2r = mean(reshape(e_abl(ok, opts.h_eval, c), [], 1).^2);
    e2u = mean(reshape(e_free(ok, opts.h_eval, c), [], 1).^2);
    rows(end + 1) = struct('cell', cells(c).label, 'role', cells(c).role, 'h', 0, ...
        'n', sum(ok), 'msfe_free', e2u, 'msfe_ablated', e2r, 'ratio', e2r / e2u, ...
        'dm_stat', hac_tstat(dp, hmax) * hln_factor(sum(ok), hmax), ...
        'cw_stat', hac_tstat(fp, hmax)); %#ok<AGROW>
end

% reference MSFEs per involved equation
ref = struct('eq', {}, 'h', {}, 'msfe_free', {}, 'msfe_fmar', {}, 'msfe_bvar', {});
for e = 1:nE
    i = eqs(e);
    for h = opts.h_eval
        ef = yact(:, h, i) - yhat_free_eq(:, h, e);
        ok = isfinite(ef);
        ref(end + 1) = struct('eq', varnames{i}, 'h', h, ...
            'msfe_free', mean(ef(ok).^2), 'msfe_fmar', mean(e_fmar(ok, h, i).^2), ...
            'msfe_bvar', mean(e_bvar(ok, h, i).^2)); %#ok<AGROW>
    end
end

% --- assemble, save, print -------------------------------------------------------
res.table    = rows;
res.ref      = ref;
res.cells    = cells;
res.origins  = origins;
res.h_eval   = opts.h_eval;
res.varnames = varnames;
res.err  = struct('free', e_free, 'ablated', e_abl, 'fmar', e_fmar, 'bvar', e_bvar);
res.yhat = struct('free', yhat_free, 'ablated', yhat_abl, 'fmar', yhat_fmar, 'bvar', yhat_bvar);
res.yact = yact;
res.tau_free   = tau_free;
res.p_gt1_free = pgt1_free;
res.cfg  = cfg;
res.opts = opts;
res.dataset_path = ds_path;
res.synthetic = isfield(ds, 'synthetic') && logical(ds.synthetic);
res.elapsed_s = toc(t_start);

if ~exist(opts.out_dir, 'dir'), mkdir(opts.out_dir); end
stem = ['ablation_' opts.out_stem];
res.files.mat = fullfile(opts.out_dir, [stem '.mat']);
res.files.csv = fullfile(opts.out_dir, [stem '.csv']);
save(res.files.mat, 'res');
write_csv(res.files.csv, rows);

if opts.verbose
    fprintf('\nBlock ablation (%s): MSFE ratio ablated/free, DM (HLN) and Clark-West statistics\n', opts.out_stem);
    fprintf('  %-18s %-8s %-6s %4s %10s %10s %7s %8s %8s\n', ...
        'cell', 'role', 'h', 'n', 'msfe_free', 'msfe_abl', 'ratio', 'DM', 'CW');
    for r = 1:numel(rows)
        if rows(r).h == 0, hs = 'pooled'; else, hs = sprintf('%d', rows(r).h); end
        fprintf('  %-18s %-8s %-6s %4d %10.4g %10.4g %7.3f %8.2f %8.2f\n', ...
            rows(r).cell, rows(r).role, hs, rows(r).n, rows(r).msfe_free, ...
            rows(r).msfe_ablated, rows(r).ratio, rows(r).dm_stat, rows(r).cw_stat);
    end
    fprintf('  reference MSFE (free / FMAR / BVAR):\n');
    for r = 1:numel(ref)
        fprintf('    %-8s h = %2d: %10.4g %10.4g %10.4g\n', ref(r).eq, ref(r).h, ...
            ref(r).msfe_free, ref(r).msfe_fmar, ref(r).msfe_bvar);
    end
    fprintf('  saved %s and .csv (%.0f s)\n', res.files.mat, res.elapsed_s);
end
end

% =====================================================================
function s = setdef(s, f, v)
if ~isfield(s, f) || (isempty(s.(f)) && ~iscell(s.(f))), s.(f) = v; end
end

function base = merge_struct(base, over)
% Recursive merge: fields of `over` replace / extend those of `base`.
f = fieldnames(over);
for k = 1:numel(f)
    if isfield(base, f{k}) && isstruct(base.(f{k})) && isstruct(over.(f{k}))
        base.(f{k}) = merge_struct(base.(f{k}), over.(f{k}));
    else
        base.(f{k}) = over.(f{k});
    end
end
end

function cells = resolve_cells(list, varnames, role)
% {equation, block} pairs (names or indices) -> struct array with i, g.
cells = struct('i', {}, 'g', {}, 'eq_name', {}, 'block_name', {}, 'role', {}, 'label', {});
if isempty(list), return; end
assert(iscell(list) && size(list, 2) == 2, ...
    'run_block_ablation: %s cells must be an (n x 2) cell array {equation, block}.', role);
for r = 1:size(list, 1)
    i = name_to_index(list{r, 1}, varnames);
    g = name_to_index(list{r, 2}, varnames);
    cells(end + 1) = struct('i', i, 'g', g, 'eq_name', varnames{i}, ...
        'block_name', varnames{g}, 'role', role, ...
        'label', sprintf('%s_on_%s', varnames{i}, varnames{g})); %#ok<AGROW>
end
end

function k = name_to_index(v, varnames)
if ischar(v) || isstring(v)
    k = find(strcmp(varnames, char(v)));
    assert(isscalar(k), 'run_block_ablation: unknown variable ''%s''.', char(v));
else
    k = v;
    assert(isscalar(k) && k >= 1 && k <= numel(varnames) && k == round(k), ...
        'run_block_ablation: variable index out of range.');
end
end

function Xf = var_continue(X, c, A, H)
% Continue the VAR(p) recursion x(t) = c + sum_j A_j x(t-j) for H steps
% beyond the last row of X (no shocks).  Used for the deterministic
% component (X = trend) and for the iterated BVAR forecast (X = levels).
[T0, K] = size(X);  p = size(A, 3);
Xa = [X; zeros(H, K)];
for t = T0 + 1:T0 + H
    v = c(:);
    for j = 1:p
        v = v + A(:, :, j) * Xa(t - j, :)';
    end
    Xa(t, :) = v';
end
Xf = Xa(T0 + 1:end, :);
end

function [r, d, f] = cell_stats(e_r, e_u, yhat_r, yhat_u, h)
% Restricted (ablated) vs unrestricted (free) forecast-error series.
ok = isfinite(e_r) & isfinite(e_u);
d = nan(size(e_r));  f = nan(size(e_r));
d(ok) = e_r(ok).^2 - e_u(ok).^2;
f(ok) = e_r(ok).^2 - (e_u(ok).^2 - (yhat_r(ok) - yhat_u(ok)).^2);
r.n = sum(ok);
r.msfe_r = mean(e_r(ok).^2);
r.msfe_u = mean(e_u(ok).^2);
r.ratio  = r.msfe_r / r.msfe_u;
r.dm = hac_tstat(d(ok), h) * hln_factor(r.n, h);
r.cw = hac_tstat(f(ok), h);
end

function t = hac_tstat(v, L)
% t-statistic of mean(v) with a Newey-West (Bartlett) long-run variance,
% truncation L: S = g0 + 2 sum_{l=1}^{L} (1 - l/(L+1)) g_l.
v = v(isfinite(v));  n = numel(v);
if n < 2, t = NaN; return; end
vc = v - mean(v);
S = sum(vc.^2) / n;
for l = 1:min(L, n - 1)
    S = S + 2 * (1 - l / (L + 1)) * sum(vc(1 + l:end) .* vc(1:end - l)) / n;
end
if S <= 0, t = NaN; return; end
t = mean(v) / sqrt(S / n);
end

function k = hln_factor(n, h)
% Harvey, Leybourne & Newbold (1997) small-sample correction of the DM
% statistic for h-step forecasts.
k = sqrt(max((n + 1 - 2 * h + h * (h - 1) / n) / n, 0));
end

function write_csv(fname, rows)
fid = fopen(fname, 'w');
assert(fid > 0, 'run_block_ablation: cannot write %s', fname);
fprintf(fid, 'cell,role,h,n,msfe_free,msfe_ablated,ratio,dm_stat,cw_stat\n');
for r = 1:numel(rows)
    if rows(r).h == 0, hs = 'pooled'; else, hs = sprintf('%d', rows(r).h); end
    fprintf(fid, '%s,%s,%s,%d,%.10g,%.10g,%.10g,%.10g,%.10g\n', ...
        rows(r).cell, rows(r).role, hs, rows(r).n, rows(r).msfe_free, ...
        rows(r).msfe_ablated, rows(r).ratio, rows(r).dm_stat, rows(r).cw_stat);
end
fclose(fid);
end
