function prot = ea_apply_protocol(sets, nl, varnames, h_early, held_blocks, out_csv, opts)
% EA_APPLY_PROTOCOL  The pre-registered tau reading protocol against a null.
%
% PURPOSE
% -------
% Applies the reading protocol of DESIGN_V1_LEGACY.md Sec. 5.1 / DESIGN.md Sec. 7
% to one or two scale estimators, using a stored null calibration
% (run_null_calibration).  This is the local apply_protocol of
% RUN_EMPIRICAL made callable from anywhere, with three additions:
%   * a DESIGN-KEY CHECK (ea_design_key): a null simulated under a
%     different design is refused unless opts.force is set;
%   * per cell, tau_bar relative to its null q95 (ratio_q95), the null
%     percentile of tau_bar among the R null draws, the naive one-sided
%     p-value p_cell = 1 - null_pct (which can be exactly 0), and the
%     MONTE CARLO p-value p_mc = (r + 1) / (R + 1) with r the number of
%     null draws at or above tau_bar (Davison & Hinkley 1997), whose
%     smallest attainable value 1/(R + 1) makes the resolution of the
%     null explicit;
%   * per cell, a FAMILY-WISE p-value p_fwer from the max statistic: in
%     every null draw take the largest tau_bar / q95 over the free cells;
%     p_fwer(i, g) = (1 + #{draws whose max >= ratio(i, g)}) / (R + 1).
%     This is the single-step max-T adjustment of Westfall & Young, using
%     the joint null distribution the calibration already provides, so it
%     controls the family-wise error rate without Holm's independence
%     bound and without needing p-values below 0.05 / m;
%   * per estimator, the escape COUNTS: cells above q95, cells rejected by
%     Holm's step-down at 5% over the free cells (on p_mc), cells with
%     p_fwer <= 5%, and equations whose max_g tau_bar exceeds the null
%     max-statistic q95.
% With R null draws Holm can only ever reject when 1/(R + 1) <= 0.05 / m,
% i.e. R >= 20 m - 1 (R >= 319 for the 16 free cells of a K = 4 system);
% a Holm count of 0 at R = 200 is therefore uninformative and the max-
% statistic p_fwer is the family-wise test to read.
% Blocks held at tau = 1 carry no signal and are excluded from every
% comparison and every count.
%
% INPUTS
% ------
% sets        : n x 2 cell, {name, est}: name in {'block', 'pooled'}, est a
%               struct with .tau_mean and .p_tau_gt1 (K x G x H)
% nl          : null struct, either the current layout (.block / .pooled
%               substructs) or the pre-revision flat layout (fields at the
%               top level, treated as 'block')
% varnames    : cellstr, K names (equations and blocks)
% h_early     : [h1 h2] window of the tau_bar statistic
% held_blocks : indices of blocks held at tau = 1 ([] = none)
% out_csv     : path of the long-format csv; '' = do not write
% opts        : optional struct
%     .design_key  key of the real-data run (ea_design_key); '' = no check
%     .force       false; true turns a key mismatch into a warning
%     .alpha       0.05 (Holm level)
%     .quiet       false; true suppresses the printed lines
%
% OUTPUTS
% -------
% prot : struct
%     .key_status   'match' | 'mismatch-forced' | 'legacy' | 'unchecked'
%     .design_key_run, .design_key_null
%     .varnames, .h_early, .held_blocks, .free (1 x K), .n_free_cells
%     .names        the estimators actually processed
%     .<name>       per estimator: .tau_bar, .pgt1_bar, .ratio_q95,
%                   .null_pct, .p_cell, .p_mc, .p_fwer, .escape_q95,
%                   .holm_reject, .fwer_reject (K x G), .eq_max, .eq_argmax,
%                   .eq_flag (K x 1), the counts .n_escape_q95,
%                   .n_escape_holm, .n_escape_fwer, .n_eq_flag,
%                   .n_free_cells, and .p_fwer_min, .null_R, .null_method
%
% NOTES
% -----
% * csv columns: those of the v1 tau_protocol csv, followed by
%   ratio_q95,null_pct,p_cell,holm_reject,p_mc,p_fwer,fwer_reject.
%   null_pct is the fraction of null draws STRICTLY below tau_bar; it,
%   p_mc and p_fwer need nl.<name>.tau_bar_draws (R x K x G) and are NaN
%   when the null does not carry the draws (Holm and FWER counts then NaN).
% * Holm: the free cells (all equations x free blocks) are ranked by
%   p_mc; cell j is rejected if every cell with a smaller or equal
%   p-value satisfies p <= alpha / (m - rank + 1).
% * The printed lines are those of RUN_EMPIRICAL step E, plus one count
%   line per estimator.

if nargin < 7 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'design_key'), opts.design_key = ''; end
if ~isfield(opts, 'force'), opts.force = false; end
if ~isfield(opts, 'alpha'), opts.alpha = 0.05; end
if ~isfield(opts, 'quiet'), opts.quiet = false; end
if nargin < 6, out_csv = ''; end
if ~iscell(varnames), varnames = cellstr(varnames); end
vn = reshape(varnames, 1, []);
K = numel(vn);
say = @(varargin) ifprint(~opts.quiet, varargin{:});

% --- design-key check ---------------------------------------------------------
null_key = '';
if isfield(nl, 'design_key') && ~isempty(nl.design_key), null_key = nl.design_key; end
if isempty(null_key)
    key_status = 'legacy';
    say('  WARNING: legacy null without design key (no design check possible)\n');
elseif isempty(opts.design_key)
    key_status = 'unchecked';
elseif strcmp(opts.design_key, null_key)
    key_status = 'match';
else
    msg = sprintf(['ea_apply_protocol: design key of the null differs from the run''s.\n' ...
                   '    run : %s\n    null: %s'], opts.design_key, null_key);
    if opts.force
        warning('ea_apply_protocol:design_key', '%s\n  (opts.force = true: proceeding anyway)', msg);
        key_status = 'mismatch-forced';
    else
        error('ea_apply_protocol:design_key', '%s\n  Pass opts.force = true to override.', msg);
    end
end

% --- layout of the null -------------------------------------------------------
n_rep = getf(nl, 'n_rep', NaN);  p_null = getf(nl, 'p', NaN);

free = true(1, K);  free(held_blocks) = false;
n_free_cells = K * sum(free);

prot = struct('key_status', key_status, 'design_key_run', opts.design_key, ...
              'design_key_null', null_key, 'varnames', {vn}, 'h_early', h_early, ...
              'held_blocks', held_blocks, 'free', free, 'n_free_cells', n_free_cells, ...
              'names', {{}});

fid = -1;
if ~isempty(out_csv)
    fid = fopen(out_csv, 'w');
    assert(fid > 0, 'ea_apply_protocol: cannot write %s', out_csv);
    fprintf(fid, ['scale_estimator,equation,block,block_held_at_1,tau_bar,null_mean,null_q90,null_q95,' ...
                  'null_q99,escape_q95,pgt1_bar,null_pgt1_q95,is_argmax,null_argmax_freq,' ...
                  'eq_maxstat,eq_maxstat_q95,eq_flag,ratio_q95,null_pct,p_cell,holm_reject,p_mc,p_fwer,fwer_reject\n']);
end

for e = 1:size(sets, 1)
    nm = sets{e, 1};  est = sets{e, 2};
    if isfield(nl, nm)
        S = nl.(nm);
    elseif strcmp(nm, 'block') && isfield(nl, 'q95')
        S = nl;                                  % pre-revision flat layout
    else
        say('  [%s] skipped: the null carries no %s estimator\n', nm, nm);
        continue
    end
    tb = mean(est.tau_mean(:, :, h_early(1):h_early(2)), 3);        % K x G
    pb = mean(est.p_tau_gt1(:, :, h_early(1):h_early(2)), 3);
    meth = 'resample';
    if isfield(nl, 'sim_method') && ~isempty(nl.sim_method), meth = nl.sim_method; end
    say('  [%s] reading protocol (null R = %d, p = %d, innovations ''%s'', window h = %d..%d):\n', ...
        nm, n_rep, p_null, meth, h_early(1), h_early(2));

    % per-cell statistics
    ratio = tb ./ S.q95;
    npct = nan(K, K);  pmc = nan(K, K);  pfw = nan(K, K);  R = NaN;
    if isfield(S, 'tau_bar_draws') && ~isempty(S.tau_bar_draws)
        D = S.tau_bar_draws;                     % R x K x G
        R = size(D, 1);
        % null distribution of the family-wise MAX statistic: in each null
        % draw, the largest ratio to the cell's own q95 over the free cells
        Dr = D;
        for i = 1:K, for g = 1:K, Dr(:, i, g) = D(:, i, g) / S.q95(i, g); end, end
        Dr(:, :, ~free) = -Inf;
        mx = max(reshape(Dr, R, []), [], 2);
        for i = 1:K
            for g = 1:K
                npct(i, g) = mean(D(:, i, g) < tb(i, g));
                pmc(i, g)  = (sum(D(:, i, g) >= tb(i, g)) + 1) / (R + 1);
                pfw(i, g)  = (sum(mx >= ratio(i, g)) + 1) / (R + 1);
            end
        end
    end
    pcell = 1 - npct;                            % naive r/R (kept: can be exactly 0)
    esc = false(K, K);  eq_flag = false(K, 1);  eq_max = zeros(K, 1);  eq_am = zeros(K, 1);
    for i = 1:K
        tbi = tb(i, :);  tbi(~free) = -inf;
        [mx, am] = max(tbi);
        eq_max(i) = mx;  eq_am(i) = am;
        eq_flag(i) = mx > S.maxstat_q95(i);
        esc(i, :) = free & tb(i, :) > S.q95(i, :);
    end
    holm = holm_stepdown(pmc, free, opts.alpha);          % Monte Carlo p-values
    fwer = free & (pfw <= opts.alpha);                    % max-statistic test

    for i = 1:K
        mx = eq_max(i);  am = eq_am(i);
        say('    eq %-6s max tau_bar %.3f (null q95 %.3f) %s | argmax block %s (null freq %.2f)', ...
            vn{i}, mx, S.maxstat_q95(i), tern(eq_flag(i), 'ESCAPE', 'quiet '), ...
            vn{am}, S.argmax_freq(i, am));
        ei = find(esc(i, :));
        if isempty(ei), say(' | no cell above q95\n');
        else, say(' | cells above q95: %s\n', strjoin(vn(ei), ', ')); end
        if fid > 0
            for g = 1:K
                pq = NaN;
                if isfield(S, 'pgt1_q95'), pq = S.pgt1_q95(i, g); end
                fprintf(fid, '%s,%s,%s,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%d,%.4f,%.4f,%d,%.4f,%.5f,%.5f,%d,%.5f,%.4f,%.4f,%d,%.4f,%.4f,%d\n', ...
                        nm, vn{i}, vn{g}, ~free(g), tb(i, g), S.tau_bar_mean(i, g), ...
                        S.q90(i, g), S.q95(i, g), S.q99(i, g), esc(i, g), ...
                        pb(i, g), pq, g == am, S.argmax_freq(i, g), mx, S.maxstat_q95(i), eq_flag(i), ...
                        ratio(i, g), npct(i, g), pcell(i, g), holm(i, g), pmc(i, g), pfw(i, g), fwer(i, g));
            end
        end
    end
    n_holm = sum(holm(:));  n_fwer = sum(fwer(:));
    pf = pmc(:, free);  pfree = pfw(:, free);
    p_fwer_min = NaN;
    if all(isnan(pf(:))), n_holm = NaN;  n_fwer = NaN; else, p_fwer_min = min(pfree(:)); end
    say('    counts: %d free cells | %d above q95 | %d Holm-rejected (%.0f%%, p_mc; min attainable %.4f) | %d FWER max-stat rejected (min p_fwer %.3f) | %d of %d equations flagged\n', ...
        n_free_cells, sum(esc(:)), n_holm, 100 * opts.alpha, 1 / (R + 1), n_fwer, p_fwer_min, sum(eq_flag), K);
    if ~isnan(R) && 1 / (R + 1) > opts.alpha / n_free_cells
        say('    note: with R = %d null draws Holm cannot reject at %.0f%% over %d cells (needs R >= %d); read p_fwer.\n', ...
            R, 100 * opts.alpha, n_free_cells, ceil(n_free_cells / opts.alpha) - 1);
    end

    E = struct('tau_bar', tb, 'pgt1_bar', pb, 'null_q95', S.q95, 'ratio_q95', ratio, ...
               'null_pct', npct, 'p_cell', pcell, 'p_mc', pmc, 'p_fwer', pfw, ...
               'escape_q95', esc, 'holm_reject', holm, 'fwer_reject', fwer, ...
               'eq_max', eq_max, 'eq_argmax', eq_am, 'eq_flag', eq_flag, ...
               'n_free_cells', n_free_cells, 'n_escape_q95', sum(esc(:)), ...
               'n_escape_holm', n_holm, 'n_escape_fwer', n_fwer, 'p_fwer_min', p_fwer_min, ...
               'n_eq_flag', sum(eq_flag), 'null_R', R, 'null_method', meth);
    prot.(nm) = E;
    prot.names{end + 1} = nm;
end
if fid > 0
    fclose(fid);
    say('  wrote %s\n', out_csv);
end
end

% -------------------------------------------------------------------------
function rej = holm_stepdown(pcell, free, alpha)
% Holm step-down at level alpha over the free cells (all rows x free
% columns).  NaN p-values are never rejected and do not count towards m.
K = size(pcell, 1);
rej = false(size(pcell));
cand = repmat(free, K, 1) & ~isnan(pcell);
idx = find(cand);
m = numel(idx);
if m == 0, return; end
[ps, ord] = sort(pcell(idx));
thr = alpha ./ (m - (1:m) + 1);
ok = ps(:)' <= thr;
n_rej = find(~ok, 1) - 1;
if isempty(n_rej), n_rej = m; end
rej(idx(ord(1:n_rej))) = true;
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function ifprint(on, varargin)
if on, fprintf(varargin{:}); end
end
