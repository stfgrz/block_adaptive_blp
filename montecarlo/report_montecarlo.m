function report_montecarlo(src, fid)
% PURPOSE
% -------
% Print the headline tables of one Monte Carlo experiment: the two
% integrated-RMSE conventions, the bias-variance decomposition, interval
% performance under both band constructions, and the tau diagnostics.
% This is the single place the write-up's numbers are read off, so that
% a table in the thesis and a table in a log cannot drift apart.
%
% USAGE
% -----
%   report_montecarlo('results/mc_final_sparse.mat')   % a saved result
%   report_montecarlo(s)                               % a summary struct
%   report_montecarlo(src, fid)                        % write to a file
%
% WHAT IS PRINTED, AND HOW TO READ IT
% -----------------------------------
% 1. INTEGRATED RMSE.  Two columns per response: the LEGACY average over
%    h = 1..H and the PREFERRED average over h = 2..H.  Quote the second.
%    Under the FMAR convention (cfg.fmar.h1_mode = 'bvar') the global
%    baseline reports the Bayesian VAR at h = 1 while the adaptive
%    estimators report a local projection, so the legacy metric mixes
%    "VAR vs LP" into "global vs adaptive"; the preferred metric does
%    not.  Ratios below 1 mean the adaptive estimator beats the global
%    baseline.
% 2. BIAS-VARIANCE.  Mean |bias|, mean variance and RMSE over h = 2..H,
%    and separately over the early and late windows.  This is the table
%    that answers "does adaptation trade bias for variance, and where?".
%    MSE = bias^2 + variance exactly, by construction.
% 3. INTERVALS.  Coverage and average length of the PRIMARY bands (the
%    quasi-Bayesian Newey-West sandwich) and, for the sampled
%    estimators, of the posterior-quantile bands.  A gap between the two
%    is the price of the quasi-Bayesian construction
%    (docs/APPROXIMATIONS.md, item 3).
% 4. TAU DIAGNOSTICS.  Localisation of the misspecified block, or -- on
%    the correctly specified DGP -- the false-positive rates of the same
%    rules.  Read the two together: a detection rate means nothing
%    without the rate at which the same rule fires when nothing is
%    wrong.
%
% MONTE CARLO PRECISION.  The standard error of an RMSE is about
% rmse/sqrt(2R); the header prints it so no difference is read as real
% when it is not.  Differences BETWEEN estimators are more precise than
% that, because every estimator sees the same simulated samples.
%
% INPUTS
% ------
% src : path to a .mat containing `s` (or `mc`), or a summary struct.
% fid : optional file identifier (default 1, the terminal).

if nargin < 2 || isempty(fid), fid = 1; end
if ischar(src)
    L = load(src);
    % A summary saved before the h = 2..H metric and the bias-variance
    % fields existed is STALE: recompute it from the stored replication
    % output rather than reporting a subset of the tables (or failing on
    % a missing field).  This is what lets the earlier R = 500 results be
    % re-read under the preferred metric without re-running them.
    if isfield(L, 's') && isfield(L.s, 'irmse_h2')
        s = L.s;
        if isfield(L, 'mc'), s.mc_ref = L.mc; end
    elseif isfield(L, 'mc')
        s = summarize_montecarlo(L.mc);
        s.mc_ref = L.mc;
        fprintf(1, ['[report_montecarlo] the stored summary predates the ' ...
                    'current metrics; recomputed from mc.]\n']);
    else
        s = L.s;
    end
else
    s = src;
end

nE = numel(s.est_names);
K = s.K;  H = s.H;  R = s.R;
ig = find(strcmp(s.est_names, 'BLP-FMAR'), 1);
if isempty(ig), ig = find(strcmp(s.est_names, 'BLP-glob'), 1); end

fprintf(fid, '\n========================================================================\n');
fprintf(fid, ' Monte Carlo: DGP = %s', s.dgp_name);
if isfield(s, 'meta')
    m = s.meta;
    fprintf(fid, ' | mode = %s | R = %d | T = %d | p = %d | H = %d', ...
            m.mode, R, m.T, m.p, H);
    if isfield(m, 'preset'), fprintf(fid, ' | preset = %s', m.preset); end
    if isfield(m, 'fmar') && isfield(m.fmar, 'h1_mode')
        fprintf(fid, ' | h1 = %s', m.fmar.h1_mode);
    end
end
fprintf(fid, '\n');
if isfield(s, 'description'), fprintf(fid, ' %s\n', s.description); end
if isfield(s, 'meta') && isfield(s.meta, 'is_complete') && ~s.meta.is_complete
    fprintf(fid, ' *** INCOMPLETE: %d of %d replications ***\n', ...
            s.meta.R_done, s.meta.R);
end
fprintf(fid, ' Monte Carlo se of an RMSE ~ rmse/sqrt(2R) = %.1f%% of it.\n', ...
        100 / sqrt(2 * R));
fprintf(fid, '========================================================================\n');

% --- 1. integrated RMSE --------------------------------------------------
fprintf(fid, '\n[1] INTEGRATED RMSE   (legacy h=1..%d | PREFERRED h=2..%d)\n', H, H);
fprintf(fid, '    %-12s', 'estimator');
for i = 1:K, fprintf(fid, '     y_%d       ', i); end
fprintf(fid, '     mean      ratio\n');
for e = 1:nE
    fprintf(fid, '    %-12s', s.est_names{e});
    for i = 1:K
        fprintf(fid, ' %5.3f|%5.3f ', s.irmse(e, i), s.irmse_h2(e, i));
    end
    mh2 = mean(s.irmse_h2(e, :));
    fprintf(fid, ' %5.3f|%5.3f ', mean(s.irmse(e, :)), mh2);
    if ~isempty(ig)
        fprintf(fid, '  %6.3f', mh2 / mean(s.irmse_h2(ig, :)));
    end
    fprintf(fid, '\n');
end
if ~isempty(ig)
    fprintf(fid, '    (ratio = preferred IRMSE relative to %s; < 1 = better)\n', ...
            s.est_names{ig});
end

% --- 1b. paired comparisons against the global baseline ------------------
% A 1-2% gap in integrated RMSE cannot be judged against the naive
% Monte Carlo se of a single RMSE (rmse/sqrt(2R)); every estimator sees
% the SAME simulated samples, so the DIFFERENCE is estimated far more
% precisely than either level.  This is that difference, with its own
% standard error.
if ~isempty(ig) && H >= 2
    fprintf(fid, '\n[1b] PAIRED COMPARISON against %s, h = 2..%d (same samples)\n', ...
            s.est_names{ig}, H);
    fprintf(fid, '     %-12s %11s %10s %7s %14s %8s\n', 'estimator', ...
            'dMSE', 'se', 't', 'RMSE ratio', 'corr');
    for e = 1:nE
        if e == ig, continue; end
        try
            P = paired_comparison(mc_for(s), s.est_names{e}, s.est_names{ig}, 2:H);
        catch
            continue
        end
        fprintf(fid, '     %-12s %11.2e %10.2e %7.2f  %.3f [%.3f,%.3f] %6.2f\n', ...
                s.est_names{e}, P.dMSE, P.se, P.t, P.rmse_ratio, ...
                P.rmse_ratio_lo, P.rmse_ratio_hi, P.correlation);
    end
    fprintf(fid, ['     dMSE < 0 means lower MSE than the baseline; the bracket is a\n' ...
                  '     one-standard-error band on the RMSE ratio.  |t| < 2 means the\n' ...
                  '     difference is not resolved even with %d paired replications.\n' ...
                  '     (This RMSE ratio is sqrt of the ratio of integrated MSEs; the\n' ...
                  '     ratio in table [1] averages per-horizon RMSEs, so the two differ\n' ...
                  '     slightly by Jensen. Both are reported; neither is wrong.)\n'], R);

    % --- THE table: the ranking depends on which horizons are averaged --
    % Adaptation buys bias reduction early and pays variance late, so the
    % integrated verdict is decided by how many late horizons the average
    % happens to include.  Reporting a single integrated number without
    % this decomposition invites reading an arbitrary choice of H as a
    % property of the estimator.
    wins = {2:min(6, H), 2:min(12, H), 2:H, min(7, H):H};
    lbl  = {'early', 'to h=12', 'ALL h>=2', 'late'};
    fprintf(fid, '\n[1c] RMSE RATIO BY HORIZON WINDOW (paired; < 1 = adaptation helps)\n');
    fprintf(fid, '     %-12s', 'estimator');
    for w = 1:numel(wins)
        fprintf(fid, ' %16s', sprintf('%s (%d-%d)', lbl{w}, wins{w}(1), wins{w}(end)));
    end
    fprintf(fid, '\n');
    for e = 1:nE
        if e == ig, continue; end
        if isempty(strfind(s.est_names{e}, 'BLP')), continue; end  %#ok<STREMP>
        fprintf(fid, '     %-12s', s.est_names{e});
        for w = 1:numel(wins)
            try
                P = paired_comparison(mc_for(s), s.est_names{e}, s.est_names{ig}, wins{w});
                fprintf(fid, ' %8.3f (t%5.1f)', P.rmse_ratio, P.t);
            catch
                fprintf(fid, ' %16s', '-');
            end
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, ['     A sign flip across the columns means the integrated ranking is an\n' ...
                  '     artefact of the reporting horizon, not a property of the estimator.\n']);
end

% --- 2. bias-variance ----------------------------------------------------
fprintf(fid, '\n[2] BIAS-VARIANCE DECOMPOSITION (averaged over responses)\n');
if s.early_H >= H
    late_lbl = '(no late window)';
else
    late_lbl = sprintf('h = %d..%d (late)', s.early_H + 1, H);
end
fprintf(fid, '    %-12s %-22s %-22s %-22s\n', '', ...
        sprintf('h = 2..%d', H), sprintf('h = 2..%d (early)', min(s.early_H, H)), ...
        late_lbl);
fprintf(fid, '    %-12s %7s %7s %6s %7s %7s %6s %7s %7s %6s\n', 'estimator', ...
        '|bias|', 'var', 'rmse', '|bias|', 'var', 'rmse', '|bias|', 'var', 'rmse');
h_e = 2:min(s.early_H, H);
h_l = min(s.early_H, H) + 1:H;
for e = 1:nE
    fprintf(fid, '    %-12s', s.est_names{e});
    for w = {2:H, h_e, h_l}
        idx = w{1};
        if isempty(idx)
            fprintf(fid, ' %7s %7s %6s', '-', '-', '-');
        else
            b = mean(mean(abs(s.bias(e, :, idx))));
            v = mean(mean(s.variance(e, :, idx)));
            rm = mean(mean(s.rmse(e, :, idx)));
            fprintf(fid, ' %7.4f %7.4f %6.3f', b, v, rm);
        end
    end
    fprintf(fid, '\n');
end
fprintf(fid, '    (variance is the Monte Carlo variance of the estimator;\n');
fprintf(fid, '     MSE = bias^2 + variance exactly, asserted in summarize_montecarlo)\n');

% --- per-response bias/variance where the misspecification bites --------
if ~isempty(s.misspec_block)
    fprintf(fid, '\n    Per-response detail, h = 2..%d:\n', H);
    fprintf(fid, '    %-12s', 'estimator');
    for i = 1:K, fprintf(fid, '   y_%d |bias| var  ', i); end
    fprintf(fid, '\n');
    for e = 1:nE
        fprintf(fid, '    %-12s', s.est_names{e});
        for i = 1:K
            fprintf(fid, '  %6.4f %6.4f ', ...
                    mean(abs(s.bias(e, i, 2:H))), mean(s.variance(e, i, 2:H)));
        end
        fprintf(fid, '\n');
    end
end

% --- 3. intervals ---------------------------------------------------------
lvl = 0.90;
if isfield(s, 'meta') && isfield(s.meta, 'ci_level'), lvl = s.meta.ci_level; end
fprintf(fid, '\n[3] INTERVALS at nominal %.0f%%, averaged over responses and h = 2..%d\n', ...
        100 * lvl, H);
fprintf(fid, '    %-12s %10s %10s %14s %12s\n', 'estimator', 'coverage', ...
        'avg length', 'cover (post)', 'len (post)');
for e = 1:nE
    cp = NaN;  ap = NaN;
    if isfield(s, 'coverage_post')
        v = reshape(s.coverage_post(e, :, 2:H), [], 1);  v = v(~isnan(v));
        if ~isempty(v), cp = mean(v); end
        v = reshape(s.avg_len_post(e, :, 2:H), [], 1);   v = v(~isnan(v));
        if ~isempty(v), ap = mean(v); end
    end
    fprintf(fid, '    %-12s %10.3f %10.3f %14s %12s\n', s.est_names{e}, ...
            mean(mean(s.coverage(e, :, 2:H))), mean(mean(s.avg_len(e, :, 2:H))), ...
            num_or_dash(cp), num_or_dash(ap));
end

% --- 4. tau diagnostics ---------------------------------------------------
if isfield(s, 'tau_stats') && isstruct(s.tau_stats)
    names = fieldnames(s.tau_stats);
    for n = 1:numel(names)
        d = s.tau_stats.(names{n});
        fprintf(fid, '\n[4%s] TAU DIAGNOSTICS -- %s scales\n', char('a' + n - 1), d.est);
        if isfield(d, 'no_unique_block') && d.no_unique_block
            fprintf(fid, ['    This DGP IS misspecified but has NO single block to find\n' ...
                          '    (dense/VARMA): the flags below are neither detections nor\n' ...
                          '    false positives.\n']);
        elseif d.false_positive
            if isfield(d, 'fitted_p_nests_truth') && d.fitted_p_nests_truth
                fprintf(fid, ['    The FITTED lag order nests this DGP''s truth: nothing is\n' ...
                              '    misspecified here, so every flag below is a FALSE POSITIVE.\n']);
            else
                fprintf(fid, '    Correctly specified DGP: every flag below is a FALSE POSITIVE.\n');
            end
        else
            fprintf(fid, '    Truly misspecified block: g* = %d\n', d.misspec_block);
        end
        fprintf(fid, '    mean posterior tau by block (averaged over equations, h = 1..%d):\n', d.H);
        fprintf(fid, '      ');
        for g = 1:d.G
            fprintf(fid, 'block %d: %5.2f   ', g, mean(mean(d.tau_bar(:, g, :))));
        end
        fprintf(fid, '\n');
        fprintf(fid, '    mean posterior P(tau > 1) by block:\n      ');
        for g = 1:d.G
            fprintf(fid, 'block %d: %5.2f   ', g, mean(mean(d.p_gt1_bar(:, g, :))));
        end
        fprintf(fid, '\n');
        fprintf(fid, '    argmax frequency by block (aggregate): ');
        fprintf(fid, '%5.2f ', d.argmax_freq);
        fprintf(fid, '   (chance = %.2f)\n', 1 / d.G);
        % The per-equation winner map is the object the diagnostic
        % reading rests on.  On the correct DGP it IS the per-equation
        % false-positive map, so print it in both cases rather than only
        % when there is a true block to compare against.
        fprintf(fid, '    argmax frequency by equation x block (all horizons):\n');
        for i = 1:d.K
            fprintf(fid, '      eq %d: ', i);
            fprintf(fid, '%5.2f ', d.argmax_freq_by_eq(i, :));
            fprintf(fid, '\n');
        end
        fprintf(fid, '    the same on early horizons (h <= %d):\n', d.early_H);
        for i = 1:d.K
            fprintf(fid, '      eq %d: ', i);
            fprintf(fid, '%5.2f ', d.argmax_freq_by_eq_early(i, :));
            fprintf(fid, '\n');
        end
        if isfield(d, 'max_argmax_freq')
            fprintf(fid, ['    CONCENTRATION (max over equation x block of the win ' ...
                          'rate; nominal chance = %.2f, but the honest null is\n' ...
                          '     the same statistic on the CORRECT-DGP run -- it is a ' ...
                          'maximum over %d cells):\n'], d.argmax_chance, d.K * d.G);
            fprintf(fid, '      all horizons : %.2f  at (eq %d, block %d)\n', ...
                    d.max_argmax_freq, d.argmax_cell(1), d.argmax_cell(2));
            fprintf(fid, '      early (h<=%d) : %.2f  at (eq %d, block %d)', ...
                    d.early_H, d.max_argmax_freq_early, ...
                    d.argmax_cell_early(1), d.argmax_cell_early(2));
            if d.false_positive
                fprintf(fid, '   <- FALSE-POSITIVE baseline\n');
            elseif ~isempty(d.misspec_block)
                if d.argmax_cell_early(2) == d.misspec_block
                    fprintf(fid, '   <- the TRUE block\n');
                else
                    fprintf(fid, '   <- NOT the true block (g* = %d)\n', d.misspec_block);
                end
            else
                fprintf(fid, '\n');
            end
        end
        if ~isnan(d.detect_prob)
            fprintf(fid, '    P(block g* ranked first): aggregate %.2f;  by equation:', d.detect_prob);
            fprintf(fid, ' %.2f', d.detect_prob_by_eq);  fprintf(fid, '\n');
            fprintf(fid, '    early horizons (h <= %d), by equation:      ', d.early_H);
            fprintf(fid, ' %.2f', d.detect_prob_by_eq_early);  fprintf(fid, '\n');
            fprintf(fid, '    mean rank of g* : %.2f (chance = %.2f, best = 1)\n', ...
                    d.mean_rank_true, d.rank_chance);
        end
        fprintf(fid, ['    NOTE: the aggregate flag rules below average tau over ' ...
                  'equations before\n          comparing blocks, which dilutes a ' ...
                  'signal confined to one equation.\n          Compare their rates ' ...
                  'against the correct-DGP run before using them.\n']);
        fprintf(fid, '    flag rate at ratio thresholds');
        fprintf(fid, ' %.2f:', d.flag_thresholds);  fprintf(fid, '  ');
        fprintf(fid, '%.2f ', d.flag_ratio);  fprintf(fid, '\n');
        if ~all(isnan(d.flag_prob))
            fprintf(fid, '    flag rate at P(tau>1) thresholds');
            fprintf(fid, ' %.2f:', d.prob_thresholds);  fprintf(fid, '  ');
            fprintf(fid, '%.2f ', d.flag_prob);  fprintf(fid, '\n');
        end
        fprintf(fid, '    chain: mean ESS(beta) %s, ESS(log tau) %s, acceptance %s, tau clips/rep %s\n', ...
                num_or_dash(d.ess_beta_mean), num_or_dash(d.ess_logtau_mean), ...
                num_or_dash(d.acc_rate_mean), num_or_dash(d.n_tau_clip_mean));
    end
end
fprintf(fid, '\n');
end

% =====================================================================
function s = num_or_dash(v)
if isnan(v), s = '-'; else, s = sprintf('%.3f', v); end
end

function mc = mc_for(s)
% The paired comparison needs the raw replication output.  It is carried
% on the summary as .mc_ref when the summary was produced from a file or
% handed one; without it the paired table is skipped rather than faked.
if isfield(s, 'mc_ref') && ~isempty(s.mc_ref)
    mc = s.mc_ref;
else
    error('report_montecarlo: no replication output available for pairing.');
end
end
