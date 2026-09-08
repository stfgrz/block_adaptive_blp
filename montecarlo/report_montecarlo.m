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
    if isfield(L, 's'), s = L.s; else, s = summarize_montecarlo(L.mc); end
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

% --- 2. bias-variance ----------------------------------------------------
fprintf(fid, '\n[2] BIAS-VARIANCE DECOMPOSITION (averaged over responses)\n');
fprintf(fid, '    %-12s %-22s %-22s %-22s\n', '', ...
        sprintf('h = 2..%d', H), sprintf('h = 2..%d (early)', s.early_H), ...
        sprintf('h = %d..%d (late)', s.early_H + 1, H));
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
        if d.false_positive
            fprintf(fid, '    Correctly specified DGP: every flag below is a FALSE POSITIVE.\n');
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
        fprintf(fid, '%5.2f ', d.argmax_freq);  fprintf(fid, '\n');
        if ~isnan(d.detect_prob)
            fprintf(fid, '    P(block g* ranked first): aggregate %.2f;  by equation:', d.detect_prob);
            fprintf(fid, ' %.2f', d.detect_prob_by_eq);  fprintf(fid, '\n');
            fprintf(fid, '    early horizons (h <= %d), by equation:      ', d.early_H);
            fprintf(fid, ' %.2f', d.detect_prob_by_eq_early);  fprintf(fid, '\n');
            fprintf(fid, '    mean rank of g* : %.2f (chance = %.2f, best = 1)\n', ...
                    d.mean_rank_true, d.rank_chance);
        end
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
