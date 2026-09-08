function files = export_montecarlo_csv(s, out_stem)
% PURPOSE
% -------
% Write a Monte Carlo summary to TIDY CSV files so that results can be
% inspected, plotted or tabulated outside MATLAB/Octave without anyone
% having to reload a .mat file and remember its array layout.  Together
% with mc.meta this is what makes a stored result self-describing.
%
% FILES WRITTEN (out_stem = ".../mc_fmar_sparse")
% -----------------------------------------------
%   <stem>_metrics.csv      one row per (estimator, variable, horizon):
%       dgp, mode, R, T, p, H, estimator, variable, horizon,
%       bias, variance, mse, rmse, coverage, avg_len,
%       coverage_post, avg_len_post   (the posterior-quantile bands of
%       the sampled estimators; NaN for the closed-form ones)
%   <stem>_integrated.csv   one row per (estimator, variable):
%       dgp, estimator, variable, irmse_legacy_h1, irmse_h2,
%       irmse_early, irmse_late, ibias_h2, ivar_h2
%   <stem>_tau.csv          one row per (scale estimator, equation,
%                           block, horizon):
%       dgp, scale_estimator, equation, block, horizon,
%       tau_mean, tau_median, p_tau_gt1
%   <stem>_tau_localization.csv  one row per (scale estimator, metric):
%       dgp, scale_estimator, metric, index, value
%       (argmax frequencies, detection probabilities, mean ranks,
%        flag/false-positive rates, chain diagnostics)
%   <stem>_meta.txt         the reproducibility record in key = value
%                           form (DGP parameters, seed rule, R/T/p/H,
%                           estimator settings, sampler sizes).
%
% ALL numbers are printed with %.10g, i.e. losslessly enough to
% reproduce every table in the write-up from the CSVs alone.
%
% INPUTS
% ------
% s        : output of summarize_montecarlo (must carry .meta to write
%            the meta file; the other files need only the metrics).
% out_stem : path stem WITHOUT extension.  Directories are created.
%
% OUTPUTS
% -------
% files : cellstr of the paths written.
%
% NOTES
% -----
% Uses only fopen/fprintf, so it needs no toolbox and behaves
% identically in MATLAB and Octave.

assert(isstruct(s) && isfield(s, 'rmse'), ...
    'export_montecarlo_csv: expects a summarize_montecarlo output.');
outdir = fileparts(out_stem);
if ~isempty(outdir) && exist(outdir, 'dir') ~= 7
    [ok, msg] = mkdir(outdir);
    assert(ok == 1, 'export_montecarlo_csv: cannot create %s (%s)', outdir, msg);
end

[nE, K, H] = size(s.rmse);
dgp = s.dgp_name;
if isfield(s, 'meta')
    md = s.meta;
else
    md = struct('mode', 'unknown', 'R', s.R, 'T', NaN, 'p', NaN, 'H', H);
end
mode_str = getfield_default(md, 'mode', 'unknown');
Rv = getfield_default(md, 'R', s.R);
Tv = getfield_default(md, 'T', NaN);
pv = getfield_default(md, 'p', NaN);

files = {};

% --- per-horizon metrics ------------------------------------------------
f1 = [out_stem '_metrics.csv'];
fid = fopen(f1, 'w');
assert(fid > 0, 'export_montecarlo_csv: cannot write %s', f1);
fprintf(fid, ['dgp,mode,R,T,p,H,estimator,variable,horizon,bias,variance,' ...
              'mse,rmse,coverage,avg_len,coverage_post,avg_len_post\n']);
for e = 1:nE
    for i = 1:K
        for h = 1:H
            cp = NaN;  ap = NaN;
            if isfield(s, 'coverage_post')
                cp = s.coverage_post(e, i, h);  ap = s.avg_len_post(e, i, h);
            end
            fprintf(fid, '%s,%s,%g,%g,%g,%g,%s,%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', ...
                dgp, mode_str, Rv, Tv, pv, H, s.est_names{e}, i, h, ...
                s.bias(e, i, h), s.variance(e, i, h), s.mse(e, i, h), ...
                s.rmse(e, i, h), s.coverage(e, i, h), s.avg_len(e, i, h), cp, ap);
        end
    end
end
fclose(fid);  files{end + 1} = f1;

% --- integrated metrics --------------------------------------------------
f2 = [out_stem '_integrated.csv'];
fid = fopen(f2, 'w');
assert(fid > 0, 'export_montecarlo_csv: cannot write %s', f2);
fprintf(fid, ['dgp,mode,R,T,p,H,estimator,variable,irmse_legacy_h1,irmse_h2,' ...
              'irmse_early,irmse_late,ibias_h2,ivar_h2\n']);
for e = 1:nE
    for i = 1:K
        fprintf(fid, '%s,%s,%g,%g,%g,%g,%s,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', ...
            dgp, mode_str, Rv, Tv, pv, H, s.est_names{e}, i, ...
            s.irmse(e, i), s.irmse_h2(e, i), s.irmse_early(e, i), ...
            s.irmse_late(e, i), s.ibias_h2(e, i), s.ivar_h2(e, i));
    end
end
fclose(fid);  files{end + 1} = f2;

% --- tau surfaces and localisation --------------------------------------
if isfield(s, 'tau_stats') && isstruct(s.tau_stats) && ...
        ~isempty(fieldnames(s.tau_stats))
    names = fieldnames(s.tau_stats);

    f3 = [out_stem '_tau.csv'];
    fid = fopen(f3, 'w');
    assert(fid > 0, 'export_montecarlo_csv: cannot write %s', f3);
    fprintf(fid, 'dgp,scale_estimator,equation,block,horizon,tau_mean,tau_median,p_tau_gt1\n');
    for n = 1:numel(names)
        d = s.tau_stats.(names{n});
        for i = 1:d.K
            for g = 1:d.G
                for h = 1:d.H
                    fprintf(fid, '%s,%s,%d,%d,%d,%.10g,%.10g,%.10g\n', ...
                        dgp, d.est, i, g, h, d.tau_bar(i, g, h), ...
                        d.tau_med_bar(i, g, h), d.p_gt1_bar(i, g, h));
                end
            end
        end
    end
    fclose(fid);  files{end + 1} = f3;

    f4 = [out_stem '_tau_localization.csv'];
    fid = fopen(f4, 'w');
    assert(fid > 0, 'export_montecarlo_csv: cannot write %s', f4);
    fprintf(fid, 'dgp,scale_estimator,misspec_block,is_false_positive_setting,metric,index,value\n');
    for n = 1:numel(names)
        d = s.tau_stats.(names{n});
        gs = d.misspec_block;
        if isempty(gs), gs = NaN; end
        w = @(metric, idx, val) fprintf(fid, '%s,%s,%g,%d,%s,%s,%.10g\n', ...
            dgp, d.est, gs, d.false_positive, metric, idx, val);
        for g = 1:d.G
            w('argmax_freq', sprintf('block%d', g), d.argmax_freq(g));
        end
        for i = 1:d.K
            for g = 1:d.G
                w('argmax_freq_by_eq', sprintf('eq%d_block%d', i, g), ...
                  d.argmax_freq_by_eq(i, g));
                w('argmax_freq_by_eq_early', sprintf('eq%d_block%d', i, g), ...
                  d.argmax_freq_by_eq_early(i, g));
            end
        end
        w('detect_prob', 'aggregate', d.detect_prob);
        for i = 1:d.K
            w('detect_prob_by_eq', sprintf('eq%d', i), d.detect_prob_by_eq(i));
            w('detect_prob_by_eq_early', sprintf('eq%d', i), d.detect_prob_by_eq_early(i));
            w('mean_rank_true_by_eq', sprintf('eq%d', i), d.mean_rank_true_by_eq(i));
        end
        if isfield(d, 'max_argmax_freq')
            w('max_argmax_freq', 'all_horizons', d.max_argmax_freq);
            w('max_argmax_freq_early', 'early', d.max_argmax_freq_early);
            w('argmax_cell', 'equation_all', d.argmax_cell(1));
            w('argmax_cell', 'block_all', d.argmax_cell(2));
            w('argmax_cell_early', 'equation_early', d.argmax_cell_early(1));
            w('argmax_cell_early', 'block_early', d.argmax_cell_early(2));
            w('argmax_chance', 'per_cell', d.argmax_chance);
        end
        w('mean_rank_true', 'aggregate', d.mean_rank_true);
        w('rank_chance', 'aggregate', d.rank_chance);
        for k = 1:numel(d.flag_thresholds)
            w('flag_ratio', sprintf('thr%.2f', d.flag_thresholds(k)), d.flag_ratio(k));
            w('joint_detect_ratio', sprintf('thr%.2f', d.flag_thresholds(k)), ...
              d.joint_detect_ratio(k));
        end
        for k = 1:numel(d.prob_thresholds)
            w('flag_prob', sprintf('q%.2f', d.prob_thresholds(k)), d.flag_prob(k));
            w('joint_detect_prob', sprintf('q%.2f', d.prob_thresholds(k)), ...
              d.joint_detect_prob(k));
        end
        w('ess_beta_mean', 'chain', d.ess_beta_mean);
        w('ess_logtau_mean', 'chain', d.ess_logtau_mean);
        w('lag1_beta_mean', 'chain', d.lag1_beta_mean);
        w('acc_rate_mean', 'chain', d.acc_rate_mean);
        w('n_tau_clip_mean', 'chain', d.n_tau_clip_mean);
    end
    fclose(fid);  files{end + 1} = f4;
end

% --- reproducibility record ----------------------------------------------
if isfield(s, 'meta')
    f5 = [out_stem '_meta.txt'];
    fid = fopen(f5, 'w');
    assert(fid > 0, 'export_montecarlo_csv: cannot write %s', f5);
    write_struct(fid, s.meta, '');
    fclose(fid);  files{end + 1} = f5;
end
end

% =====================================================================
function v = getfield_default(s, f, dflt)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function write_struct(fid, s, prefix)
% Flatten a (possibly nested) struct into "key = value" lines.
f = fieldnames(s);
for k = 1:numel(f)
    key = f{k};
    if ~isempty(prefix), key = [prefix '.' key]; end
    v = s.(f{k});
    if isstruct(v) && numel(v) == 1
        write_struct(fid, v, key);
    elseif ischar(v)
        fprintf(fid, '%s = %s\n', key, v);
    elseif iscellstr(v)  %#ok<ISCLSTR>
        fprintf(fid, '%s = %s\n', key, strjoin(v(:)', ', '));
    elseif islogical(v) || isnumeric(v)
        if numel(v) <= 40
            fprintf(fid, '%s = %s\n', key, num2str(double(v(:)'), '%.10g '));
        else
            fprintf(fid, '%s = [%s array, first 5: %s ...]\n', key, ...
                    mat2str(size(v)), num2str(double(v(1:5)'), '%.10g '));
        end
    else
        fprintf(fid, '%s = <%s>\n', key, class(v));
    end
end
end
