function T = collect_grid_results(in_dir, out_csv)
% PURPOSE
% -------
% Turn a folder of exploratory grid results into ONE comparison table
% answering the question the grid was built for: in which cells does
% block adaptation lower integrated RMSE, in which does it only help at
% early horizons, and in which is it useful only as a diagnostic?
%
% THE THREE VERDICTS
% ------------------
% For each cell, with r_all  = IRMSE(adaptive, h = 2..H) / IRMSE(global),
%                    r_early = IRMSE(adaptive, h = 2..early_H) / IRMSE(global):
%   'improves RMSE'     r_all   <= 1 - margin
%   'early horizons'    r_all   >  1 - margin  and  r_early <= 1 - margin
%   'diagnostic only'   neither ratio improves, but tau still localises
%                       the misspecified block better than chance
%   'no signal'         no RMSE gain and no localisation
% margin defaults to 0.02 (a 2% RMSE gap), which is BELOW the Monte Carlo
% noise of the exploratory preset (about 11% at R = 40 for a single RMSE,
% less for a ratio computed under common random numbers).  The verdict
% column is therefore a SORTING device for choosing R = 500 cells, not a
% finding; the table also carries the raw ratios so a reader can apply
% their own threshold.
%
% INPUTS
% ------
% in_dir  : folder of grid_*.mat files (default <repo>/results/grid).
% out_csv : where to write the table (default <in_dir>/grid_summary.csv).
%
% OUTPUTS
% -------
% T : struct array, one entry per cell, with the fields written to the csv.

if nargin < 1 || isempty(in_dir)
    root = fileparts(fileparts(mfilename('fullpath')));
    in_dir = fullfile(root, 'results', 'grid');
end
if nargin < 2 || isempty(out_csv)
    out_csv = fullfile(in_dir, 'grid_summary.csv');
end
margin = 0.02;

d = dir(fullfile(in_dir, 'grid_*.mat'));
assert(~isempty(d), 'collect_grid_results: no grid_*.mat in %s', in_dir);

T = struct('name', {}, 'axis', {}, 'dgp', {}, 'label', {}, 'R', {}, ...
           'T_obs', {}, 'p', {}, 'H', {}, 'irmse_global', {}, ...
           'irmse_block', {}, 'irmse_pooled', {}, 'ratio_block', {}, ...
           'ratio_pooled', {}, 'ratio_block_early', {}, ...
           'ratio_pooled_early', {}, 'detect_prob', {}, ...
           'detect_prob_best_eq', {}, 'flag_rate', {}, 'verdict_block', {}, ...
           'verdict_pooled', {});

for k = 1:numel(d)
    L = load(fullfile(in_dir, d(k).name));
    s = L.s;  cs = L.case_spec;
    ig = find(strcmp(s.est_names, 'BLP-FMAR'), 1);
    ib = find(strcmp(s.est_names, 'BLP-block'), 1);
    ip = find(strcmp(s.est_names, 'BLP-pooled'), 1);
    assert(~isempty(ig) && ~isempty(ib), ...
        'collect_grid_results: %s lacks the global or block estimator', d(k).name);

    g_all = mean(s.irmse_h2(ig, :));   g_early = mean(s.irmse_early(ig, :));
    b_all = mean(s.irmse_h2(ib, :));   b_early = mean(s.irmse_early(ib, :));
    if ~isempty(ip)
        p_all = mean(s.irmse_h2(ip, :));  p_early = mean(s.irmse_early(ip, :));
    else
        p_all = NaN;  p_early = NaN;
    end

    dg = [];
    if isfield(s, 'tau_stats') && isfield(s.tau_stats, 'block')
        dg = s.tau_stats.block;
    end
    [det_prob, det_best, flagr, localises] = localisation(dg);

    e = numel(T) + 1;
    T(e).name = cs.name;  T(e).axis = cs.axis;  T(e).dgp = cs.dgp;
    T(e).label = cs.label;
    T(e).R = s.R;  T(e).T_obs = s.meta.T;  T(e).p = s.meta.p;  T(e).H = s.H;
    T(e).irmse_global = g_all;  T(e).irmse_block = b_all;  T(e).irmse_pooled = p_all;
    T(e).ratio_block  = b_all / g_all;
    T(e).ratio_pooled = p_all / g_all;
    T(e).ratio_block_early  = b_early / g_early;
    T(e).ratio_pooled_early = p_early / g_early;
    T(e).detect_prob = det_prob;
    T(e).detect_prob_best_eq = det_best;
    T(e).flag_rate = flagr;
    T(e).verdict_block  = verdict(b_all / g_all, b_early / g_early, localises, margin);
    T(e).verdict_pooled = verdict(p_all / g_all, p_early / g_early, localises, margin);
end

[~, ord] = sort({T.name});
T = T(ord);

fid = fopen(out_csv, 'w');
assert(fid > 0, 'collect_grid_results: cannot write %s', out_csv);
fprintf(fid, ['case,axis,dgp,R,T,p,H,irmse_global_h2,irmse_block_h2,' ...
              'irmse_pooled_h2,ratio_block,ratio_pooled,ratio_block_early,' ...
              'ratio_pooled_early,detect_prob,detect_prob_best_eq,flag_rate,' ...
              'verdict_block,verdict_pooled,label\n']);
for e = 1:numel(T)
    fprintf(fid, '%s,%s,%s,%g,%g,%g,%g,%.6g,%.6g,%.6g,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f,%.3f,%s,%s,"%s"\n', ...
        T(e).name, T(e).axis, T(e).dgp, T(e).R, T(e).T_obs, T(e).p, T(e).H, ...
        T(e).irmse_global, T(e).irmse_block, T(e).irmse_pooled, ...
        T(e).ratio_block, T(e).ratio_pooled, T(e).ratio_block_early, ...
        T(e).ratio_pooled_early, T(e).detect_prob, T(e).detect_prob_best_eq, ...
        T(e).flag_rate, T(e).verdict_block, T(e).verdict_pooled, T(e).label);
end
fclose(fid);

fprintf('\n%-22s %-12s %7s %7s %7s %7s  %-16s\n', 'case', 'axis', ...
        'r_blk', 'r_pool', 'rE_blk', 'rE_pool', 'verdict (block)');
for e = 1:numel(T)
    fprintf('%-22s %-12s %7.3f %7.3f %7.3f %7.3f  %-16s\n', T(e).name, T(e).axis, ...
            T(e).ratio_block, T(e).ratio_pooled, T(e).ratio_block_early, ...
            T(e).ratio_pooled_early, T(e).verdict_block);
end
fprintf('\ncollect_grid_results: wrote %s (%d cells)\n', out_csv, numel(T));
fprintf(['Ratios are adaptive / global integrated RMSE; < 1 means adaptation ' ...
         'helps.\nVerdicts use a %.0f%% margin and are a sorting device at ' ...
         'R = %g, not a finding.\n'], 100 * margin, T(1).R);
end

% =====================================================================
function [det_prob, det_best, flagr, localises] = localisation(dg)
det_prob = NaN;  det_best = NaN;  flagr = NaN;  localises = false;
if isempty(dg), return; end
flagr = dg.flag_rate;
if ~isnan(dg.detect_prob)
    det_prob = dg.detect_prob;
    det_best = max(dg.detect_prob_by_eq);
    % "Better than chance" = the true block is ranked first in some
    % equation more often than a uniform pick over the G blocks would.
    localises = det_best > 1.5 / dg.G;
end
end

function v = verdict(r_all, r_early, localises, margin)
if isnan(r_all)
    v = 'not run';
elseif r_all <= 1 - margin
    v = 'improves RMSE';
elseif r_early <= 1 - margin
    v = 'early horizons';
elseif localises
    v = 'diagnostic only';
else
    v = 'no signal';
end
end
