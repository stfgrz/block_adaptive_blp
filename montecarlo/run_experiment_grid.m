function files = run_experiment_grid(case_idx, out_dir, preset)
% PURPOSE
% -------
% Run cells of the exploratory experiment grid defined in
% montecarlo/experiment_grid_spec.m, save one .mat and one set of CSVs
% per cell, and print the headline comparison as it goes.
%
% USAGE
% -----
%   run_experiment_grid                 % every cell, serially
%   run_experiment_grid(3:5)            % cells 3 to 5 (for parallel runs)
%   run_experiment_grid([], out_dir)    % choose the output folder
%
% Each cell is independent, so N processes can each take a slice of the
% index range; nothing is shared and nothing is appended to.
%
% INPUTS
% ------
% case_idx : indices into experiment_grid_spec (default: all).
% out_dir  : output folder (default <repo>/results/grid).
% preset   : mc_preset name (default 'grid').
%
% OUTPUTS
% -------
% files : cellstr of the .mat files written.  Each contains `mc`, `s`
%         and `case_spec`; montecarlo/collect_grid_results.m turns the
%         whole folder into one comparison table.
%
% NOTES
% -----
% The printed lines report the PREFERRED integrated RMSE (h = 2..H) for
% every estimator plus the adaptive/global ratios, which is what decides
% whether a cell is worth an R = 500 run.

cases = experiment_grid_spec();
if nargin < 1 || isempty(case_idx), case_idx = 1:numel(cases); end
if nargin < 2 || isempty(out_dir)
    root = fileparts(fileparts(mfilename('fullpath')));
    out_dir = fullfile(root, 'results', 'grid');
end
if nargin < 3 || isempty(preset), preset = 'grid'; end
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end

files = {};
for k = case_idx(:)'
    cs = cases(k);
    fprintf('\n=== grid cell %d/%d: %s (%s) ===\n', k, numel(cases), ...
            cs.name, cs.label);
    cfg = mc_preset(preset, cs.overrides);
    mc = run_montecarlo(cfg, cs.dgp);
    mc.meta.preset = preset;
    mc.meta.grid_case = cs;
    s = summarize_montecarlo(mc);

    f = fullfile(out_dir, sprintf('grid_%s.mat', cs.name));
    case_spec = cs;                                     %#ok<NASGU>
    cfg = cfg_to_savable(cfg);      % handles cannot be written to a .mat
    save(f, 'mc', 's', 'case_spec', 'cfg', '-v7');
    export_montecarlo_csv(s, fullfile(out_dir, sprintf('grid_%s', cs.name)));
    files{end + 1} = f;  %#ok<AGROW>

    print_cell_summary(s, cs);
end
fprintf('\nrun_experiment_grid: %d cells written to %s\n', numel(files), out_dir);
end

% =====================================================================
function print_cell_summary(s, cs)
ig = find(strcmp(s.est_names, 'BLP-FMAR'), 1);
ib = find(strcmp(s.est_names, 'BLP-block'), 1);
ip = find(strcmp(s.est_names, 'BLP-pooled'), 1);
fprintf('  integrated RMSE, h = 2..%d (preferred metric), per response:\n', s.H);
fprintf('    %-12s', 'estimator');
for i = 1:s.K, fprintf('    y_%d ', i); end
fprintf('   mean\n');
for e = 1:numel(s.est_names)
    fprintf('    %-12s', s.est_names{e});
    fprintf(' %6.3f', s.irmse_h2(e, :));
    fprintf('  %6.3f\n', mean(s.irmse_h2(e, :)));
end
if ~isempty(ig) && ~isempty(ib)
    fprintf('  adaptive/global IRMSE ratio (h>=2): block %.3f', ...
            mean(s.irmse_h2(ib, :)) / mean(s.irmse_h2(ig, :)));
    if ~isempty(ip)
        fprintf(', pooled %.3f', mean(s.irmse_h2(ip, :)) / mean(s.irmse_h2(ig, :)));
    end
    fprintf('   (< 1 = adaptation helps)\n');
    fprintf('  early horizons only (h = 2..%d): block %.3f', s.early_H, ...
            mean(s.irmse_early(ib, :)) / mean(s.irmse_early(ig, :)));
    if ~isempty(ip)
        fprintf(', pooled %.3f', mean(s.irmse_early(ip, :)) / mean(s.irmse_early(ig, :)));
    end
    fprintf('\n');
end
if isfield(s, 'tau_stats') && isfield(s.tau_stats, 'block')
    d = s.tau_stats.block;
    if ~isnan(d.detect_prob)
        fprintf('  tau localisation: P(true block %d ranked first) = %.2f aggregate, %.2f in eq %d\n', ...
                d.misspec_block, d.detect_prob, ...
                max(d.detect_prob_by_eq), find(d.detect_prob_by_eq == max(d.detect_prob_by_eq), 1));
    else
        fprintf('  tau false-positive rate (flag ratio > 1.5): %.2f\n', d.flag_rate);
    end
end
fprintf('  [%s axis: %s]\n', cs.axis, cs.label);
end
