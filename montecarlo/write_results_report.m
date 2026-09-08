function out_file = write_results_report(out_file, results_dir)
% PURPOSE
% -------
% Write ONE plain-text report containing every headline table of the
% study: the full report_montecarlo output for each stored Monte Carlo
% result, followed by the exploratory grid summary.
%
% WHY
% ---
% Numbers that are transcribed by hand into a document drift from the
% numbers in the files.  This regenerates the whole report from the
% stored results in one command, so `docs/RESULTS.md` can quote a few
% headline figures and point at a file that is always current rather
% than restating dozens of numbers that might not be.
%
% USAGE
% -----
%   write_results_report                       % -> results/REPORT.txt
%   write_results_report('somewhere/else.txt')
%
% INPUTS
% ------
% out_file    : destination (default <repo>/results/REPORT.txt).
% results_dir : where to look for mc_*.mat and grid/ (default
%               <repo>/results).
%
% OUTPUTS
% -------
% out_file : the path written.
%
% NOTES
% -----
% Results are reported in a fixed order (final runs first, by DGP, then
% anything else) so successive versions of the file are comparable line
% by line.  A summary saved before the current metrics existed is
% recomputed from its stored replication output.

root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(out_file)
    out_file = fullfile(root, 'results', 'REPORT.txt');
end
if nargin < 2 || isempty(results_dir)
    results_dir = fullfile(root, 'results');
end

d = dir(fullfile(results_dir, 'mc_*.mat'));
names = {d.name};
% Fixed order: the final runs first, by DGP, then everything else.
pref = {'mc_final_sparse.mat', 'mc_final_correct.mat', ...
        'mc_final_intermediate.mat', 'mc_final_dense.mat'};
ordered = {};
for k = 1:numel(pref)
    if any(strcmp(names, pref{k})), ordered{end + 1} = pref{k}; end %#ok<AGROW>
end
rest = sort(setdiff(names, ordered));
% Skip throwaway smoke output.
rest = rest(cellfun(@(n) isempty(strfind(n, 'smoke')), rest));  %#ok<STREMP>
ordered = [ordered, rest];

outdir = fileparts(out_file);
if ~isempty(outdir) && exist(outdir, 'dir') ~= 7, mkdir(outdir); end
fid = fopen(out_file, 'w');
assert(fid > 0, 'write_results_report: cannot write %s', out_file);

fprintf(fid, 'BLOCK-ADAPTIVE BLP -- FULL RESULTS REPORT\n');
fprintf(fid, 'generated %s by montecarlo/write_results_report.m\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));  %#ok<TNOW1,DATST>
fprintf(fid, ['Read docs/RESULTS.md for how to interpret these tables, and\n' ...
              'docs/APPROXIMATIONS.md for what the remaining simplifications cost.\n']);

for k = 1:numel(ordered)
    f = fullfile(results_dir, ordered{k});
    fprintf(fid, '\n\n########################################################################\n');
    fprintf(fid, '# %s\n', ordered{k});
    fprintf(fid, '########################################################################\n');
    try
        report_montecarlo(f, fid);
    catch err
        fprintf(fid, '  [could not report: %s]\n', err.message);
    end
    fprintf('write_results_report: %s\n', ordered{k});
end

% --- exploratory grid ----------------------------------------------------
gcsv = fullfile(results_dir, 'grid', 'grid_summary.csv');
if exist(gcsv, 'file') == 2
    fprintf(fid, '\n\n########################################################################\n');
    fprintf(fid, '# EXPLORATORY GRID (R = 40 -- a sorting device, not a finding)\n');
    fprintf(fid, '########################################################################\n');
    fprintf(fid, ['Ratios are adaptive / global integrated RMSE over h = 2..H;\n' ...
                  '< 1 means adaptation helps.  At R = 40 the Monte Carlo se of a\n' ...
                  'single RMSE is about 11%% of it.\n\n']);
    g = fopen(gcsv, 'r');
    while true
        l = fgetl(g);
        if ~ischar(l), break; end
        fprintf(fid, '%s\n', l);
    end
    fclose(g);
end

% --- sensitivity ---------------------------------------------------------
scsv = fullfile(results_dir, 'sensitivity_approximations.csv');
if exist(scsv, 'file') == 2
    fprintf(fid, '\n\n########################################################################\n');
    fprintf(fid, '# SENSITIVITY OF THE METHODOLOGICAL APPROXIMATIONS\n');
    fprintf(fid, '########################################################################\n');
    s = fopen(scsv, 'r');
    while true
        l = fgetl(s);
        if ~ischar(l), break; end
        fprintf(fid, '%s\n', l);
    end
    fclose(s);
end

fclose(fid);
fprintf('write_results_report: wrote %s\n', out_file);
end
