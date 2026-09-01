function P = ea_paths()
% EA_PATHS  Absolute locations used by the Ch. 7 empirical package.
%
% PURPOSE
% -------
% Every file of the empirical package resolves its inputs and outputs
% through this function, so the pipeline runs identically whatever the
% current directory is.  The locations are derived from THIS FILE's own
% path, the same idiom the rest of the repo already uses (RUN_ME_FIRST.m,
% RUN_FMAR_DEMO.m, tests/run_all_tests.m).
%
% The package deliberately lives in its own <repo>/empirical/ subtree
% rather than being unpacked over the repo root: the root already has a
% dgp/ and an estimators/ folder holding the SIMULATION machinery, and
% mixing the euro-area data-construction code into them makes both harder
% to read.  Results still land in the repo-wide <repo>/results/ folder,
% next to the Monte Carlo output.
%
% INPUTS
% ------
% (none)
%
% OUTPUTS
% -------
% P : struct of absolute paths
%     .root       repository root
%     .empirical  <root>/empirical
%     .data       <root>/empirical/data
%     .raw        <root>/empirical/data/raw       (downloaded/shipped csv)
%     .derived    <root>/empirical/data/derived   (built shock series)
%     .docs       <root>/empirical/docs
%     .results    <root>/results                  (figures, .mat, .csv)
%     .dataset    <root>/empirical/data/ea_dataset.mat  (baseline dataset)
%
% NOTES
% -----
% * .results and the data folders are NOT created here; each writer
%   creates the directory it is about to write to.
% * Every entry point of the package (RUN_EMPIRICAL, SMOKE_TEST_EMPIRICAL,
%   build_shock_series, fetch_outcome_data, assemble_dataset,
%   run_null_calibration) starts with a two-line bootstrap that puts the
%   repo on the path if this function is not visible yet, so any of them
%   can be invoked directly from any directory.  The bootstrap cannot
%   simply call a shared helper: that helper would have the same
%   visibility problem.

here = fileparts(mfilename('fullpath'));      % <repo>/empirical

P.empirical = here;
P.root      = fileparts(here);
P.data      = fullfile(here, 'data');
P.raw       = fullfile(P.data, 'raw');
P.derived   = fullfile(P.data, 'derived');
P.docs      = fullfile(here, 'docs');
P.results   = fullfile(P.root, 'results');
P.dataset   = fullfile(P.data, 'ea_dataset.mat');
end
