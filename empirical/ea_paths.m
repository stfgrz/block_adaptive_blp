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
%     .legacy     <root>/empirical/legacy_v1          (the 2026-09-12 internal-
%                 instrument design: drivers and data builders, kept runnable)
%     .results    <root>/results/empirical            (current design: iv_*,
%                 nulls, ablations; figures in .figures)
%     .figures    <root>/results/empirical/figures
%     .results_legacy, .figures_legacy   the same for the v1 outputs
%                 (<root>/results/empirical/legacy_v1[/figures])
%     .results_sim  <root>/results/simulation         (Monte Carlo output;
%                 written by the root-level drivers, listed here for reference)
%     .dataset    <root>/empirical/data/ea_dataset_v2_ois4.mat  (the headline
%                 system of the current design; built by assemble_dataset_v2)
%     .dataset_legacy  <root>/empirical/data/ea_dataset.mat  (v1 five-variable
%                 system with the surprise inside the VAR)
%
% NOTES
% -----
% * .results and the data folders are NOT created here; each writer
%   creates the directory it is about to write to.
% * Every entry point of the package (RUN_EMPIRICAL_IV, the data builders,
%   run_null_calibration, run_block_ablation and the legacy_v1 drivers)
%   starts with a two-line bootstrap that puts the
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
P.legacy    = fullfile(here, 'legacy_v1');
P.results   = fullfile(P.root, 'results', 'empirical');
P.figures   = fullfile(P.results, 'figures');
P.results_legacy = fullfile(P.results, 'legacy_v1');
P.figures_legacy = fullfile(P.results_legacy, 'figures');
P.results_sim    = fullfile(P.root, 'results', 'simulation');
P.dataset        = fullfile(P.data, 'ea_dataset_v2_ois4.mat');
P.dataset_legacy = fullfile(P.data, 'ea_dataset.mat');
end
