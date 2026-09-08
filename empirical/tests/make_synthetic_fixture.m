function ds = make_synthetic_fixture(opts)
% MAKE_SYNTHETIC_FIXTURE  Build a clearly-labelled SYNTHETIC euro-area
% dataset, for INTERFACE TESTING ONLY.
%
% WHY THIS EXISTS, AND WHAT IT MUST NEVER BE USED FOR
% ---------------------------------------------------
% The Chapter 7 pipeline needs four licensed outcome series that are not
% in the repository, so on a clean clone (or an offline machine) there is
% no ea_dataset.mat and NOTHING in the empirical package can be executed
% -- not even to check that the estimator interfaces still line up after
% a refactor.  This function produces a dataset with the right SHAPE
% (T x 5, the documented variable order, the documented isrw vector) from
% a simulated VAR, so the plumbing can be exercised.
%
% The numbers are simulated.  They are not euro-area data, they are not
% calibrated to euro-area data, and no result computed from them is an
% empirical finding.  Three safeguards keep them from leaking into the
% chapter:
%   1. every field is stamped ds.synthetic = true;
%   2. the default output path is
%      empirical/data/SYNTHETIC_ea_dataset_FIXTURE.mat, never the
%      P.dataset name that RUN_EMPIRICAL reads;
%   3. RUN_EMPIRICAL and SMOKE_TEST_EMPIRICAL refuse a dataset carrying
%      ds.synthetic = true (see the check in each).
% This package already rejected four fabricated placeholder csvs once
% (empirical/data/raw/placeholder_rejected/README.md); the fixture exists
% so that "we need something runnable" never becomes a reason to fabricate
% inputs again.
%
% MODEL
% -----
% A stable VAR(2) in the five variables, with a white-noise-ish first
% variable (the "surprise") and four persistent levels, so the system has
% the same qualitative shape as the real one:
%     y(t) = c + A1 y(t-1) + A2 y(t-2) + B0 eta(t).
% The surprise has near-zero persistence and loads on nothing; the four
% level variables are persistent and respond to it.  Nothing here is
% calibrated: it only has to be stable and non-degenerate.
%
% INPUTS (opts, all optional)
% ---------------------------
%   opts.T        number of months (default 240, the real sample length)
%   opts.seed     RNG seed (default 4242; the state is restored on exit)
%   opts.out_mat  where to save (default the SYNTHETIC_... path above);
%                 pass '' to skip saving.
%
% OUTPUTS
% -------
%   ds : struct with the same fields as assemble_dataset's output plus
%        .synthetic = true and .fixture_note.
%
% NOTES
% -----
% Deterministic given opts.seed.

if nargin < 1, opts = struct(); end
if exist('ea_paths', 'file') ~= 2
    addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
P = ea_paths();
if ~isfield(opts, 'T')      || isempty(opts.T),      opts.T = 240;   end
if ~isfield(opts, 'seed')   || isempty(opts.seed),   opts.seed = 4242; end
if ~isfield(opts, 'out_mat')
    opts.out_mat = fullfile(P.data, 'SYNTHETIC_ea_dataset_FIXTURE.mat');
end

K = 5;
A1 = [ 0.05  0.00  0.00  0.00  0.00;
       0.30  0.85  0.03  0.02  0.01;
       0.20  0.08  0.70  0.04  0.03;
       0.10  0.04  0.03  0.80  0.01;
       0.50  0.03  0.02  0.02  0.10];
A2 = zeros(5);
A2(2, 2) = 0.05;  A2(3, 3) = 0.12;  A2(4, 4) = 0.08;  A2(5, 5) = 0.05;
% max |companion eigenvalue| = 0.953: persistent like the real system but
% comfortably stable (asserted below).
c  = [0; 0.02; 0.05; 0.03; 0.04];
B0 = diag([0.05 0.15 0.60 0.20 3.00]);
B0(3, 1) = 0.10;  B0(5, 1) = 0.50;      % lower triangular: shock loads on levels

F = [A1, A2; eye(K), zeros(K)];
assert(max(abs(eig(F))) < 0.98, ...
    'make_synthetic_fixture: fixture VAR is not comfortably stable.');

T = opts.T;  burn = 200;
st = rng();  cleanup = onCleanup(@() rng(st));   %#ok<NASGU>
rng(opts.seed, 'twister');
Yfull = zeros(T + burn, K);
for t = 3:(T + burn)
    Yfull(t, :) = (c + A1 * Yfull(t-1, :)' + A2 * Yfull(t-2, :)' + ...
                   B0 * randn(K, 1))';
end
Y = Yfull(burn + 1:end, :);
% Put the level variables on plausible-looking scales (still synthetic).
Y(:, 2) = Y(:, 2) + 2;          % "rate", percent
Y(:, 3) = Y(:, 3) + 450;        % "100*log IP"
Y(:, 4) = Y(:, 4) + 460;        % "100*log HICP"
Y(:, 5) = Y(:, 5) + 800;        % "100*log STOXX"

ds.Y = Y;
ds.varnames = {'mps', 'i1y', 'ip', 'hicp', 'stoxx'};
ds.ym = (12 * 2000 + 1 : 12 * 2000 + T)';
ds.isrw = [0 1 1 1 1];
ds.shock_variant = 'SYNTHETIC_FIXTURE';
ds.synthetic = true;
ds.fixture_note = ['SIMULATED DATA -- interface testing only. Not euro-area ' ...
    'data and not calibrated to it. No result computed from this file is ' ...
    'an empirical finding. Built by empirical/tests/make_synthetic_fixture.m.'];
ds.meta = struct('generator', 'make_synthetic_fixture', 'seed', opts.seed, ...
                 'T', T, 'K', K);

if ~isempty(opts.out_mat)
    outdir = fileparts(opts.out_mat);
    if ~isempty(outdir) && exist(outdir, 'dir') ~= 7, mkdir(outdir); end
    save(opts.out_mat, '-struct', 'ds');
    fprintf('make_synthetic_fixture: wrote SYNTHETIC fixture %s (T = %d)\n', ...
            opts.out_mat, T);
end
end
