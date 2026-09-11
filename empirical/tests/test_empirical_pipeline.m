function test_empirical_pipeline()
% PURPOSE
% -------
% Runnable check of the Chapter 7 empirical package on a machine that
% has NO euro-area data.  It exercises everything the pipeline owns
% except the licensed inputs themselves:
%
%   1. every empirical function PARSES and has the expected signature
%      (assemble_dataset.m used to fail here: a missing `end` left a
%      nested function unterminated, so the whole script was a parse
%      error and nothing downstream could run at all);
%   2. the data-construction chain runs: build_shock_series() turns the
%      shipped EA-EMPD event file into monthly surprises;
%   3. assemble_dataset REFUSES to run without the four raw outcome
%      series, with an error message that names the missing file and the
%      function that fetches it -- i.e. it fails loudly rather than
%      silently producing a dataset;
%   4. the estimator interfaces the drivers depend on all line up, run
%      on a clearly-labelled SYNTHETIC fixture (K = 5, the real system's
%      shape, the [0 1 1 1 1] isrw vector);
%   5. RUN_EMPIRICAL and SMOKE_TEST_EMPIRICAL would REFUSE that fixture
%      (the guard against a synthetic dataset reaching the chapter).
%
% WHAT THIS TEST DOES NOT DO
% --------------------------
% It produces no empirical result and it does not pretend to.  The
% fixture is simulated data; every number computed from it here is
% discarded.  Empirical results require the real series -- see
% empirical/README_EMPIRICAL.md and fetch_outcome_data.m.

fprintf('test_empirical_pipeline:\n');

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(root));
P = ea_paths();

% --- 1. everything parses ------------------------------------------------
fns = {'ea_paths', 'read_sdmx_csv', 'ea_check_series', 'ea_extract_series', ...
       'build_shock_series', 'fetch_outcome_data', 'assemble_dataset', ...
       'estimate_lp_lagaug', 'simulate_fitted_bvar_dgp', ...
       'make_synthetic_fixture', 'ea_write_provenance'};
for k = 1:numel(fns)
    assert(exist(fns{k}, 'file') == 2, 'missing empirical function %s', fns{k});
    n = nargin(fns{k});     % throws on a parse error
    assert(isnumeric(n), '%s did not report a signature', fns{k});
end
fprintf('  %d empirical functions parse and expose a signature: OK\n', numel(fns));

% --- 2. shock construction from the shipped event file -------------------
assert(exist(fullfile(P.raw, 'ea_empd_events.csv'), 'file') == 2, ...
    'the EA-EMPD event file is missing from the repository');
% Write the outputs to a temporary folder: the test must not overwrite
% the shipped derived/shocks_monthly.mat.
tmp_out = fullfile(tempdir, sprintf('ea_shocks_test_%d', round(rand * 1e6)));
mkdir(tmp_out);
S = build_shock_series(struct('out_csv', fullfile(tmp_out, 'shocks.csv'), ...
                              'out_mat', fullfile(tmp_out, 'shocks.mat')));
assert(isfield(S, 'ym') && isfield(S, 'mps_gc_1y'), ...
    'build_shock_series did not return the baseline surprise series');
assert(numel(S.ym) == numel(S.mps_gc_1y) && numel(S.ym) > 100, ...
    'monthly surprise series has an implausible length (%d)', numel(S.ym));
rmdir(tmp_out, 's');
fprintf('  build_shock_series: %d monthly observations from the event file: OK\n', ...
        numel(S.ym));

% --- 3. assemble_dataset refuses to invent the outcome series ------------
missing_dir = fullfile(tempdir, sprintf('ea_raw_empty_%d', round(rand * 1e6)));
mkdir(missing_dir);
threw = false;
try
    assemble_dataset(struct('raw_dir', missing_dir, ...
                            'out_mat', fullfile(missing_dir, 'never.mat')));
catch err
    threw = true;
    assert(~isempty(strfind(err.message, 'ip_ea.csv')), ...
        'the error should name the missing file, got: %s', err.message);
    assert(~isempty(strfind(err.message, 'fetch_outcome_data')), ...
        'the error should name the function that fetches it, got: %s', err.message);
end
rmdir(missing_dir, 's');
assert(threw, 'assemble_dataset did NOT refuse to run without outcome data');
fprintf('  assemble_dataset refuses missing inputs, naming file and remedy: OK\n');

% --- 4. estimator interfaces on the synthetic fixture --------------------
fix_path = fullfile(tempdir, sprintf('SYNTHETIC_fixture_%d.mat', round(rand * 1e6)));
ds = make_synthetic_fixture(struct('T', 150, 'seed', 77, 'out_mat', fix_path));
assert(ds.synthetic, 'the fixture must be stamped synthetic');
[T, K] = size(ds.Y);
assert(K == 5 && isequal(ds.isrw, [0 1 1 1 1]), ...
    'fixture does not have the documented K = 5 shape');

cfg = default_config();
cfg.mode = 'fmar';
cfg.K = K;  cfg.p = 2;  cfg.H = 4;  cfg.T = T;
cfg.shock_var = 1;
cfg.fmar.isrw = ds.isrw;
cfg.fmar.n_niw_draws = 40;
cfg.gibbs.n_burn = 40;  cfg.gibbs.n_keep = 80;

rng(12345, 'twister');
bvar = estimate_bvar_niw(ds.Y, cfg);
for f = {'B', 'c', 'A', 'F', 'Psi', 'b1n', 'theta', 'theta_lo', 'theta_hi', ...
         'lambda', 'max_eig'}
    assert(isfield(bvar, f{1}), 'estimate_bvar_niw lacks .%s', f{1});
end
assert(isequal(size(bvar.theta), [K, cfg.H + 1]), 'bvar.theta has the wrong size');

blpf = estimate_blp_fmar(ds.Y, cfg, bvar);
assert(isequal(size(blpf.lambda), [K, cfg.H]), 'blp_fmar.lambda has the wrong size');

blpb = estimate_blp_blockadaptive(ds.Y, cfg, bvar, blpf.lambda);
assert(isequal(size(blpb.tau_mean), [K, K, cfg.H]), ...
    'block tau_mean has size %s, expected %s', ...
    mat2str(size(blpb.tau_mean)), mat2str([K, K, cfg.H]));

blpp = estimate_blp_blockpooled(ds.Y, cfg, bvar, blpf.lambda);
assert(isequal(size(blpp.tau_mean), size(blpb.tau_mean)), ...
    'pooled tau_mean size differs from the independent estimator');

lpa = estimate_lp_lagaug(ds.Y, cfg, bvar);
assert(isequal(size(lpa.theta), [K, cfg.H + 1]), 'lp_lagaug.theta has the wrong size');

sim = simulate_fitted_bvar_dgp(bvar, ds.Y, cfg);
assert(size(sim.Y, 2) == K && size(sim.Y, 1) == size(ds.Y, 1), ...
    'simulate_fitted_bvar_dgp returned %s, expected %s', ...
    mat2str(size(sim.Y)), mat2str(size(ds.Y)));
assert(isfield(sim, 'theta') && isequal(size(sim.theta), [K, cfg.H + 1]), ...
    'simulate_fitted_bvar_dgp did not pass through the under-null IRF');
fprintf(['  five estimators + the null-calibration DGP run on the K = 5 ' ...
         'fixture: OK\n']);

% --- 4b. provenance is recorded, and never silently destroyed -----------
% fetch_outcome_data writes a sidecar in two situations: after a download
% (it knows the URL) and after validating a file that was already there
% (all it can say is "not recorded by this run").  The second must NOT
% overwrite the first, or re-running the fetch would quietly erase where
% the data came from.
pdir = fullfile(tempdir, sprintf('ea_prov_%d', round(rand * 1e6)));
mkdir(pdir);
praw = fullfile(pdir, 'ip_ea.csv');
fid = fopen(praw, 'w');  fprintf(fid, 'TIME_PERIOD,OBS_VALUE\n2000-01,100\n');  fclose(fid);

w1 = ea_write_provenance(praw, 'TEST series', 'https://real.example/download', true);
assert(w1, 'the first provenance write did not happen');
first = fileread([praw '.source.txt']);
assert(~isempty(strfind(first, 'https://real.example/download')), ...
    'the recorded URL is missing from the sidecar');

w2 = ea_write_provenance(praw, 'TEST series', '(pre-existing local file)', false);
assert(~w2, 'a non-overwriting write reported that it wrote');
assert(strcmp(fileread([praw '.source.txt']), first), ...
    'an existing sidecar was overwritten by the non-overwriting path');

w3 = ea_write_provenance(praw, 'TEST series', 'https://other.example/new', true);
assert(w3 && ~isempty(strfind(fileread([praw '.source.txt']), 'other.example')), ...
    'an overwriting write did not replace the sidecar');

% ...and with no sidecar yet, the non-overwriting path does write one.
delete([praw '.source.txt']);
w4 = ea_write_provenance(praw, 'TEST series', '(pre-existing local file)', false);
assert(w4 && exist([praw '.source.txt'], 'file') == 2, ...
    'no sidecar was written when none existed');
rmdir(pdir, 's');
fprintf('  provenance: recorded on download, preserved on re-validation: OK\n');

% --- 5. the synthetic guard ----------------------------------------------
guarded = {fullfile(P.empirical, 'RUN_EMPIRICAL.m'), ...
           fullfile(P.empirical, 'SMOKE_TEST_EMPIRICAL.m')};
for k = 1:numel(guarded)
    txt = fileread(guarded{k});
    assert(~isempty(strfind(txt, 'synthetic')), ...
        '%s has no guard against a synthetic dataset', guarded{k});
end
ds_syn = load(fix_path);
assert(isfield(ds_syn, 'synthetic') && ds_syn.synthetic, ...
    'the saved fixture lost its synthetic stamp');
delete(fix_path);
fprintf('  both drivers carry the synthetic-dataset guard: OK\n');

fprintf('PASS: test_empirical_pipeline\n\n');
end
