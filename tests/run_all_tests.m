function run_all_tests()
% PURPOSE
% -------
% Run every test in tests/ in a fixed order.  Any failure raises an
% assertion error immediately (nothing is caught or hidden).
%
% INPUTS / OUTPUTS
% ----------------
% None; prints one PASS line per test and a final summary.

% Add project folders (absolute paths, so this works from any directory).
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'config'), fullfile(root, 'dgp'), ...
        fullfile(root, 'estimators'), fullfile(root, 'priors'), ...
        fullfile(root, 'samplers'), fullfile(root, 'montecarlo'), ...
        fullfile(root, 'plots'), fullfile(root, 'tests'), ...
        fullfile(root, 'utils'));
addpath(genpath(fullfile(root, 'empirical')));   % Ch. 7 package + its test

fprintf('==============================================\n');
fprintf(' Running all tests (block_adaptive_blp)\n');
fprintf('==============================================\n\n');

test_var_prior_mapping();
test_true_irf();
test_estimator_consistency();
test_sampler_sanity();
test_nesting();
test_isrw_vector();    % vectorised cfg.fmar.isrw prior mean (Ch. 7 needs it)
test_fmar_port();      % vs original FMAR code (skips if not installed)
test_fmar_nesting();   % block-adaptive(tau=1) == FMAR posterior mean
                       % (sampling-based, kept as the original check)

% --- tests added with the benchmark-fairness / pooling / diagnostics work
test_fmar_nesting_exact();  % tolerance-free nesting, 4 DGPs x 2 seeds x
                            % every equation and horizon, both h1 modes,
                            % both adaptive estimators
test_pooled_estimator();    % the horizon-pooling prior and its two limits
test_mc_metrics();          % metric definitions, both IRMSE conventions,
                            % reproducibility record and csv export
test_tau_diagnostics();     % localisation / false-positive accounting
test_chunk_merge();         % parallel chunking is exact; cfg round trip
test_approximations();      % the claims docs/APPROXIMATIONS.md rests on
test_psi_floor();           % the optional prior-scale floor: inert on the
                            % simulations, active + nesting-exact where
                            % the euro-area data need it
test_fixed_tau_blocks();    % holding chosen blocks at tau = 1 (the
                            % instrument block in the empirical chapter)
test_empirical_pipeline();  % Ch. 7 package on a synthetic fixture

fprintf('==============================================\n');
fprintf(' ALL TESTS PASSED\n');
fprintf('==============================================\n');
end
