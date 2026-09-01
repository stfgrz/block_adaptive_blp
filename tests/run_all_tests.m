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

fprintf('==============================================\n');
fprintf(' ALL TESTS PASSED\n');
fprintf('==============================================\n');
end
