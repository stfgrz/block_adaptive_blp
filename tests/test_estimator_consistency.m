function test_estimator_consistency()
% PURPOSE
% -------
% Close the one gap the other tests leave open: they verify the TRUTH
% against itself (analytic vs simulated) and the ESTIMATORS against
% each other (nesting), but never the DATA against the TRUTH.  A
% mismatch between the simulation engine and the parameters used for
% the analytic IRF (e.g. a lag-ordering or shock-normalisation bug)
% would slip through all of them.
%
% The check: on a VERY LARGE sample from the CORRECTLY specified DGP,
% the VAR and LP estimators are consistent, so their IRFs must
% converge to the analytic truth.  Failure = systematic bug; success =
% any single-dataset discrepancy at T = 200 is sampling error (which
% includes the well-known downward small-sample bias of estimated
% persistence).
%
% MODEL / EQUATIONS
% -----------------
% With T = 20000 the OLS standard error of a VAR coefficient is about
% 1/sqrt(T) ~ 0.007; compounding over 20 horizons keeps the IRF error
% well below the 0.05 tolerance used here.  LP is noisier at long
% horizons, so its tolerance is looser (0.08).  The impact vector b1n
% must also match the true normalised impact column (tolerance 0.05).
%
% As a bonus, the same large-T logic quantifies the POPULATION prior
% bias of the sparse DGP: the gap between the best-fitting VAR(2) IRF
% (= the BLP prior centre) and the true VAR(3) IRF.  This gap is what
% the block-adaptive mechanism exists to fix; the test only PRINTS it
% (no assertion), so the researcher can judge whether the calibration
% makes the experiment winnable.
%
% INPUTS / OUTPUTS
% ----------------
% None; prints PASS or raises an assertion error.  Runtime: a few
% seconds in MATLAB, up to ~1 minute in Octave.

fprintf('--- test_estimator_consistency ---\n');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'config'), fullfile(root, 'dgp'), ...
        fullfile(root, 'estimators'), fullfile(root, 'priors'), ...
        fullfile(root, 'samplers'), fullfile(root, 'utils'));

cfg = default_config();
cfg.T = 20000;                          % large sample: consistency regime
cfg.var.n_ci_sim = 50;                  % CIs irrelevant here, keep cheap
rng(20260820, 'twister');

% ---- Correct DGP: estimators must converge to the analytic truth ----
dgp = simulate_var_dgp(cfg);
v   = estimate_var(dgp.Y, cfg);

err_b1 = max(abs(v.b1n - dgp.theta_true(:, 1)));
assert(err_b1 < 0.05, 'b1n far from true impact column: %.3f', err_b1);
fprintf('  impact vector b1n vs truth: max err %.4f  OK\n', err_b1);

err_var = max(abs(v.theta(:) - dgp.theta_true(:)));
assert(err_var < 0.05, ...
    ['VAR IRF does not converge to the analytic truth (max gap %.3f): ' ...
     'systematic data/truth mismatch.'], err_var);
fprintf('  VAR IRF vs truth (T = %d): max err %.4f  OK\n', cfg.T, err_var);

lp = estimate_lp(dgp.Y, cfg, v.b1n);
err_lp = max(abs(lp.theta(:) - dgp.theta_true(:)));
assert(err_lp < 0.08, ...
    'LP IRF does not converge to the analytic truth (max gap %.3f).', err_lp);
fprintf('  LP  IRF vs truth (T = %d): max err %.4f  OK\n', cfg.T, err_lp);

% ---- Sparse DGP: PRINT the population prior bias (no assertion) ------
dgp2 = simulate_sparse_misspec_dgp(cfg);
v2   = estimate_var(dgp2.Y, cfg);       % best-fitting VAR(2) = prior centre
gap  = abs(v2.theta - dgp2.theta_true); % (K x (H+1)) prior bias by horizon
fprintf(['  [info] sparse DGP, population VAR(2)-prior bias vs true IRF:\n' ...
         '         max over horizons, per variable: %.4f  %.4f  %.4f\n' ...
         '         (this is the signal the block-adaptive prior must find;\n' ...
         '          compare it to single-sample noise of roughly 0.05-0.10\n' ...
         '          at T = 200 before interpreting detection results)\n'], ...
        max(gap(1, :)), max(gap(2, :)), max(gap(3, :)));

fprintf('PASS: test_estimator_consistency\n\n');
end
