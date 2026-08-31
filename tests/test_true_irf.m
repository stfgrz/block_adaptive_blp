function test_true_irf()
% PURPOSE
% -------
% Cross-validate dgp/compute_true_irf.m: for each DGP, the ANALYTIC
% true IRF (MA recursion) must coincide with the SIMULATION-BASED
% true IRF (shocked minus unshocked paths with common random numbers).
% For linear DGPs each simulated pair difference is exact, so the two
% methods must agree to near machine precision -- a strong joint check
% on the MA recursion, the companion algebra and the path simulator.
%
% Also checks that the baseline VAR(2) parameters duplicated inside
% simulate_sparse_misspec_dgp.m agree with simulate_var_dgp.m (lags 1
% and 2 of the sparse DGP must match the correct DGP exactly).
%
% INPUTS / OUTPUTS
% ----------------
% None; prints PASS or raises an assertion error.

fprintf('--- test_true_irf ---\n');
cfg = default_config();
cfg.true_irf.n_paths = 50;              % linear DGP: any count is exact

names = {'correct', 'sparse', 'dense'};
sims  = {@simulate_var_dgp, @simulate_sparse_misspec_dgp, ...
         @simulate_dense_misspec_dgp};

dgps = cell(1, 3);
for k = 1:3
    rng(cfg.seed + k, 'twister');
    dgp = sims{k}(cfg);
    dgps{k} = dgp;
    th_ana = compute_true_irf(dgp, cfg, 'analytic');
    rng(1, 'twister');                  % seed for the simulation method
    th_sim = compute_true_irf(dgp, cfg, 'simulation');
    err = max(abs(th_ana(:) - th_sim(:)));
    fprintf('  DGP %-8s analytic vs simulated: max err %.2e\n', names{k}, err);
    assert(err < 1e-8, 'true-IRF mismatch for DGP %s: %.2e', names{k}, err);
    assert(isequal(size(th_ana), [cfg.K, cfg.H + 1]), 'true-IRF dims.');
end

% Baseline parameters shared between DGP 1 and DGP 2.
errA = max(abs(dgps{1}.A(:) - reshape(dgps{2}.A(:, :, 1:2), [], 1)));
errc = max(abs(dgps{1}.c - dgps{2}.c));
errS = max(abs(dgps{1}.Sigma(:) - dgps{2}.Sigma(:)));
assert(max([errA, errc, errS]) < 1e-12, ...
    'DGP 1 / DGP 2 baseline parameter copies disagree.');
fprintf('  DGP1/DGP2 shared baseline parameters agree\n');
fprintf('PASS: test_true_irf\n\n');
end
