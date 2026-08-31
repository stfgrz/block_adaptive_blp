function test_sampler_sanity()
% PURPOSE
% -------
% Exercise samplers/gibbs_block_horseshoe.m on a SIMPLE Gaussian
% regression whose posterior behaviour is easy to reason about, and
% assert (not merely eyeball) the basics:
%   1. output dimensions;
%   2. all draws finite;
%   3. positive variances (sigma2, tau2);
%   4. exact reproducibility given the same seed;
%   5. posterior mean of beta close to OLS/truth when the data are
%      informative and the prior is loose;
%   6. shrinkage direction: with a very TIGHT prior centred at mu0 the
%      posterior mean must move towards mu0 relative to OLS;
%   7. no numerical explosions (bounded draws).
%
% MODEL / EQUATIONS
% -----------------
% y = Z beta0 + e, e ~ N(0, 0.5^2), N = 400, m = 5 (intercept + 2
% blocks of 2).  With d = 1, lambda = 2, half-Cauchy tau, the prior is
% loose, so the posterior should sit close to OLS: |post - OLS| within
% a few posterior standard deviations, and within 0.15 of the truth.
%
% INPUTS / OUTPUTS
% ----------------
% None; prints PASS or raises an assertion error.

fprintf('--- test_sampler_sanity ---\n');
rng(4242, 'twister');

N = 400;  m = 5;
beta0 = [0.5; 1.0; -0.8; 0.0; 0.6];
Z = [ones(N, 1), randn(N, m - 1)];
y = Z * beta0 + 0.5 * randn(N, 1);

prior.mu       = zeros(m, 1);
prior.d        = [100; ones(m - 1, 1)];   % loose intercept, unit blocks
prior.block_id = [0; 1; 1; 2; 2];
prior.lambda   = 2.0;                     % loose global scale
prior.a0 = 0.01;  prior.b0 = 0.01;

opts = struct('n_burn', 500, 'n_keep', 3000, 'sample_tau', true, ...
              'tau2_min', 1e-10, 'tau2_max', 1e8, 'seed', 99);

out1 = gibbs_block_horseshoe(y, Z, prior, opts);
out2 = gibbs_block_horseshoe(y, Z, prior, opts);   % same seed

% 1. Dimensions.
assert(isequal(size(out1.beta_draws), [3000, m]), 'beta_draws dims.');
assert(isequal(size(out1.sig2_draws), [3000, 1]), 'sig2_draws dims.');
assert(isequal(size(out1.tau2_draws), [3000, 2]), 'tau2_draws dims.');
fprintf('  dimensions OK\n');

% 2. Finite values.
assert(out1.diag.all_finite, 'non-finite draws.');
fprintf('  finiteness OK\n');

% 3. Positive variances.
assert(all(out1.sig2_draws > 0), 'sigma2 draws not positive.');
assert(all(out1.tau2_draws(:) > 0), 'tau2 draws not positive.');
fprintf('  positive variances OK\n');

% 4. Reproducibility.
assert(isequal(out1.beta_draws, out2.beta_draws) && ...
       isequal(out1.tau2_draws, out2.tau2_draws), 'not reproducible.');
fprintf('  reproducibility OK\n');

% 5. Posterior mean near OLS and truth.
beta_ols = Z \ y;
post_sd = std(out1.beta_draws, 0, 1)';
gap_ols = abs(out1.beta_mean - beta_ols);
assert(all(gap_ols < 4 * post_sd + 0.02), ...
    'posterior mean far from OLS under a loose prior.');
assert(max(abs(out1.beta_mean - beta0)) < 0.15, ...
    'posterior mean far from truth.');
fprintf('  posterior location OK (max |post-OLS| = %.3f)\n', max(gap_ols));

% 6. Shrinkage direction with a tight prior towards mu0 = 2.
prior_tight = prior;
prior_tight.mu = 2 * ones(m, 1);
prior_tight.lambda = 0.01;              % very tight around mu0
opts_fix = opts;  opts_fix.sample_tau = false;  opts_fix.tau_fixed = 1;
out3 = gibbs_block_horseshoe(y, Z, prior_tight, opts_fix);
moved = abs(out3.beta_mean(2:end) - 2) < abs(beta_ols(2:end) - 2) - 1e-6;
assert(all(moved), 'tight prior did not pull the posterior towards mu0.');
% sigma2 must blow up to absorb the now-unexplained variation:
assert(mean(out3.sig2_draws) > mean(out1.sig2_draws), ...
    'sigma2 did not rise under the badly centred tight prior.');
fprintf('  shrinkage direction OK\n');

% 7. No explosions.
assert(max(abs(out1.beta_draws(:))) < 10, 'beta draws exploded.');
assert(max(out1.sig2_draws) < 100, 'sigma2 draws exploded.');
fprintf('PASS: test_sampler_sanity\n\n');
end
