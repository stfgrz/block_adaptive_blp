function dgp = simulate_var_dgp(cfg)
% PURPOSE
% -------
% DGP 1: simulate data from a CORRECTLY SPECIFIED, stable VAR(2).
% Under this DGP the estimated VAR(p) prior centre is (asymptotically)
% correct for every block, so the block-adaptive prior should learn
% small local scales tau_g and no estimator should be badly biased.
%
% MODEL / EQUATIONS
% -----------------
%   y(t) = c + A_1 y(t-1) + A_2 y(t-2) + e(t),  e(t) = B0 eta(t),
%   eta(t) ~ N(0, I_K),  Sigma = B0 B0' (B0 lower triangular).
% Coefficients are chosen for moderate, economically plausible
% persistence with the companion matrix comfortably inside the unit
% circle (max |eigenvalue| checked by assertion below).
%
% INPUTS
% ------
% cfg : configuration struct (uses cfg.K, cfg.T, cfg.burnin, cfg.H,
%       cfg.shock_var).  The RNG state is NOT set here: the CALLER
%       controls the seed so Monte Carlo replications are reproducible.
%
% OUTPUTS
% -------
% dgp : struct with fields
%   .name          : 'correct'
%   .type          : 'var'
%   .A             : (K x K x 2) true AR matrices
%   .M             : [] (no MA part)
%   .c             : (K x 1) intercept
%   .Sigma, .B0    : (K x K) innovation covariance and its Cholesky
%   .Y             : (T x K) simulated data (burn-in discarded)
%   .eta           : (K x T) structural shocks kept with the sample
%   .theta_true    : (K x (H+1)) true IRF to the unit shock (from
%                    compute_true_irf, analytic; NEVER from an estimator)
%   .misspec_block : [] (nothing is misspecified)
%   .description   : plain-language description
%
% DIMENSIONS
% ----------
% K = 3 in the default configuration; the matrices below are written
% for K = 3 and an assertion enforces that.
%
% NOTES
% -----
% This DGP is the benchmark "the VAR prior is right" case.

assert(cfg.K == 3, ...
    'simulate_var_dgp: default true parameters are written for K = 3.');

% --- True parameters -------------------------------------------------
A1 = [0.55 0.10 0.05;
      0.10 0.50 0.10;
      0.05 0.10 0.45];
A2 = [0.15 0.00 0.00;
      0.00 0.10 0.05;
      0.00 0.05 0.15];
A          = cat(3, A1, A2);
c          = [0.20; 0.10; 0.10];
B0         = [1.00 0.00 0.00;
              0.30 0.90 0.00;
              0.15 0.25 0.80];
Sigma      = B0 * B0';

% Stability check of the TRUE process (companion eigenvalues).
F = companion_from_A(A);
max_eig = max(abs(eig(F)));
assert(max_eig < 0.95, ...
    'simulate_var_dgp: true VAR not comfortably stable (max|eig|=%.3f).', max_eig);

% --- Assemble DGP struct and simulate --------------------------------
dgp.name          = 'correct';
dgp.type          = 'var';
dgp.A             = A;
dgp.M             = [];
dgp.c             = c;
dgp.Sigma         = Sigma;
dgp.B0            = B0;
dgp.misspec_block = [];
dgp.description   = 'Correctly specified stable VAR(2).';

[dgp.Y, dgp.eta] = simulate_linear_dgp(dgp, cfg);
dgp.theta_true   = compute_true_irf(dgp, cfg, 'analytic');
end

% =====================================================================
function F = companion_from_A(A)
% Companion matrix of a VAR with AR matrices A (K x K x p).
[K, ~, p] = size(A);
F = zeros(K * p);
for j = 1:p
    F(1:K, (j-1)*K + (1:K)) = A(:, :, j);
end
if p > 1
    F(K+1:end, 1:K*(p-1)) = eye(K * (p - 1));
end
end
