function Mu = var_implied_lp_prior(var_est, h)
% PURPOSE
% -------
% Map an estimated VAR(p) into the LOCAL-PROJECTION COEFFICIENTS it
% implies at horizon h.  This matrix is the PRIOR CENTRE mu_h^VAR of
% both Bayesian LP estimators.  The mapping is exact linear algebra,
% not an approximation, and is verified numerically against iterated
% VAR forecasts in tests/test_var_prior_mapping.m.
%
% MODEL / EQUATIONS
% -----------------
% Stack the state s(t) = [y(t); y(t-1); ...; y(t-p+1)]  (Kp x 1) and
% write the VAR(p) in companion form
%     s(t+1) = C + F s(t) + [e(t+1); 0],
% with C = [c; 0; ...; 0] and F the (Kp x Kp) companion matrix.
% Iterating h steps and taking conditional expectations,
%
%     E[ y(t+h) | s(t) ] = d_h + (J F^h) s(t),          J = [I_K, 0],
%     d_h = ( sum_{j=0}^{h-1} Psi_j ) c,   Psi_j = J F^j J'.
%
% Because the conditional expectation is EXACTLY LINEAR in s(t), the
% population regression of y(t+h) on z(t) = [1; s(t)] under the VAR has
% coefficient matrix
%
%     Mu = [ d_h' ;  (J F^h)' ]        (m x K),  m = 1 + K p,
%
% i.e. column i of Mu holds the LP coefficients of equation i implied
% by the VAR.  Two useful identities used by the tests:
%   * the coefficient block on y(t) (rows 2 : 1+K of Mu, transposed)
%     equals the reduced-form IRF matrix Psi_h;
%   * Z * Mu reproduces the h-step-ahead iterated VAR forecasts from
%     every state in the sample, to machine precision.
%
% *** Prototype caveat *** ------------------------------------------
% This is the internally consistent conditional-mean mapping.  The
% published BLP of Ferreira, Miranda-Agrippino and Ricco specifies its
% prior moments (centre AND scale) through its own machinery; when
% their replication code is adopted, THIS FUNCTION is the substitution
% point: replace its output while keeping the (m x K) interface.
%
% INPUTS
% ------
% var_est : output struct of estimate_var (uses .A, .c, .F, .Psi).
% h       : horizon, integer >= 0.
%
% OUTPUTS
% -------
% Mu : (m x K) prior-centre matrix; Mu(:, i) is mu_{h}^VAR for LP
%      equation i, ordered like z(t) = [1; y(t); ...; y(t-p+1)].
%
% DIMENSIONS
% ----------
% K = size(var_est.A, 1); p = size(var_est.A, 3); m = 1 + K*p.
%
% NOTES
% -----
% h = 0 returns the trivial identity mapping ([0; I; 0...]); it is not
% used in estimation (the h = 0 step is the identification step).

K = size(var_est.A, 1);
p = size(var_est.A, 3);
m = 1 + K * p;
assert(h >= 0 && h == round(h), 'var_implied_lp_prior: h must be >= 0.');

J  = [eye(K), zeros(K, K * (p - 1))];
Fh = var_est.F^h;                       % (Kp x Kp); p is small, this is fine

% Intercept: d_h = (sum_{j=0}^{h-1} Psi_j) c.  Psi_j is recomputed here
% from F powers (independent of var_est.Psi) to keep this file
% self-contained; the tests check consistency with var_est.Psi.
d = zeros(K, 1);
Fj = eye(K * p);
for j = 0:h-1
    d  = d + (J * Fj * J') * var_est.c;
    Fj = Fj * var_est.F;
end

Mu = zeros(m, K);
Mu(1, :)     = d';                      % intercept row
Mu(2:end, :) = (J * Fh)';               % (Kp x K): lag-coefficient rows

assert(all(isfinite(Mu(:))), 'var_implied_lp_prior: non-finite centre.');
end
