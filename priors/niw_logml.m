function [logML, betahat, sigmahat] = niw_logml(Y, X, b, omega, psi, d)
% PURPOSE
% -------
% Closed-form LOG MARGINAL LIKELIHOOD of a conjugate Normal - Inverse-
% Wishart (NIW) multivariate regression, plus the posterior mode of the
% coefficients and of the residual covariance.  This is a line-by-line
% port of the logML computation in the FMAR replication files
% logMLVAR_formin.m / logMLBLP_formin.m (which themselves follow
% Giannone, Lenza & Primiceri, 2015), specialised to the switches used
% by the FMAR IRF application: no sum-of-coefficients dummy (noc = 0),
% no single-unit-root dummy (sur = 0), psi fixed (MNpsi = 0).
%
% MODEL / EQUATIONS
% -----------------
% System of n equations with common regressors (T observations):
%     Y = X B + E,   E_t ~ N(0, Sigma),
%     Sigma ~ IW(diag(psi), d),
%     vec(B) | Sigma ~ N( vec(b), Sigma (x) Omega ),  Omega = diag(omega).
% Marginally, Y follows a matrix-variate Student-t whose density has
% the closed form (GLP 2015, eq. (A.2) of their appendix; identical to
% the FMAR code):
%
%   logML = -nT/2 log(pi)
%           + sum_{i=0}^{n-1} [ lnGamma((T+d-i)/2) - lnGamma((d-i)/2) ]
%           - T/2 * sum(log psi)
%           - n/2 * log det( I + Omega^{1/2} X'X Omega^{1/2} )
%           - (T+d)/2 * log det( I + Psi^{-1/2} R Psi^{-1/2} ),
%
% where R = Ehat'Ehat + (Bhat - b)' Omega^{-1} (Bhat - b) and
%     Bhat = (X'X + Omega^{-1})^{-1} (X'Y + Omega^{-1} b)
% is the posterior mode/mean of the coefficients (independent of Sigma
% thanks to the Kronecker prior).  Both log-determinants are computed,
% exactly as in the FMAR code, through the eigenvalues of the
% symmetrised inner matrices with tiny negative eigenvalues floored at
% zero before adding 1 (a numerical guard against round-off).
%
% The posterior mode of the residual covariance is
%     sigmahat = ( Ehat'Ehat + diag(psi) + (Bhat-b)' Omega^{-1} (Bhat-b) )
%                / (T + d + n + 1),
% i.e. the mode of IW(S_end, d + T) with
%     S_end = Ehat'Ehat + diag(psi) + (Bhat-b)' Omega^{-1} (Bhat-b).
%
% INPUTS
% ------
% Y     : (T x n) dependent variables.
% X     : (T x k) common regressors (constant in column 1).
% b     : (k x n) prior mean of the coefficients.
% omega : (k x 1) diagonal of Omega (prior coefficient covariance,
%         conditional on Sigma; units of "per Sigma_ii").
% psi   : (n x 1) diagonal of the IW prior scale.
% d     : scalar IW prior degrees of freedom (FMAR: d = n + 2).
%
% OUTPUTS
% -------
% logML    : scalar log marginal likelihood (WITHOUT any hyperprior
%            term; the caller adds log p(lambda) etc. if wanted).
% betahat  : (k x n) posterior mode/mean of the coefficients.
% sigmahat : (n x n) posterior mode of the residual covariance.
%
% DIMENSIONS
% ----------
% T observations, n equations, k regressors.
%
% NOTES
% -----
% Verified against the original FMAR logMLBLP_formin.m in
% tests/test_fmar_port.m (equality up to the hyperprior term).

[T, n] = size(Y);
k = size(X, 2);
assert(numel(omega) == k && numel(b) == k * n && numel(psi) == n, ...
    'niw_logml: dimension mismatch.');
omega = omega(:);
psi   = psi(:);
assert(all(omega > 0) && all(psi > 0), 'niw_logml: omega, psi must be > 0.');

% Posterior mode of coefficients (stable solve, no explicit inverse).
betahat = (X' * X + diag(1 ./ omega)) \ (X' * Y + diag(1 ./ omega) * b);

Ehat = Y - X * betahat;

% --- the two log-determinant terms, FMAR's eigenvalue route -----------
aaa = diag(sqrt(omega)) * (X' * X) * diag(sqrt(omega));
bbb = diag(1 ./ sqrt(psi)) * ...
      (Ehat' * Ehat + (betahat - b)' * diag(1 ./ omega) * (betahat - b)) * ...
      diag(1 ./ sqrt(psi));

eigaaa = real(eig(aaa));  eigaaa(eigaaa < 1e-12) = 0;  eigaaa = eigaaa + 1;
eigbbb = real(eig(bbb));  eigbbb(eigbbb < 1e-12) = 0;  eigbbb = eigbbb + 1;

logML = - n * T * log(pi) / 2 ...
        + sum(gammaln((T + d - (0:n-1)) / 2) - gammaln((d - (0:n-1)) / 2)) ...
        - T * sum(log(psi)) / 2 ...
        - n * sum(log(eigaaa)) / 2 ...
        - (T + d) * sum(log(eigbbb)) / 2;

sigmahat = (Ehat' * Ehat + diag(psi) + ...
            (betahat - b)' * diag(1 ./ omega) * (betahat - b)) / (T + d + n + 1);

assert(isfinite(logML), 'niw_logml: non-finite log marginal likelihood.');
end
