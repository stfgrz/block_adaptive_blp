function L = safe_chol_lower(S)
% PURPOSE
% -------
% Numerically guarded lower-triangular Cholesky factor of a covariance /
% precision matrix.  Used everywhere a Cholesky is taken so that tiny
% asymmetries or near-singularity from floating-point arithmetic do not
% crash a long Monte Carlo run.
%
% MODEL / EQUATIONS
% -----------------
% Returns L with L*L' = S (up to the added jitter).  Stabilisation:
%   1. Symmetrise: S <- (S + S')/2. Removes O(eps) asymmetry only.
%   2. If chol fails, add jitter eps_j * mean(diag(S)) * I with
%      eps_j = 1e-12, 1e-10, 1e-8, 1e-6 in turn.  Jitter inflates all
%      variances by a relatively negligible amount; if even 1e-6 fails,
%      an error is raised rather than silently returning garbage.
%
% INPUTS
% ------
% S : (n x n) symmetric positive (semi-)definite matrix.
%
% OUTPUTS
% -------
% L : (n x n) lower triangular with L*L' ~= S.
%
% DIMENSIONS
% ----------
% n = size(S,1).
%
% NOTES
% -----
% Any jitter beyond the first attempt is reported with a warning so
% numerical trouble is visible, never hidden.

assert(size(S,1) == size(S,2), 'safe_chol_lower: S must be square.');
assert(all(isfinite(S(:))), 'safe_chol_lower: S has non-finite entries.');

S = (S + S') / 2;                       % remove floating-point asymmetry
[L, flag] = chol(S, 'lower');
if flag == 0
    return;
end

scale = mean(abs(diag(S))) + eps;
for ej = [1e-12, 1e-10, 1e-8, 1e-6]
    [L, flag] = chol(S + ej * scale * eye(size(S,1)), 'lower');
    if flag == 0
        warning('safe_chol_lower:jitter', ...
            'Cholesky needed jitter %g * mean(diag).', ej);
        return;
    end
end
error('safe_chol_lower: matrix not positive definite even with jitter.');
end
