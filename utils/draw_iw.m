function X = draw_iw(S, df)
% PURPOSE
% -------
% Draw one matrix from the Inverse-Wishart distribution IW(S, df) using
% only base MATLAB (rand/randn via utils/draw_gamma.m).  Needed for the
% FMAR-style Bayesian VAR posterior, whose residual covariance is
% IW-distributed; MATLAB's iwishrnd lives in a toolbox, so this project
% provides its own transparent implementation.
%
% MODEL / EQUATIONS
% -----------------
% Parameterisation (the one used by MATLAB's iwishrnd and by the FMAR
% code): X ~ IW(S, df) has density
%     p(X) ~ |X|^{-(df+n+1)/2} exp( -tr(S X^{-1}) / 2 ),
% with E[X] = S / (df - n - 1) for df > n + 1.
%
% Sampling route:
%     X ~ IW(S, df)   <=>   X^{-1} ~ Wishart(S^{-1}, df),
% and Wishart(V, df) is drawn by the Bartlett decomposition:
%     V = L L' (lower Cholesky),
%     A lower triangular with A(i,i) = sqrt(chi2(df - i + 1)) and
%     A(i,j) ~ N(0,1) for i > j,
%     W = (L A)(L A)'  ~  Wishart(V, df).
% chi-square draws come from draw_gamma: chi2(k) = Gamma(k/2, scale 2).
% Finally X = W^{-1}, computed via a triangular solve on (L A).
%
% INPUTS
% ------
% S  : (n x n) symmetric positive definite scale matrix.
% df : scalar degrees of freedom, df > n - 1 (integer not required).
%
% OUTPUTS
% -------
% X : (n x n) symmetric positive definite IW(S, df) draw.
%
% DIMENSIONS
% ----------
% n = size(S, 1).
%
% NOTES
% -----
% Fully driven by rand/randn, hence reproducible under rng(seed).
% Inversion of W uses the triangular factor directly:
%     W = (L A)(L A)'  =>  W^{-1} = (L A)^{-T} (L A)^{-1},
% i.e. X = M' * M with M = inv(L A) obtained from a triangular solve,
% which is more stable than inv() on the assembled W.

n = size(S, 1);
assert(size(S, 2) == n, 'draw_iw: S must be square.');
assert(df > n - 1, 'draw_iw: need df > n - 1.');

% Lower Cholesky of the WISHART scale V = S^{-1}: if S = Ls Ls' then
% S^{-1} = Ls'^{-1} Ls^{-1}, whose lower Cholesky is obtained from the
% inverse of the upper factor.  Compute via triangular solve for
% stability.
Ls = safe_chol_lower(S);               % S = Ls Ls'
Linv_upper = Ls' \ eye(n);             % Ls'^{-1}, upper triangular
% V = S^{-1} = Linv_upper * Linv_upper'; its lower Cholesky is
% NOT Linv_upper (which is upper), so factorise V explicitly:
V  = Linv_upper * Linv_upper';
Lv = safe_chol_lower((V + V') / 2);    % V = Lv Lv'

% Bartlett factor A (lower triangular).
A = zeros(n);
for i = 1:n
    A(i, i) = sqrt(draw_gamma((df - i + 1) / 2, 2));   % chi2(df-i+1)
    for j = 1:i-1
        A(i, j) = randn();
    end
end

M = Lv * A;                            % W = M M' ~ Wishart(V, df)
Minv = M \ eye(n);                     % triangular-ish solve (M lower tri)
X = Minv' * Minv;                      % X = W^{-1} ~ IW(S, df)
X = (X + X') / 2;                      % enforce exact symmetry

assert(all(isfinite(X(:))), 'draw_iw: non-finite draw.');
end
