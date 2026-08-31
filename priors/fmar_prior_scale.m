function psi_h = fmar_prior_scale(x, p, h)
% PURPOSE
% -------
% Horizon-specific prior residual scales psi_j(h) for the FMAR Bayesian
% Local Projection.  For each variable j, psi_j(h) is the NEWEY-WEST
% LONG-RUN VARIANCE of the residuals of a UNIVARIATE own-lag local
% projection at horizon h.  This is FMAR's device for making the prior
% scale acknowledge that LP errors are serially correlated (MA(h-1)
% under correct specification): at longer horizons the long-run
% variance grows, which loosens the prior in absolute terms and keeps
% lambda comparable across horizons.  Port of the "use univariate local
% projection to initialize scale (NW corrected)" block of FMAR's
% IRFbayesianLocalProj.m.
%
% MODEL / EQUATIONS
% -----------------
% For variable j, run the univariate LP (data already detrended, but a
% constant is kept, as in FMAR):
%     x_j(t) = a + sum_{l=0}^{p-1} c_l x_j(t-h-l) + u_j(t),
% estimated by OLS.  With demeaned residuals u (FMAR demean them before
% the autocovariances), the Bartlett long-run variance is
%     psi_j(h) = u'u/nT + sum_{l=1}^{L} w_l * ( g_l + g_l' ),
%     g_l = u(l+1:nT)' u(1:nT-l) / (nT - l),
%     w_l = (L + 1 - l) / (L + 1),        L = h + 1.
% NOTE the FMAR scaling quirk, reproduced exactly: the l = 0 term is
% divided by nT while the autocovariances are divided by (nT - l).
% This is harmless (both are consistent long-run variance estimators)
% and is kept so this port matches their code to machine precision.
%
% INPUTS
% ------
% x : (Tx x K) data used for the horizon regressions (detrended).
% p : number of lags in the LP conditioning set.
% h : horizon (>= 1).
%
% OUTPUTS
% -------
% psi_h : (K x 1) long-run residual variances, used as
%         (i)  the diagonal of the IW prior scale at horizon h, and
%         (ii) the relative units of the coefficient prior: the prior
%              variance of a coefficient on (any lag of) variable v in
%              equation i is Sigma_ii * lambda_h^2 / psi_v(h).
%
% DIMENSIONS
% ----------
% The univariate LP at horizon h has nT = Tx - p + 1 - h observations
% and 1 + p regressors (constant + p own lags).
%
% NOTES
% -----
% Under correct specification the residual of an h-step projection is
% MA(h-1), which motivates the truncation lag L = h + 1 used by FMAR.

[Tx, K] = size(x);
nT = Tx - p + 1 - h;
assert(nT > p + 5, 'fmar_prior_scale: too few observations at h = %d.', h);

Zall = build_lp_regressors(x, p);       % ((Tx-p+1) x (1+K*p)); shared builder
psi_h = zeros(K, 1);
L = h + 1;                              % FMAR truncation rule
w = (L + 1 - (1:L)) / (L + 1);          % Bartlett weights

for j = 1:K
    % columns of Zall holding [const, x_j(t), x_j(t-1), ..., x_j(t-p+1)]
    own_cols = [1, 1 + j + (0:p-1) * K];
    Zj = Zall(1:end - h, own_cols);     % (nT x (1+p))
    yj = x(p + h:Tx, j);                % (nT x 1)

    cj = Zj \ yj;                       % OLS
    u  = yj - Zj * cj;
    u  = u - mean(u);                   % FMAR demean the residuals

    G = (u' * u) / nT;                  % l = 0 term (FMAR scaling)
    for l = 1:min(L, nT - 1)
        gl = (u(l+1:nT)' * u(1:nT-l)) / (nT - l);
        G  = G + w(l) * (gl + gl');     % gl is scalar; kept in FMAR form
    end
    psi_h(j) = G;
end

assert(all(psi_h > 0), 'fmar_prior_scale: nonpositive long-run variance.');
end
