function [psi_h, info] = fmar_prior_scale(x, p, h, floor_at_var)
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
% OPTIONAL FLOOR AT THE PLAIN RESIDUAL VARIANCE (floor_at_var = true)
% --------------------------------------------------------------------
% The Bartlett sum adds the sample autocovariances of u with POSITIVE
% weights.  For an overlapping h-step residual of a persistent variable
% those autocovariances are positive and the long-run variance exceeds
% the variance -- the case the correction was designed for.  For a
% variable whose residual is NEGATIVELY autocorrelated (a mean-reverting,
% near-white-noise series such as a monthly high-frequency policy
% surprise) the same sum DRIVES THE ESTIMATE BELOW THE VARIANCE: on the
% euro-area surprise series it is one fifth of the residual variance at
% h >= 20.  Since psi sets the UNITS of the coefficient prior (prior
% variance Sigma_ii lambda^2 / psi_v), an understated psi_v makes the
% prior on every coefficient attached to variable v far too loose, and
% the marginal likelihood then compensates by shrinking the GLOBAL
% lambda_h for every block -- on the euro-area data to ~0.01, with a
% bimodal objective at long horizons (see empirical/README_EMPIRICAL.md).
% With floor_at_var = true,
%     psi_j(h) = max( NW long-run variance, u'u / nT ),
% i.e. a long-run variance estimate is never allowed below the variance
% it is a long-run version of.  The floor does not bind when the residual
% autocorrelations are positive.  On the simulation designs that is the
% case at every h >= 2; at h = 1 the one-step residual is white and its
% sample autocovariances are noise, so the floor binds by a few percent
% in about half the samples, moving lambda_1 by < 0.05 in logs and no
% IRF value by more than 1e-3 (tests/test_psi_floor.m pins these
% numbers).  Default false = the published FMAR scale, which is what
% every stored simulation result uses.
%
% INPUTS
% ------
% x : (Tx x K) data used for the horizon regressions (detrended).
% p : number of lags in the LP conditioning set.
% h : horizon (>= 1).
% floor_at_var : OPTIONAL logical, default false (see above).
%
% OUTPUTS
% -------
% psi_h : (K x 1) long-run residual variances, used as
%         (i)  the diagonal of the IW prior scale at horizon h, and
%         (ii) the relative units of the coefficient prior: the prior
%              variance of a coefficient on (any lag of) variable v in
%              equation i is Sigma_ii * lambda_h^2 / psi_v(h).
% info  : struct with .psi_nw (the unfloored NW estimate), .psi_var (the
%         plain residual variance u'u/nT) and .floor_bound (K x 1
%         logical, true where the floor was applied).
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

if nargin < 4 || isempty(floor_at_var), floor_at_var = false; end

[Tx, K] = size(x);
nT = Tx - p + 1 - h;
assert(nT > p + 5, 'fmar_prior_scale: too few observations at h = %d.', h);

Zall = build_lp_regressors(x, p);       % ((Tx-p+1) x (1+K*p)); shared builder
psi_nw  = zeros(K, 1);
psi_var = zeros(K, 1);
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
    psi_var(j) = G;
    for l = 1:min(L, nT - 1)
        gl = (u(l+1:nT)' * u(1:nT-l)) / (nT - l);
        G  = G + w(l) * (gl + gl');     % gl is scalar; kept in FMAR form
    end
    psi_nw(j) = G;
end

if floor_at_var
    psi_h = max(psi_nw, psi_var);
else
    psi_h = psi_nw;
end
info.psi_nw = psi_nw;
info.psi_var = psi_var;
info.floor_bound = floor_at_var & (psi_var > psi_nw);

assert(all(psi_h > 0), 'fmar_prior_scale: nonpositive long-run variance.');
end
