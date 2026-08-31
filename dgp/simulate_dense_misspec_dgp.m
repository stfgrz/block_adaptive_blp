function dgp = simulate_dense_misspec_dgp(cfg)
% PURPOSE
% -------
% DGP 3: DENSE misspecification of the fitted VAR(2).  The truth is a
% VARMA(2,1): the same baseline VAR(2) plus a moving-average term that
% loads on ALL equations,
%
%     y(t) = c + A_1 y(t-1) + A_2 y(t-2) + e(t) + M e(t-1),
%
% with M having sizeable entries in every row and column.
%
% WHY THE FITTED VAR(2) IS MISSPECIFIED (EVERYWHERE)
% --------------------------------------------------
% A VARMA process has an infinite-order VAR representation
%     y(t) = c* + sum_{j>=1} Pi_j y(t-j) + e(t),
% with Pi_j decaying geometrically but nonzero at ALL lags and in ALL
% equations because M is dense.  Truncating at p = 2 lags therefore
% leaves omitted dynamics in EVERY equation and EVERY variable block:
% no single block is "the" wrong one (dgp.misspec_block = []).  This is
% the polar case to DGP 2: block-adaptive escape has no sparse target,
% so any advantage over the global prior should shrink or vanish.
%
% MODEL / EQUATIONS
% -----------------
% A_1, A_2, c, Sigma as in DGP 1;  M given below (spectral radius well
% below 1, so the MA part is invertible and the implicit Pi_j decay).
% True IRFs are analytic:
%   Psi_0 = I,  Psi_1 = A_1 + M,
%   Psi_h = A_1 Psi_{h-1} + A_2 Psi_{h-2},  h >= 2,
% handled by compute_true_irf (type 'varma').
%
% INPUTS
% ------
% cfg : configuration struct; RNG state controlled by the caller.
%
% OUTPUTS
% -------
% dgp : struct with the same fields as the other DGP files and
%       .misspec_block = [] (dense: all blocks affected).
%
% DIMENSIONS
% ----------
% Written for K = 3 (asserted).
%
% NOTES
% -----
% Stationarity depends only on the AR part, which is unchanged and
% stable; invertibility of the MA part is asserted so the VAR(inf)
% representation exists.

assert(cfg.K == 3, ...
    'simulate_dense_misspec_dgp: true parameters written for K = 3.');

A1 = [0.55 0.10 0.05;
      0.10 0.50 0.10;
      0.05 0.10 0.45];
A2 = [0.15 0.00 0.00;
      0.00 0.10 0.05;
      0.00 0.05 0.15];
c  = [0.20; 0.10; 0.10];
B0 = [1.00 0.00 0.00;
      0.30 0.90 0.00;
      0.15 0.25 0.80];

M  = [0.40 0.15 0.10;
      0.15 0.35 0.15;
      0.10 0.15 0.40];

assert(max(abs(eig(M))) < 1, ...
    'simulate_dense_misspec_dgp: MA polynomial not invertible.');

dgp.name          = 'dense';
dgp.type          = 'varma';
dgp.A             = cat(3, A1, A2);
dgp.M             = M;
dgp.c             = c;
dgp.Sigma         = B0 * B0';
dgp.B0            = B0;
dgp.misspec_block = [];
dgp.description   = ['True VARMA(2,1): baseline VAR(2) plus a dense ' ...
    'MA(1) term in all equations. The fitted VAR(2) omits dynamics ' ...
    'in every block (dense misspecification).'];

[dgp.Y, dgp.eta] = simulate_linear_dgp(dgp, cfg);
dgp.theta_true   = compute_true_irf(dgp, cfg, 'analytic');
end
