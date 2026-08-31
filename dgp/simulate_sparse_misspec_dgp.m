function dgp = simulate_sparse_misspec_dgp(cfg)
% PURPOSE
% -------
% DGP 2: SPARSE / LOCALISED misspecification of the fitted VAR(2),
% placed ON THE SHOCK-TRANSMISSION CHANNEL.  The truth is the SAME
% VAR(2) as in simulate_var_dgp.m, except for ONE extra term: variable
% 3 also responds to the THIRD lag of variable 1 (the shock variable),
%
%     y_3(t) = ... same VAR(2) terms ... + a31 * y_1(t-3) + e_3(t),
%
% i.e. the truth is a VAR(3) whose A_3 matrix is zero except entry
% (3,1) = a31.  Everything else (A_1, A_2, c, Sigma) is unchanged.
%
% WHAT IS MISSPECIFIED, AND WHY IT MATTERS FOR THE MEASURED IRF
% --------------------------------------------------------------
% * Changed: only the (3,1) entry of the lag-3 matrix, from 0 to a31.
% * The estimated model remains a VAR(2), so the omitted delayed effect
%   of y_1 on y_3 is projected onto the included lags.  Because the
%   monitored IRF is precisely the response to a y_1 shock, this
%   misspecification sits ON the transmission channel being measured:
%   the true response of y_3 has a delayed NEGATIVE correction around
%   h = 3 that no VAR(2) can reproduce (population prior bias for y_3
%   of roughly 0.2, an order of magnitude above the earlier own-lag
%   design).  The wrong prior-centre coefficients are those on LAGS OF
%   VARIABLE 1, above all in the equation for y_3, so the truly
%   misspecified block is g* = 1 (dgp.misspec_block = 1).
%
% MODEL / EQUATIONS
% -----------------
% True process: stable VAR(3) with
%   A_1, A_2 as in DGP 1;  A_3 = 0 except A_3(3,1) = a31 = -0.30.
% Stability of the true companion matrix is asserted below.
%
% INPUTS
% ------
% cfg : configuration struct; RNG state controlled by the caller.
%
% OUTPUTS
% -------
% dgp : struct with the same fields as simulate_var_dgp.m, plus
%   .misspec_block = 1 (index of the truly misspecified block).
%   The true IRF is ANALYTIC (the truth is itself a finite-order VAR),
%   computed by compute_true_irf -- never by an estimator.
%
% DIMENSIONS
% ----------
% Written for K = 3 (asserted). True lag order p* = 3 > estimated p = 2.
%
% NOTES
% -----
% a31 = -0.30 keeps max |companion eigenvalue| well below 0.95 (about
% 0.82) while producing a population VAR(2)-prior bias for the y_3
% response of roughly 0.2 -- large relative to single-sample estimation
% noise at T = 200 (roughly 0.05-0.10), so the misspecification is
% actually detectable.  The earlier design (omitted OWN third lag of
% y_3, a33 = 0.18, misspec_block = 3) produced a bias of only ~0.03,
% smaller than the noise; it is kept below, commented out, for
% comparison experiments.

assert(cfg.K == 3, ...
    'simulate_sparse_misspec_dgp: true parameters written for K = 3.');

base = base_var2_parameters();          % same A1, A2, c, B0 as DGP 1

a31 = -0.30;                            % the single misspecified entry
A3  = zeros(3);
A3(3, 1) = a31;
% --- alternative weak-signal design (kept for comparison): ----------
% a33 = 0.18;  A3 = zeros(3);  A3(3, 3) = a33;   % misspec_block = 3
A   = cat(3, base.A1, base.A2, A3);

F = companion_from_A(A);
max_eig = max(abs(eig(F)));
assert(max_eig < 0.95, ...
    'simulate_sparse_misspec_dgp: true VAR(3) unstable (max|eig|=%.3f).', max_eig);

dgp.name          = 'sparse';
dgp.type          = 'var';
dgp.A             = A;
dgp.M             = [];
dgp.c             = base.c;
dgp.Sigma         = base.B0 * base.B0';
dgp.B0            = base.B0;
dgp.misspec_block = 1;
dgp.description   = sprintf(['True VAR(3): baseline VAR(2) plus an omitted ' ...
    'delayed effect of the shock variable on y_3 (A3(3,1) = %.2f). ' ...
    'Fitted VAR(2) prior is misspecified mainly in the "lags of ' ...
    'variable 1" block of the y_3 equation.'], a31);

[dgp.Y, dgp.eta] = simulate_linear_dgp(dgp, cfg);
dgp.theta_true   = compute_true_irf(dgp, cfg, 'analytic');
end

% =====================================================================
function base = base_var2_parameters()
% Baseline VAR(2) parameters, identical to those in simulate_var_dgp.m.
% Duplicated deliberately in a single local function per DGP file so
% each DGP file is self-contained and its documentation self-reading;
% test_true_irf.m checks the two copies agree.
base.A1 = [0.55 0.10 0.05;
           0.10 0.50 0.10;
           0.05 0.10 0.45];
base.A2 = [0.15 0.00 0.00;
           0.00 0.10 0.05;
           0.00 0.05 0.15];
base.c  = [0.20; 0.10; 0.10];
base.B0 = [1.00 0.00 0.00;
           0.30 0.90 0.00;
           0.15 0.25 0.80];
end

function F = companion_from_A(A)
[K, ~, p] = size(A);
F = zeros(K * p);
for j = 1:p
    F(1:K, (j-1)*K + (1:K)) = A(:, :, j);
end
if p > 1
    F(K+1:end, 1:K*(p-1)) = eye(K * (p - 1));
end
end
