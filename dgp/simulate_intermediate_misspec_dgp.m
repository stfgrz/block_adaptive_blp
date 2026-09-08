function dgp = simulate_intermediate_misspec_dgp(cfg)
% PURPOSE
% -------
% DGP 4: INTERMEDIATE misspecification of the fitted VAR(2) -- halfway
% between the sparse design (one wrong coefficient, one equation) and
% the dense design (every block of every equation wrong).
%
% The truth is the same baseline VAR(2) plus a lag-3 matrix whose only
% nonzero entries sit in COLUMN 1, i.e. EVERY equation omits a delayed
% effect of the shock variable y_1:
%
%     y(t) = c + A_1 y(t-1) + A_2 y(t-2) + A_3 y(t-3) + e(t),
%     A_3  = [a * l_1, 0, 0;  a * l_2, 0, 0;  a * l_3, 0, 0].
%
% WHAT IS MISSPECIFIED, AND WHY THIS IS "INTERMEDIATE"
% ----------------------------------------------------
% * Sparse DGP:       one entry  (3,1)  -> block 1 wrong in ONE equation.
% * Intermediate:     one column (:,1)  -> block 1 wrong in ALL equations.
% * Dense DGP:        VARMA(2,1) with dense M -> EVERY block wrong in
%                     every equation, so there is no block to escape to.
% The block-adaptive prior's target (the "lags of variable 1" block,
% g* = 1) is therefore still well defined here, but the signal is
% spread over K equations instead of concentrated in one.  This is the
% design that separates "block-sparse in the coefficient space" from
% "block-sparse in the equation space": if the adaptive layer only ever
% helps when the conflict is confined to a single equation, this DGP
% will show it.
%
% MODEL / EQUATIONS
% -----------------
% A_1, A_2, c, Sigma as in simulate_var_dgp.m; A_3 as above with
%     a    = cfg.dgp.interm_scale (default -0.30)
%     l    = cfg.dgp.interm_load  (default [0.5 0.7 1.0], so the y_3
%            equation carries the same -0.30 entry as the sparse DGP
%            and the other two equations carry weaker versions).
% Stability of the true companion matrix is asserted.
%
% INPUTS
% ------
% cfg : configuration struct; RNG state controlled by the caller.
%       Optional fields cfg.dgp.interm_scale, cfg.dgp.interm_load.
%
% OUTPUTS
% -------
% dgp : same fields as the other DGP files, with
%       .misspec_block = 1 and .params recording the design constants.
%
% DIMENSIONS
% ----------
% Written for K = 3 (asserted). True lag order p* = 3 > estimated p = 2.
%
% NOTES
% -----
% The true IRF is ANALYTIC (a finite-order VAR), computed by
% compute_true_irf -- never by an estimator.

assert(cfg.K == 3, ...
    'simulate_intermediate_misspec_dgp: true parameters written for K = 3.');

a = dgp_param(cfg, 'interm_scale', -0.30);
l = dgp_param(cfg, 'interm_load',  [0.5; 0.7; 1.0]);
l = l(:);
assert(numel(l) == 3, ...
    'simulate_intermediate_misspec_dgp: cfg.dgp.interm_load needs 3 elements.');

base = base_var2_parameters();
A3 = zeros(3);
A3(:, 1) = a * l;
A  = cat(3, base.A1, base.A2, A3);

F = companion_from_A(A);
max_eig = max(abs(eig(F)));
assert(max_eig < 0.95, ...
    'simulate_intermediate_misspec_dgp: true VAR(3) unstable (max|eig|=%.3f).', ...
    max_eig);

dgp.name          = 'intermediate';
dgp.type          = 'var';
dgp.A             = A;
dgp.M             = [];
dgp.c             = base.c;
dgp.Sigma         = base.B0 * base.B0';
dgp.B0            = base.B0;
dgp.misspec_block = 1;
dgp.params        = struct('interm_scale', a, 'interm_load', l');
dgp.description   = sprintf(['True VAR(3): baseline VAR(2) plus an omitted ' ...
    'delayed effect of the shock variable in EVERY equation ' ...
    '(A3(:,1) = %.2f * [%.2f %.2f %.2f]). The fitted VAR(2) prior is ' ...
    'misspecified in the "lags of variable 1" block of all K equations ' ...
    '(intermediate between the sparse and dense designs).'], ...
    a, l(1), l(2), l(3));

[dgp.Y, dgp.eta] = simulate_linear_dgp(dgp, cfg);
dgp.theta_true   = compute_true_irf(dgp, cfg, 'analytic');
end

% =====================================================================
function v = dgp_param(cfg, name, default)
% Read cfg.dgp.<name> if present, otherwise the documented default.
% Keeps every earlier configuration struct valid.
v = default;
if isfield(cfg, 'dgp') && isstruct(cfg.dgp) && isfield(cfg.dgp, name) ...
        && ~isempty(cfg.dgp.(name))
    v = cfg.dgp.(name);
end
end

function base = base_var2_parameters()
% Baseline VAR(2) parameters, identical to simulate_var_dgp.m.
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
