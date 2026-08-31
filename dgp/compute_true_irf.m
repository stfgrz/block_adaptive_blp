function theta_true = compute_true_irf(dgp, cfg, method)
% PURPOSE
% -------
% Compute the TRUE structural impulse responses of a simulation DGP.
%
%   *** This function uses ONLY the true DGP parameters. ***
%   *** It never touches data or any estimator, so the truth ***
%   *** cannot be contaminated by the estimators under study. ***
%
% MODEL / EQUATIONS
% -----------------
% All prototype DGPs are linear:
%
%   'var'   : y(t) = c + sum_{j=1}^{p*} A_j y(t-j) + e(t)
%   'varma' : y(t) = c + sum_{j=1}^{p*} A_j y(t-j) + e(t) + M e(t-1)
%
% with e(t) = B0 * eta(t), eta(t) ~ N(0, I), B0 = chol(Sigma, 'lower').
%
% Reduced-form MA coefficients Psi_h (dy(t+h)/de(t)') satisfy
%   Psi_0 = I,
%   Psi_h = sum_{j=1}^{min(h,p*)} A_j Psi_{h-j} + 1{h==1, varma} * M.
%
% The structural shock is the recursive (Cholesky) shock to variable s =
% cfg.shock_var, NORMALISED to raise y_s by one unit on impact:
%   b1n = B0(:,s) / B0(s,s),
%   theta_true(:, h+1) = Psi_h * b1n,  h = 0..H.
%
% method = 'simulation' cross-checks the analytic answer: simulate pairs
% of paths with COMMON RANDOM NUMBERS, add the impact vector b1n to the
% reduced-form innovation at one date in the shocked path, average the
% path differences.  For linear DGPs each pair difference is exact; the
% averaging is kept so the same code covers future nonlinear DGPs.
%
% INPUTS
% ------
% dgp    : struct with fields
%          .type  : 'var' or 'varma'
%          .A     : (K x K x p*) true autoregressive matrices
%          .M     : (K x K) MA(1) matrix ('varma' only, else ignored)
%          .c     : (K x 1) intercept
%          .Sigma : (K x K) innovation covariance
% cfg    : configuration struct (uses cfg.H, cfg.shock_var,
%          cfg.true_irf.* for the simulation method).
% method : 'analytic' (default) or 'simulation'.
%
% OUTPUTS
% -------
% theta_true : (K x (H+1)) matrix; column h+1 = response at horizon h.
%
% DIMENSIONS
% ----------
% K = size(dgp.A,1); p* = size(dgp.A,3); H = cfg.H.
%
% NOTES
% -----
% Normalisation to a unit impact on the shocked variable is used
% consistently by every estimator in the project.

if nargin < 3
    method = 'analytic';
end
K  = size(dgp.A, 1);
ps = size(dgp.A, 3);
H  = cfg.H;
s  = cfg.shock_var;

B0  = safe_chol_lower(dgp.Sigma);
b1n = B0(:, s) / B0(s, s);              % unit impact on variable s

switch method
    case 'analytic'
        % --- MA recursion for Psi_h --------------------------------
        Psi = zeros(K, K, H + 1);
        Psi(:, :, 1) = eye(K);          % Psi_0
        for h = 1:H
            acc = zeros(K, K);
            for j = 1:min(h, ps)
                acc = acc + dgp.A(:, :, j) * Psi(:, :, h - j + 1);
            end
            if strcmp(dgp.type, 'varma') && h == 1
                acc = acc + dgp.M;
            end
            Psi(:, :, h + 1) = acc;
        end
        theta_true = zeros(K, H + 1);
        for h = 0:H
            theta_true(:, h + 1) = Psi(:, :, h + 1) * b1n;
        end

    case 'simulation'
        % --- Shocked vs unshocked paths, common random numbers -----
        n_paths = cfg.true_irf.n_paths;
        Tsim    = H + cfg.true_irf.t_max_extra + ps + 1;
        t0      = ps + 1;               % date the shock hits
        acc     = zeros(K, H + 1);
        for ip = 1:n_paths
            Eta = randn(K, Tsim);       % COMMON random numbers
            E   = B0 * Eta;             % reduced-form innovations
            y_base  = simulate_path(dgp, E, Tsim, K, ps);
            E_shk        = E;
            E_shk(:, t0) = E_shk(:, t0) + b1n;   % add impact vector
            y_shk   = simulate_path(dgp, E_shk, Tsim, K, ps);
            acc = acc + (y_shk(:, t0:t0 + H) - y_base(:, t0:t0 + H));
        end
        theta_true = acc / n_paths;

    otherwise
        error('compute_true_irf: unknown method "%s".', method);
end

assert(all(isfinite(theta_true(:))), 'compute_true_irf: non-finite IRF.');
end

% =====================================================================
function y = simulate_path(dgp, E, Tsim, K, ps)
% Simulate one path of the linear DGP given innovation matrix E (K x Tsim),
% starting from zeros (initial conditions cancel in the shocked-minus-
% unshocked difference of a linear model).
y = zeros(K, Tsim);
for t = (ps + 1):Tsim
    v = dgp.c + E(:, t);
    for j = 1:ps
        v = v + dgp.A(:, :, j) * y(:, t - j);
    end
    if strcmp(dgp.type, 'varma')
        v = v + dgp.M * E(:, t - 1);
    end
    y(:, t) = v;
end
end
