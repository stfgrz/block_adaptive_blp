function [lambda_star, betahat, sigmahat, info] = ...
    select_lambda_fmar(Y, X, b, psi, h, cfg)
% PURPOSE
% -------
% FMAR selection of the global prior tightness lambda_h at ONE horizon:
% maximise [ NIW log marginal likelihood + log Gamma hyperprior ] over
% lambda.  This replaces the prototype's select_global_lambda.m when
% cfg.mode = 'fmar' and is the published global-tightness procedure
% (Giannone-Lenza-Primiceri empirical Bayes with a hyperprior, as
% wired in FMAR's maxMLikelihoodBLP.m / logMLBLP_formin.m).
%
% MODEL / EQUATIONS
% -----------------
% Prior at horizon h (system form, common regressors, k = 1 + K*p):
%     Sigma_h ~ IW(diag(psi), d),   d = K + 2   (so d - K - 1 = 1),
%     vec(B_h) | Sigma_h ~ N( vec(b), Sigma_h (x) Omega(lambda) ),
%     omega_1 = Vc  (loose constant),
%     omega_j = lambda^2 / psi_v        (coefficient on ANY lag of
%                                        variable v; NO lag decay --
%                                        FMAR shut the Minnesota decay
%                                        off under the 'VAR' prior,
%                                        because deviations from a
%                                        VAR-implied centre have no
%                                        reason to shrink with the lag).
% Objective (FMAR, hyperpriors = 1):
%     Q(lambda) = logML(lambda) + log Gamma( lambda; k_h, theta_h ),
% where the Gamma hyperprior has mode 0.4 and a HORIZON-DEPENDENT
% standard deviation
%     sd(h) = 0.1 + 0.4 / (1 + exp(-0.3 * (h - 12)))
% (FMAR's logistic rule: the hyperprior loosens with the horizon, from
% sd ~ 0.11 at h = 1 to sd ~ 0.5 at long horizons, allowing the data to
% pull the BLP further from the VAR where the VAR is least reliable).
% For h <= 1 FMAR instead use the plain GLP hyperprior sd = 0.2; that
% case is handled by estimate_bvar_niw.m, but this function also
% accepts h = 1 (used by the block-adaptive estimator) with the same
% logistic rule.
%
% Maximisation: FMAR optimise a logistic transform of lambda with
% csminwel.  Since lambda is the ONLY free hyperparameter here (psi
% fixed, alpha inactive, sur = noc = 0), a one-dimensional GOLDEN-
% SECTION search on log(lambda) over [1e-4, 5] (the FMAR bounds) is
% used instead: derivative-free, deterministic, and transparent.  The
% objective is smooth and, in all cases inspected, unimodal on this
% interval; a coarse pre-grid guards against a bad bracket.
%
% INPUTS
% ------
% Y   : (nT x K) horizon-h dependent block (detrended data).
% X   : (nT x k) regressors (constant first).
% b   : (k x K) prior mean (VAR-implied centre; constant row = 0).
% psi : (K x 1) prior residual scales from fmar_prior_scale.m.
% h   : horizon (used only through the hyperprior sd rule).
% cfg : uses cfg.fmar.Vc, cfg.fmar.lambda_min/max,
%       cfg.fmar.hyper_mode, cfg.fmar.hyper_sd_rule (function handle).
%
% OUTPUTS
% -------
% lambda_star : maximiser of Q.
% betahat     : (k x K) posterior mean at lambda_star.
% sigmahat    : (K x K) posterior mode of Sigma at lambda_star.
% info        : struct (.Q_star, .n_eval, .at_bound flag).
%
% DIMENSIONS
% ----------
% k = 1 + K*p regressors, K equations, nT observations.
%
% NOTES
% -----
% Verified against FMAR's maxMLikelihoodBLP.m in tests/test_fmar_port.m.

[nT, K] = size(Y);  %#ok<ASGLU>
k = size(X, 2);
d = K + 2;

lo = log(cfg.fmar.lambda_min);
hi = log(cfg.fmar.lambda_max);

% Horizon-dependent Gamma hyperprior (FMAR rule).
sd_h = cfg.fmar.hyper_sd_rule(h);
[gk, gtheta] = gamma_coef(cfg.fmar.hyper_mode, sd_h);

    function omega = make_omega(lam)
        omega = zeros(k, 1);
        omega(1) = cfg.fmar.Vc;
        for lag = 1:(k - 1) / K
            omega(1 + (lag-1)*K + (1:K)) = lam^2 ./ psi;   % no lag decay
        end
    end

    function Q = objective(loglam)
        lam = exp(loglam);
        lml = niw_logml(Y, X, b, make_omega(lam), psi, d);
        Q = lml + (gk - 1) * log(lam) - lam / gtheta ...
                - gk * log(gtheta) - gammaln(gk);
    end

% Coarse pre-grid to bracket the maximum, then golden-section refine.
n_grid = 25;
grid = linspace(lo, hi, n_grid);
Qg = zeros(n_grid, 1);
for i = 1:n_grid
    Qg(i) = objective(grid(i));
end
[~, i0] = max(Qg);
a = grid(max(i0 - 1, 1));
c = grid(min(i0 + 1, n_grid));

gr = (sqrt(5) - 1) / 2;                 % golden ratio fraction
x1 = c - gr * (c - a);  x2 = a + gr * (c - a);
f1 = objective(x1);     f2 = objective(x2);
n_eval = n_grid + 2;
while (c - a) > 1e-6
    if f1 > f2
        c = x2;  x2 = x1;  f2 = f1;
        x1 = c - gr * (c - a);  f1 = objective(x1);
    else
        a = x1;  x1 = x2;  f1 = f2;
        x2 = a + gr * (c - a);  f2 = objective(x2);
    end
    n_eval = n_eval + 1;
end
loglam_star = (a + c) / 2;
lambda_star = exp(loglam_star);

[~, betahat, sigmahat] = niw_logml(Y, X, b, make_omega(lambda_star), psi, d);

info.Q_star  = objective(loglam_star);
info.n_eval  = n_eval;
info.at_bound = (abs(loglam_star - lo) < 1e-3) || (abs(loglam_star - hi) < 1e-3);
end
