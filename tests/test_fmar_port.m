function test_fmar_port()
% PURPOSE
% -------
% Validate the FMAR port against the ORIGINAL replication code of
% Ferreira, Miranda-Agrippino & Ricco (github.com/leonardonferreira/BLP)
% by calling their functions and this project's ports ON IDENTICAL
% INPUTS and asserting agreement:
%
%   A. priors/niw_logml.m           vs their logMLBLP_formin.m
%      (log marginal likelihood + hyperprior objective, posterior mean
%      betahat and posterior-mode sigmahat, at several lambda values
%      and horizons);
%   B. priors/select_lambda_fmar.m  vs their maxMLikelihoodBLP.m
%      (the selected lambda* and the posterior mean at the optimum);
%   C. utils/gamma_coef.m           vs their GammaCoef.m.
%
% The FMAR code is located through the environment variable FMAR_PATH
% (pointing at their demo_IRFs folder) or a small set of default
% locations.  IF THE FMAR CODE IS NOT FOUND THE TEST IS SKIPPED WITH A
% MESSAGE rather than failed: the port's internal consistency is
% covered by the other tests; this one certifies fidelity to the
% original and is meant to be run at least once on any machine that has
% the replication package.
%
% NOTES
% -----
% * Their maxMLikelihoodBLP.m calls `addpath([cd '/subroutines'])`, so
%   the working directory is temporarily changed to FMAR_PATH.
% * Their optimiser is csminwel (quasi-Newton on a logistic transform);
%   this port uses golden-section on log(lambda).  Both maximise the
%   same smooth 1-D objective, so lambda* is compared with a loose
%   tolerance (1e-3 relative) while logML values at FIXED lambda are
%   compared to near machine precision (1e-8).
% * Inputs are synthetic data from the project's own sparse DGP so the
%   test is self-contained.

fprintf('test_fmar_port: ');

% --- locate the FMAR replication code -----------------------------------
cand = {getenv('FMAR_PATH'), ...
        '/home/claude/work/fmar_blp/demo_IRFs', ...
        fullfile(fileparts(mfilename('fullpath')), '..', 'external', 'fmar', 'demo_IRFs')};
fmar_dir = '';
for i = 1:numel(cand)
    if ~isempty(cand{i}) && exist(fullfile(cand{i}, 'subroutines', 'logMLBLP_formin.m'), 'file')
        fmar_dir = cand{i};
        break
    end
end
if isempty(fmar_dir)
    fprintf('SKIPPED (FMAR replication code not found; set FMAR_PATH).\n');
    return
end

cfg = default_config();
cfg.mode = 'fmar';
rng(2026, 'twister');
dgp = simulate_sparse_misspec_dgp(cfg);
Y = dgp.Y;
[T, K] = size(Y);
p = cfg.p;
m = 1 + K * p;
d = K + 2;

% Build one horizon problem exactly as estimate_blp_fmar does.
bvar = estimate_bvar_niw(Y, cfg);
dt = var_deterministic_trend(Y, bvar.B, p);
x  = dt.x;  Tx = size(x, 1);
Zall = build_lp_regressors(x, p);
J = [eye(K), zeros(K, K * (p - 1))];

old_dir = cd(fmar_dir);
cleanup = onCleanup(@() cd(old_dir));
addpath(fullfile(fmar_dir, 'subroutines'));

tol_ml  = 1e-8;   % identical formula, fixed inputs
tol_lam = 1e-3;   % different optimisers, same objective

for h = [2, 5, 12, 20]
    Zh = Zall(1:end - h, :);
    Yh = x(p + h:Tx, :);
    Th = size(Yh, 1);
    Fh = bvar.F^h;
    b_center = [zeros(1, K); (J * Fh)'];
    psi_h = fmar_prior_scale(x, p, h);

    % ---- A. logML at fixed lambda vs logMLBLP_formin -------------------
    % Their function works on a logistic transform of lambda with
    % bounds [1e-4, 5], adds the horizon-dependent Gamma hyperprior,
    % and returns MINUS the objective.
    MINs.lambda = cfg.fmar.lambda_min;  MAXs.lambda = cfg.fmar.lambda_max;
    MINs.psi = psi_h / 100;  MAXs.psi = psi_h * 100;   % unused (mn.psi=0)
    mn.psi = 0;  mn.alpha = 0;
    sd_h = cfg.fmar.hyper_sd_rule(h);
    priorcoef.lambda = GammaCoef(cfg.fmar.hyper_mode, sd_h, 0);
    priorcoef.priorType = 'VAR';

    for lam = [0.1, 0.4, 1.5]
        par = -log((MAXs.lambda - lam) / (lam - MINs.lambda)); % logistic^-1
        [negQ_fmar, beta_fmar, sig_fmar] = logMLBLP_formin(par, Yh, Zh, p, ...
            Th, K, b_center, MINs, MAXs, psi_h, cfg.fmar.Vc, [], mn, 0, 0, ...
            mean(x(1:p, :)), 1, priorcoef, h);

        omega = zeros(m, 1);  omega(1) = cfg.fmar.Vc;
        for lag = 1:p, omega(1 + (lag-1)*K + (1:K)) = lam^2 ./ psi_h; end
        [lml, beta_port, sig_port] = niw_logml(Yh, Zh, b_center, omega, psi_h, d);
        [gk, gtheta] = gamma_coef(cfg.fmar.hyper_mode, sd_h);
        Q_port = lml + (gk - 1) * log(lam) - lam / gtheta ...
                     - gk * log(gtheta) - gammaln(gk);

        assert(abs(-negQ_fmar - Q_port) < tol_ml * max(1, abs(Q_port)), ...
            'objective mismatch at h=%d lam=%.2f: %.10g vs %.10g', ...
            h, lam, -negQ_fmar, Q_port);
        assert(max(abs(beta_fmar(:) - beta_port(:))) < 1e-9, ...
            'betahat mismatch at h=%d lam=%.2f', h, lam);
        assert(max(abs(sig_fmar(:) - sig_port(:))) < 1e-9, ...
            'sigmahat mismatch at h=%d lam=%.2f', h, lam);
    end

    % ---- B. lambda* vs maxMLikelihoodBLP --------------------------------
    hyperPars.lambda = 0.4;  hyperPars.lambdaC = cfg.fmar.Vc;
    hyperPars.lambdaP = 0.4; hyperPars.miu = 1;  hyperPars.theta = 2;
    hyperPars.alpha = 2;     hyperPars.isrw = false(1, K);
    HPO.hyperpriors = true;  HPO.Vc = cfg.fmar.Vc;
    HPO.pos = find(~hyperPars.isrw);
    HPO.MNalpha = [];  HPO.MNpsi = false;  HPO.noc = false;  HPO.sur = false;
    HPO.Fcast = false; HPO.hz = cfg.H;     HPO.mcmc = false;
    HPO.Ndraws = 100;  HPO.Ndrawsdiscard = 50;  HPO.MCMCconst = 1;
    HPO.MCMCfcast = false;  HPO.MCMCstorecoeff = false;
    HPO.initialValues = hyperPars;  HPO.priorType = 'VAR';

    pars = maxMLikelihoodBLP(Yh, Zh, b_center, psi_h, mean(x(1:p, :)), HPO, p, h);
    lam_fmar = pars.postmax.lambda;

    [lam_port, beta_at_opt] = select_lambda_fmar(Yh, Zh, b_center, psi_h, h, cfg);

    assert(abs(lam_port - lam_fmar) < tol_lam * max(1, lam_fmar), ...
        'lambda* mismatch at h=%d: %.6f (FMAR) vs %.6f (port)', ...
        h, lam_fmar, lam_port);
    assert(max(abs(beta_at_opt(:) - pars.postmax.betahat(:))) < 1e-4, ...
        'betahat at optimum mismatch at h=%d', h);
end

% ---- C. gamma_coef vs GammaCoef ----------------------------------------
for mv = [0.4, 1.0]
    for sv = [0.1, 0.2, 0.5]
        r = GammaCoef(mv, sv, 0);
        [k2, t2] = gamma_coef(mv, sv);
        assert(abs(r.k - k2) < 1e-12 && abs(r.theta - t2) < 1e-12, ...
            'gamma_coef mismatch');
    end
end

fprintf('PASSED (port matches original FMAR code)\n');
end
