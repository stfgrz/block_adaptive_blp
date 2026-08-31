function blp = estimate_blp_global(Y, cfg, var_est)
% PURPOSE
% -------
% GLOBAL VAR-centred Bayesian Local Projection (prototype).  For every
% horizon h and every equation i, a Gaussian LP regression is combined
% with a Gaussian prior centred on the VAR-implied LP coefficients:
%
%     beta_h - mu_h^VAR ~ N(0, lambda_h^2 * D),
%
% i.e. the block-adaptive model with ALL local scales tau_g fixed at 1.
% One tightness lambda_h governs how strongly every coefficient is
% pulled towards the VAR.
%
% *** Prototype status *** -------------------------------------------
% This is a transparent Gaussian Bayesian regression, NOT a replication
% of Ferreira, Miranda-Agrippino & Ricco: their prior scaling, global
% tightness procedure and quasi-Bayesian/HAC likelihood adjustments are
% NOT implemented here.  The interfaces (var_implied_lp_prior,
% select_global_lambda, this function) are the substitution points for
% the exact published implementation.
%
% MODEL / EQUATIONS
% -----------------
% Equation i, horizon h (N_h observations):
%     y_i(t+h) = z(t)' beta + u(t+h),  u ~ N(0, sigma2 I),
%     beta ~ N(mu, W),  W = diag(w),
%     w_j = lambda_h^2 d_j (blocks), w_1 = d_1 (intercept).
% Posterior sampled with samplers/gibbs_block_horseshoe.m with
% opts.sample_tau = false, tau_fixed = 1 (semi-conjugate Gibbs over
% beta and sigma2; conditionals documented there).
%
% Structural response draws: theta_i(h) = b1n' * beta(2:1+K) per draw;
% posterior mean/median and equal-tailed credible interval at level
% cfg.ci_level are reported.  h = 0 is fixed at b1n (identification
% step shared by all estimators).
%
% INPUTS
% ------
% Y       : (T x K) data.
% cfg     : configuration struct.
% var_est : output of estimate_var (prior centre, b1n).
%
% OUTPUTS
% -------
% blp : struct with fields
%   .theta_mean/.theta_med  (K x (H+1)) posterior mean/median IRFs
%   .lo/.hi                 (K x (H+1)) credible bands (NaN at h = 0)
%   .lambda                 (K x H) tightness used per equation/horizon
%   .prior_theta            (K x (H+1)) IRF implied by the prior centre
%                           (equals the VAR IRF -- returned as a check)
%   .diag                   struct of sampler diagnostics
%   .draws                  {K x H} theta draws if cfg.blp.return_draws
%
% DIMENSIONS
% ----------
% m = 1 + K*p regressors; N_h = T - p + 1 - h observations at horizon h.
%
% NOTES
% -----
% The RNG stream is controlled by the caller (set rng once per dataset
% or per replication); no seeds are set inside.

[T, K] = size(Y);
p = cfg.p;  H = cfg.H;
alpha = 1 - cfg.ci_level;
b1n = var_est.b1n;

Zall = build_lp_regressors(Y, p);
bp   = build_block_prior(Y, cfg);

% Gibbs settings (Monte Carlo callers may have overridden these fields).
gopts.n_burn     = cfg.gibbs.n_burn;
gopts.n_keep     = cfg.gibbs.n_keep;
gopts.sample_tau = false;
gopts.tau_fixed  = 1;
gopts.tau2_min   = cfg.gibbs.tau2_min;
gopts.tau2_max   = cfg.gibbs.tau2_max;

theta_mean = zeros(K, H + 1);  theta_med = zeros(K, H + 1);
lo = nan(K, H + 1);  hi = nan(K, H + 1);
theta_mean(:, 1) = b1n;  theta_med(:, 1) = b1n;
lambda_used  = zeros(K, H);
prior_theta  = zeros(K, H + 1);
prior_theta(:, 1) = b1n;
lag1 = zeros(K, H);
draws = cell(K, H);

for h = 1:H
    Zh = Zall(1:end - h, :);
    Mu = var_implied_lp_prior(var_est, h);   % (m x K) prior centre
    prior_theta(:, h + 1) = Mu(2:1 + K, :)' * b1n;
    for i = 1:K
        yh = Y(p + h:T, i);
        prior.mu       = Mu(:, i);
        prior.d        = bp.d_all(:, i);
        prior.block_id = bp.block_id;
        prior.a0       = cfg.gibbs.a0;
        prior.b0       = cfg.gibbs.b0;
        switch cfg.blp.lambda_mode
            case 'fixed'
                prior.lambda = cfg.blp.lambda_fixed;
            case 'grid'
                prior.lambda = select_global_lambda(yh, Zh, prior.mu, ...
                                   prior.d, prior.block_id, cfg);
            otherwise
                error('estimate_blp_global: unknown lambda_mode.');
        end
        lambda_used(i, h) = prior.lambda;

        out = gibbs_block_horseshoe(yh, Zh, prior, gopts);

        % Map beta draws to structural-response draws.
        th_draws = out.beta_draws(:, 2:1 + K) * b1n;   % (n_keep x 1)
        theta_mean(i, h + 1) = mean(th_draws);
        theta_med(i, h + 1)  = empirical_quantile(th_draws, 0.5);
        q = empirical_quantile(th_draws, [alpha/2, 1 - alpha/2]);
        lo(i, h + 1) = q(1);  hi(i, h + 1) = q(2);
        lag1(i, h) = out.diag.lag1_acorr_beta;
        if cfg.blp.return_draws
            draws{i, h} = th_draws;
        end
    end
end

blp.theta_mean  = theta_mean;
blp.theta_med   = theta_med;
blp.lo          = lo;
blp.hi          = hi;
blp.lambda      = lambda_used;
blp.prior_theta = prior_theta;
blp.diag.lag1_acorr = lag1;
if cfg.blp.return_draws
    blp.draws = draws;
end

assert(all(isfinite(theta_mean(:))), 'estimate_blp_global: non-finite.');
end
