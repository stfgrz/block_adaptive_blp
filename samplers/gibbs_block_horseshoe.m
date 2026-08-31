function out = gibbs_block_horseshoe(y, Z, prior, opts)
% PURPOSE
% -------
% Gibbs sampler for ONE Gaussian LP regression with a VAR-centred prior
% and (optionally) GROUP-HORSESHOE local scales on the block deviations
% from the prior centre.  With opts.sample_tau = false and tau fixed at
% 1 this file IS the global VAR-centred BLP sampler, which is exactly
% the nesting property checked by tests/test_nesting.m.
%
% MODEL / EQUATIONS
% -----------------
% Likelihood (N observations, m regressors):
%     y = Z beta + u,          u ~ N(0, sigma2 I_N).
%
% Prior, written on the deviation delta = beta - mu from the VAR centre
% mu (m x 1).  Coefficient j belongs to block g(j) in {0, 1, ..., G};
% block 0 is the intercept:
%
%     delta_j ~ N(0, w_j),
%     w_j = lambda^2 * tau_{g(j)}^2 * d_j     if g(j) >= 1,
%     w_j = d_j                               if g(j) = 0 (intercept).
%
% Group-horseshoe hierarchy on the block scales (Makalic & Schmidt,
% 2016, auxiliary-variable representation of the half-Cauchy):
%
%     tau_g   ~ Cauchy+(0, 1)   represented as
%     tau_g^2 | nu_g ~ IG(1/2, 1/nu_g),    nu_g ~ IG(1/2, 1),
%
% where IG(a, b) has density  p(x) = b^a/Gamma(a) x^{-a-1} exp(-b/x).
% (Marginalising nu_g returns the half-Cauchy density for tau_g; see
% Makalic & Schmidt 2016, "A simple sampler for the horseshoe".)
%
% Error variance:  sigma2 ~ IG(a0, b0)   (weak; a0 = b0 = 0.01).
%
% TWO PRIOR-SCALING VARIANTS (opts.prior_scales_with_sigma2):
%   false (prototype default): the prior variance of delta does NOT
%       scale with sigma2 (semi-conjugate), as in the original project
%       brief.  Conditionals (1)-(4) below.
%   true  (FMAR / NIW-nesting variant): the prior variance of EVERY
%       coefficient (intercept included) scales with sigma2,
%           delta_j | sigma2, tau ~ N(0, sigma2 * w_j),
%       exactly like the Normal-Inverse-Wishart prior of Ferreira,
%       Miranda-Agrippino & Ricco, in which the coefficient prior
%       variance is Sigma_ii * omega_j.  With tau_g = 1 the conditional
%       posterior mean of beta is then
%           (Z'Z + W^{-1})^{-1} (Z'y + W^{-1} mu)
%       for ANY sigma2 -- identical to the FMAR closed-form posterior
%       mean, which gives an EXACT nesting of the block-adaptive
%       estimator in the FMAR baseline (tests/test_fmar_nesting.m).
%       Conditionals (1')-(3') below replace (1)-(3); (4) is unchanged.
%
% FULL CONDITIONAL POSTERIORS (each drawn exactly; none approximate)
% ------------------------------------------------------------------
% Notation: W = diag(w_1..w_m) built from current (lambda, tau, d).
%
% (1) beta | sigma2, tau, y:
%     Standard Gaussian linear-model update with Gaussian prior
%     N(mu, W):
%         P  = Z'Z / sigma2 + W^{-1}              (m x m precision)
%         mp = P^{-1} ( Z'y / sigma2 + W^{-1} mu ) (m x 1 mean)
%         beta | . ~ N(mp, P^{-1}).
%     Code: with lower Cholesky L (P = L L'), solve mp = L' \ (L \ rhs)
%     and draw beta = mp + L' \ z, z ~ N(0, I_m), because
%     Cov(L'^{-1} z) = L'^{-1} L^{-1} = P^{-1}.
%
% (2) sigma2 | beta, y:
%     Prior IG(a0, b0), likelihood (sigma2)^{-N/2} exp(-SSR/(2 sigma2)),
%     SSR = ||y - Z beta||^2, gives
%         sigma2 | . ~ IG( a0 + N/2,  b0 + SSR/2 ).
%
% (3) tau_g^2 | delta_g, nu_g, lambda   (only if opts.sample_tau):
%     Block g has p_g coefficients; let S_g = sum_{j in g} delta_j^2/d_j.
%     Likelihood contribution: (tau_g^2)^{-p_g/2} exp( -S_g/(2 lambda^2 tau_g^2) ).
%     Prior IG(1/2, 1/nu_g): (tau_g^2)^{-3/2} exp( -1/(nu_g tau_g^2) ).
%     Product is an inverse-gamma kernel:
%         tau_g^2 | . ~ IG( (p_g + 1)/2,  1/nu_g + S_g/(2 lambda^2) ).
%
% (4) nu_g | tau_g^2:
%     Prior nu_g ~ IG(1/2, 1) (the Makalic-Schmidt choice; this is what
%     makes the marginal of tau_g exactly half-Cauchy(0,1)).
%     Two nu_g-dependent factors:
%       from tau_g^2 | nu_g ~ IG(1/2, 1/nu_g):
%           (1/nu_g)^{1/2} exp( -1/(nu_g tau_g^2) )
%       from the prior:
%           nu_g^{-3/2} exp( -1/nu_g )
%     Product:  nu_g^{-2} exp( -(1 + 1/tau_g^2)/nu_g ),  i.e.
%         nu_g | . ~ IG( 1,  1 + 1/tau_g^2 ).
%
% CONDITIONALS OF THE CONJUGATE VARIANT (prior_scales_with_sigma2 = true)
% -----------------------------------------------------------------------
% Notation: Wt = diag(wt), wt_j = lambda^2 tau_{g(j)}^2 d_j for blocks,
% wt_1 = d_1 for the intercept (in FMAR mode d_1 = Vc); prior is
% delta | sigma2 ~ N(0, sigma2 * Wt).
%
% (1') beta | sigma2, tau, y:
%      Both likelihood and prior precision carry 1/sigma2, so it
%      cancels from the mean:
%          Pt = Z'Z + Wt^{-1}                      (m x m)
%          mp = Pt^{-1} ( Z'y + Wt^{-1} mu )       (independent of sigma2)
%          beta | . ~ N( mp, sigma2 * Pt^{-1} ).
%      Code: L = chol(Pt) lower; beta = mp + sqrt(sigma2) * (L' \ z).
%
% (2') sigma2 | beta, tau, y:
%      The prior on delta now contributes m further "observations":
%      p(sigma2|.) ~ (sigma2)^{-N/2} e^{-SSR/(2 s2)}
%                    * (sigma2)^{-m/2} e^{-delta' Wt^{-1} delta/(2 s2)}
%                    * IG(a0, b0)
%          => sigma2 | . ~ IG( a0 + (N + m)/2,
%                              b0 + [SSR + delta' Wt^{-1} delta]/2 ).
%
% (3') tau_g^2 | delta, sigma2, nu_g, lambda:
%      delta_j | . ~ N(0, sigma2 lambda^2 tau_g^2 d_j) implies the same
%      inverse-gamma kernel as (3) with S_g replaced by
%          St_g = sum_{j in g} delta_j^2 / (sigma2 * d_j):
%          tau_g^2 | . ~ IG( (p_g + 1)/2, 1/nu_g + St_g/(2 lambda^2) ).
%
% Numerical guard: tau_g^2 draws are clipped to
% [opts.tau2_min, opts.tau2_max] (default [1e-10, 1e8]).  The horseshoe
% has heavy tails; this WIDE truncation changes the distribution only
% in regions irrelevant to the estimand but prevents overflow in long
% Monte Carlo runs.  Any clipping is counted in out.diag.n_tau_clip so
% it is visible, never silent.
%
% INPUTS
% ------
% y     : (N x 1) dependent variable.
% Z     : (N x m) regressors.
% prior : struct
%   .mu       (m x 1) VAR-implied prior centre
%   .d        (m x 1) diagonal of D
%   .block_id (m x 1) block labels, 0 = intercept
%   .lambda   scalar  global tightness lambda_h (> 0)
%   .a0, .b0  scalars IG prior for sigma2
% opts  : struct
%   .n_burn, .n_keep : burn-in and retained draws
%   .sample_tau      : true = grouped horseshoe; false = tau fixed
%   .tau_fixed       : scalar used when sample_tau = false (usually 1)
%   .tau2_min/max    : numerical clip bounds
%   .prior_scales_with_sigma2 : OPTIONAL, default false.  If true, the
%                      prior variance of every coefficient scales with
%                      sigma2 (conjugate / FMAR-NIW variant; see the
%                      primed conditionals above).
%   .seed            : OPTIONAL; if present, rng(opts.seed) is set here.
%                      Otherwise the caller controls the RNG stream.
%
% OUTPUTS
% -------
% out : struct
%   .beta_draws (n_keep x m)   retained beta draws
%   .sig2_draws (n_keep x 1)
%   .tau2_draws (n_keep x G)   (all equal tau_fixed^2 if not sampled)
%   .beta_mean  (m x 1)
%   .diag       struct: .n_tau_clip, .all_finite, .lag1_acorr_beta
%
% DIMENSIONS
% ----------
% N observations, m regressors, G = max(block_id) blocks.
%
% NOTES
% -----
% Initialisation at OLS (ridge-regularised if needed) and tau = 1.
% Every draw routes through utils/draw_gamma.m -> rand/randn, so the
% whole sampler is reproducible under a single rng(seed).

% --- Checks and setup --------------------------------------------------
[N, m] = size(Z);
assert(numel(y) == N, 'gibbs_block_horseshoe: y/Z size mismatch.');
assert(numel(prior.mu) == m && numel(prior.d) == m && ...
       numel(prior.block_id) == m, 'gibbs_block_horseshoe: prior sizes.');
assert(prior.lambda > 0, 'gibbs_block_horseshoe: lambda must be > 0.');
assert(all(prior.d > 0), 'gibbs_block_horseshoe: d must be > 0.');

if isfield(opts, 'seed') && ~isempty(opts.seed)
    rng(opts.seed);
end

y  = y(:);
G  = max(prior.block_id);
lam2 = prior.lambda^2;
is_int = (prior.block_id == 0);         % logical: intercept coefficients
if ~isfield(opts, 'prior_scales_with_sigma2')
    opts.prior_scales_with_sigma2 = false;   % prototype default
end
conj = opts.prior_scales_with_sigma2;

ZtZ = Z' * Z;                           % (m x m), precomputed once
Zty = Z' * y;                           % (m x 1)

% Initial values: (ridge) OLS for beta, its residual variance, tau = 1.
beta = (ZtZ + 1e-8 * trace(ZtZ) / m * eye(m)) \ Zty;
sig2 = max((y - Z * beta)' * (y - Z * beta) / max(N - m, 1), 1e-10);
if opts.sample_tau
    tau2 = ones(G, 1);                  % start at the global prior
else
    assert(isfield(opts, 'tau_fixed'), ...
        'gibbs_block_horseshoe: tau_fixed required when sample_tau = false.');
    tau2 = ones(G, 1) * opts.tau_fixed^2;
end
nu = ones(G, 1);

n_total = opts.n_burn + opts.n_keep;
beta_draws = zeros(opts.n_keep, m);
sig2_draws = zeros(opts.n_keep, 1);
tau2_draws = zeros(opts.n_keep, G);
n_tau_clip = 0;

for it = 1:n_total
    % Relative prior variances (shared by both variants).
    w = lam2 .* tau2(max(prior.block_id, 1)) .* prior.d;  % (m x 1)
    w(is_int) = prior.d(is_int);        % intercept: no lambda, no tau

    if ~conj
        % --- (1) beta | .  [prototype: prior does NOT scale with sig2]
        P   = ZtZ / sig2 + diag(1 ./ w);    % posterior precision (m x m)
        rhs = Zty / sig2 + prior.mu ./ w;   % (m x 1)
        L   = safe_chol_lower(P);
        mp  = L' \ (L \ rhs);               % posterior mean
        beta = mp + L' \ randn(m, 1);       % exact Gaussian draw

        % --- (2) sigma2 | . ---------------------------------------------
        resid = y - Z * beta;
        SSR   = resid' * resid;
        sig2  = draw_ig(prior.a0 + N / 2, prior.b0 + SSR / 2);
    else
        % --- (1') beta | .  [conjugate: prior variance = sigma2 * w] ----
        Pt  = ZtZ + diag(1 ./ w);           % sigma2-free precision core
        rhs = Zty + prior.mu ./ w;
        L   = safe_chol_lower(Pt);
        mp  = L' \ (L \ rhs);               % mean (independent of sigma2)
        beta = mp + sqrt(sig2) * (L' \ randn(m, 1));

        % --- (2') sigma2 | .  [prior contributes m pseudo-observations] -
        resid = y - Z * beta;
        SSR   = resid' * resid;
        delta = beta - prior.mu;
        quad  = sum(delta.^2 ./ w);         % delta' Wt^{-1} delta
        sig2  = draw_ig(prior.a0 + (N + m) / 2, ...
                        prior.b0 + (SSR + quad) / 2);
    end

    % --- (3)/(3')-(4) block scales | . -----------------------------------
    if opts.sample_tau
        delta = beta - prior.mu;        % deviation from VAR centre
        for g = 1:G
            idx = (prior.block_id == g);
            p_g = sum(idx);
            S_g = sum(delta(idx).^2 ./ prior.d(idx));
            if conj
                S_g = S_g / sig2;       % (3'): scales carry sigma2
            end
            tau2(g) = draw_ig((p_g + 1) / 2, 1 / nu(g) + S_g / (2 * lam2));
            if tau2(g) < opts.tau2_min || tau2(g) > opts.tau2_max
                tau2(g) = min(max(tau2(g), opts.tau2_min), opts.tau2_max);
                n_tau_clip = n_tau_clip + 1;
            end
            nu(g) = draw_ig(1, 1 + 1 / tau2(g));
        end
    end

    % --- Store ----------------------------------------------------------
    if it > opts.n_burn
        keep = it - opts.n_burn;
        beta_draws(keep, :) = beta';
        sig2_draws(keep)    = sig2;
        tau2_draws(keep, :) = tau2';
    end
end

% --- Diagnostics --------------------------------------------------------
out.beta_draws = beta_draws;
out.sig2_draws = sig2_draws;
out.tau2_draws = tau2_draws;
out.beta_mean  = mean(beta_draws, 1)';
out.diag.n_tau_clip  = n_tau_clip;
out.diag.all_finite  = all(isfinite(beta_draws(:))) && ...
                       all(isfinite(sig2_draws)) && ...
                       all(isfinite(tau2_draws(:)));
out.diag.lag1_acorr_beta = mean_lag1_autocorr(beta_draws);

assert(out.diag.all_finite, 'gibbs_block_horseshoe: non-finite draws.');
end

% =====================================================================
function rho = mean_lag1_autocorr(D)
% Average (across columns) lag-1 autocorrelation of the draw matrix D
% (n_keep x m).  Simple mixing diagnostic for the Monte Carlo logs.
Dc = D - mean(D, 1);
num = sum(Dc(1:end-1, :) .* Dc(2:end, :), 1);
den = sum(Dc.^2, 1) + eps;
rho = mean(num ./ den);
end
