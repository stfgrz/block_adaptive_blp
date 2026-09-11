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
%   .sigma2_fixed    : OPTIONAL scalar.  If non-empty, sigma2 is HELD at
%                      this value instead of being drawn.  Used by
%                      montecarlo/run_sensitivity_approximations.m to ask
%                      what the per-equation treatment of Sigma costs:
%                      the equation-by-equation sampler ignores the
%                      cross-equation correlation of Sigma, and fixing
%                      sigma2 at the system NIW posterior value is the
%                      cleanest way to see whether that matters for the
%                      tau path and hence for the point estimate.  (It
%                      cannot matter for the tau = 1 posterior mean,
%                      which is sigma2-free -- see conditional (1').)
%   .tau_probs       : OPTIONAL vector of probabilities for the tau
%                      quantile summary (default [.05 .25 .5 .75 .95]).
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
%   .beta_rb_mean (m x 1)  RAO-BLACKWELLISED posterior mean: the average
%               over retained draws of the CONDITIONAL mean
%               E[beta | tau, sigma2, y] rather than of beta itself.
%               Same estimand as .beta_mean, strictly smaller Monte
%               Carlo error (Rao-Blackwell), because the exact Gaussian
%               draw noise is integrated out instead of simulated.  With
%               tau fixed it equals .beta_cond_mean exactly.
%   .beta_cond_draws (n_keep x m)  the retained conditional means, i.e.
%               the series whose average is .beta_rb_mean; its effective
%               sample size is the honest mixing diagnostic for the
%               REPORTED point estimate.
%   .beta_cond_mean (m x 1)  conditional posterior mean of beta at the
%               LAST draw's (tau, sigma2).  When tau is FIXED and
%               prior_scales_with_sigma2 = true this quantity does not
%               depend on the iteration at all (sigma2 cancels from
%               conditional (1')), so it is the EXACT closed-form NIW
%               posterior mean and carries no Monte Carlo error --
%               tests/test_fmar_nesting_exact.m uses it for a
%               tolerance-free check of the FMAR nesting.
%   .tau2_mean/.tau2_med/.tau_q/.p_tau_gt1 : block-scale summaries
%               (tau scale for .tau_q; quantiles at opts.tau_probs)
%   .diag       struct: .n_tau_clip, .all_finite, .lag1_acorr_beta,
%               .ess_beta, .ess_logtau, .lag1_acorr_logtau
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

% --- loop-invariant objects, hoisted out of the sweep -------------------
% This block changes NOTHING about the algorithm or the order in which
% random numbers are consumed; it only stops the inner loop from
% recomputing indices and rebuilding matrices K*H times per estimator
% call.  (The sweep runs n_burn + n_keep times for every equation and
% every horizon, so the interpreter overhead of the naive form dominated
% the actual arithmetic.)
blk      = max(prior.block_id, 1);      % tau index for every column
int_idx  = find(is_int);                % intercept positions
d_int    = prior.d(int_idx);
diag_idx = 1:(m + 1):(m * m);           % linear indices of diag(Pt)
% Bsel(g, j) = 1 when regressor j is in block g: turns the per-block sums
% of conditional (3)/(3') into one matrix product.
Bsel = zeros(G, m);
p_g  = zeros(G, 1);
for g = 1:G
    Bsel(g, :) = (prior.block_id(:) == g)';
    p_g(g) = sum(prior.block_id == g);
end
d_inv = 1 ./ prior.d;
tau_shape = (p_g + 1) / 2;

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

if ~isfield(opts, 'sigma2_fixed'), opts.sigma2_fixed = []; end
fix_sig2 = ~isempty(opts.sigma2_fixed);
if fix_sig2
    assert(opts.sigma2_fixed > 0, 'gibbs_block_horseshoe: sigma2_fixed must be > 0.');
    sig2 = opts.sigma2_fixed;
end
if ~isfield(opts, 'tau_probs') || isempty(opts.tau_probs)
    opts.tau_probs = [0.05 0.25 0.50 0.75 0.95];
end

n_total = opts.n_burn + opts.n_keep;
mp_last = zeros(m, 1);              % conditional mean at the last draw
mp_sum  = zeros(m, 1);              % running sum for the Rao-Blackwell mean
mp_draws = zeros(opts.n_keep, m);   % retained conditional means (RB series)
beta_draws = zeros(opts.n_keep, m);
sig2_draws = zeros(opts.n_keep, 1);
tau2_draws = zeros(opts.n_keep, G);
n_tau_clip = 0;

for it = 1:n_total
    % Relative prior variances (shared by both variants).
    w = lam2 .* tau2(blk) .* prior.d;   % (m x 1)
    w(int_idx) = d_int;                 % intercept: no lambda, no tau
    winv = 1 ./ w;

    if ~conj
        % --- (1) beta | .  [prototype: prior does NOT scale with sig2]
        P = ZtZ / sig2;                     % posterior precision (m x m)
        P(diag_idx) = P(diag_idx) + winv';  % + W^{-1}, without diag()
        rhs = Zty / sig2 + prior.mu .* winv;   % (m x 1)
        L   = chol_lower_fast(P);
        mp  = L' \ (L \ rhs);               % posterior mean
        mp_last = mp;
        beta = mp + L' \ randn(m, 1);       % exact Gaussian draw

        % --- (2) sigma2 | . ---------------------------------------------
        if ~fix_sig2
            resid = y - Z * beta;
            SSR   = resid' * resid;
            sig2  = prior.b0 + SSR / 2;
            sig2  = sig2 / draw_gamma(prior.a0 + N / 2, 1);   % = draw_ig(.,.)
        end
    else
        % --- (1') beta | .  [conjugate: prior variance = sigma2 * w] ----
        Pt = ZtZ;                           % sigma2-free precision core
        Pt(diag_idx) = Pt(diag_idx) + winv';
        rhs = Zty + prior.mu .* winv;
        L   = chol_lower_fast(Pt);
        mp  = L' \ (L \ rhs);               % mean (independent of sigma2)
        mp_last = mp;
        beta = mp + sqrt(sig2) * (L' \ randn(m, 1));

        % --- (2') sigma2 | .  [prior contributes m pseudo-observations] -
        if ~fix_sig2
            resid = y - Z * beta;
            SSR   = resid' * resid;
            delta = beta - prior.mu;
            quad  = sum(delta.^2 .* winv);  % delta' Wt^{-1} delta
            rate  = prior.b0 + (SSR + quad) / 2;
            sig2  = rate / draw_gamma(prior.a0 + (N + m) / 2, 1);
        end
    end

    % --- (3)/(3')-(4) block scales | . -----------------------------------
    if opts.sample_tau
        delta = beta - prior.mu;        % deviation from VAR centre
        Svec = Bsel * (delta.^2 .* d_inv);      % (G x 1), all blocks at once
        if conj
            Svec = Svec / sig2;         % (3'): scales carry sigma2
        end
        for g = 1:G
            % draw_ig(a, b) is exactly b / draw_gamma(a, 1); inlining the
            % one-line wrapper removes ~2*G function calls per sweep and
            % consumes the random stream in exactly the same order.
            rate_t  = 1 / nu(g) + Svec(g) / (2 * lam2);
            tau2(g) = rate_t / draw_gamma(tau_shape(g), 1);
            if tau2(g) < opts.tau2_min || tau2(g) > opts.tau2_max
                tau2(g) = min(max(tau2(g), opts.tau2_min), opts.tau2_max);
                n_tau_clip = n_tau_clip + 1;
            end
            nu(g) = (1 + 1 / tau2(g)) / draw_gamma(1, 1);
        end
    end

    % --- Store ----------------------------------------------------------
    if it > opts.n_burn
        keep = it - opts.n_burn;
        mp_sum = mp_sum + mp_last;
        mp_draws(keep, :) = mp_last';
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
out.beta_cond_mean = mp_last;
out.beta_rb_mean   = mp_sum / opts.n_keep;
out.beta_cond_draws = mp_draws;     % the series whose average is .beta_rb_mean

% --- block-scale summaries (tau scale, not tau^2) ----------------------
tau_draws     = sqrt(tau2_draws);            % (n_keep x G)
out.tau2_mean = mean(tau2_draws, 1)';
out.tau_mean  = mean(tau_draws, 1)';
out.tau_med   = zeros(G, 1);
out.tau_q     = zeros(G, numel(opts.tau_probs));
for g = 1:G
    out.tau_med(g) = empirical_quantile(tau_draws(:, g), 0.5);
    out.tau_q(g, :) = empirical_quantile(tau_draws(:, g), opts.tau_probs);
end
out.tau_probs   = opts.tau_probs(:)';
out.p_tau_gt1   = mean(tau_draws > 1, 1)';    % P(tau_g > 1 | data)

out.diag.n_tau_clip  = n_tau_clip;
out.diag.all_finite  = all(isfinite(beta_draws(:))) && ...
                       all(isfinite(sig2_draws)) && ...
                       all(isfinite(tau2_draws(:)));
out.diag.lag1_acorr_beta = mean_lag1_autocorr(beta_draws);
out.diag.ess_beta        = min_ess(beta_draws);
logtau = log(max(tau_draws, realmin));
out.diag.lag1_acorr_logtau = colwise_lag1(logtau);
out.diag.ess_logtau        = colwise_ess(logtau);

assert(out.diag.all_finite, 'gibbs_block_horseshoe: non-finite draws.');
end

% =====================================================================
function L = chol_lower_fast(P)
% Lower Cholesky with the same numerical guarantees as
% utils/safe_chol_lower.m, but taking the fast path (a bare chol) when
% the matrix is already positive definite -- which it is on essentially
% every sweep, because P = Z'Z + W^{-1} with W^{-1} > 0.  The guarded
% routine is called only when chol actually fails, so the fallback
% behaviour and its jitter warning are unchanged.
[L, flag] = chol(P, 'lower');
if flag ~= 0
    L = safe_chol_lower(P);
end
end

function e = min_ess(D)
% Smallest effective sample size across the columns of D.
e = inf;
for j = 1:size(D, 2)
    e = min(e, ess_ips(D(:, j)));
end
end

function r = colwise_lag1(D)
% Lag-1 autocorrelation of each column, returned as a column vector.
r = zeros(size(D, 2), 1);
for j = 1:size(D, 2)
    v = D(:, j) - mean(D(:, j));
    r(j) = sum(v(1:end-1) .* v(2:end)) / (sum(v.^2) + eps);
end
end

function e = colwise_ess(D)
% Effective sample size of each column, returned as a column vector.
e = zeros(size(D, 2), 1);
for j = 1:size(D, 2)
    e(j) = ess_ips(D(:, j));
end
end

function ess = ess_ips(v)
% Initial-positive-sequence effective sample size (Geyer 1992).
n = numel(v);
v = v - mean(v);
den = sum(v.^2);
if den <= 0, ess = n; return; end
maxlag = min(n - 2, 200);
rho = zeros(1, maxlag);
for l = 1:maxlag
    rho(l) = sum(v(1:end-l) .* v(l+1:end)) / den;
end
s = 0;
for l = 1:2:maxlag - 1
    pair = rho(l) + rho(l + 1);
    if pair <= 0, break; end
    s = s + pair;
end
ess = min(max(n / (1 + 2 * s), 1), n);
end

function rho = mean_lag1_autocorr(D)
% Average (across columns) lag-1 autocorrelation of the draw matrix D
% (n_keep x m).  Simple mixing diagnostic for the Monte Carlo logs.
Dc = D - mean(D, 1);
num = sum(Dc(1:end-1, :) .* Dc(2:end, :), 1);
den = sum(Dc.^2, 1) + eps;
rho = mean(num ./ den);
end
