function out = gibbs_block_pooled_horizons(hdata, prior, opts)
% PURPOSE
% -------
% Gibbs / Metropolis-within-Gibbs sampler for ONE equation's local
% projections at ALL horizons h = 1..H JOINTLY, with block escape
% scales tau_{g,h} that are SMOOTHED ACROSS HORIZONS.
%
% Motivation: in samplers/gibbs_block_horseshoe.m each tau_{g,h} is
% estimated from a single horizon-h regression, i.e. from p_g = p
% coefficient deviations.  That is very little information, so the
% independent tau path is noisy and the induced extra variance can
% swamp the bias reduction it buys.  Here the tau path is given a
% first-difference (random-walk) shrinkage prior on the LOG scale, so
% neighbouring horizons borrow strength while a genuinely localised
% escape can still develop gradually.
%
% MODEL / EQUATIONS
% -----------------
% For each horizon h (N_h observations, m regressors):
%
%     y_h        = Z_h beta_h + u_h,     u_h ~ N(0, sigma2_h I),
%     beta_h - mu_h | sigma2_h, tau ~ N(0, sigma2_h W_h),
%     W_h  = diag(w_h),
%     w_{h,j} = lambda_h^2 tau_{g(j),h}^2 d_{h,j}   (g(j) >= 1),
%     w_{h,1} = d_{h,1}                             (intercept, block 0),
%     sigma2_h ~ IG(a0_h, b0_h).
%
% This is EXACTLY the conjugate ("prior_scales_with_sigma2 = true")
% variant of gibbs_block_horseshoe applied horizon by horizon; the ONLY
% change is the prior on the tau path.  Conditional on tau, the
% horizons remain independent, so steps (1)-(2) below are the same
% per-horizon conjugate updates.
%
% HORIZON-POOLING PRIOR (the new object)
% --------------------------------------
% Write x_{g,h} = log tau_{g,h}.  For each block g:
%
%     x_{g,1}  ~ half-Cauchy(0,1) on the tau scale, i.e.
%                p(x) ∝ exp(x) / (1 + exp(2x))            [same anchor as
%                                                          the independent
%                                                          estimator]
%     x_{g,h}  = phi * x_{g,h-1} + e_{g,h},  e ~ N(0, kappa_g^2), h >= 2,
%     kappa_g^2 ~ IG(a_kappa, b_kappa).
%
% phi = 1 (default) is a RANDOM WALK on log tau, i.e. a first-difference
% shrinkage prior; phi < 1 gives a stationary AR(1) that also pulls the
% path back towards tau = 1.  kappa_g controls how much smoothing:
%   kappa -> 0    : full pooling, one tau_g for all horizons;
%   kappa -> inf  : the independent-by-horizon estimator (up to the
%                   different h = 1 anchor);
% and sampling kappa_g lets the data choose, per block.
%
% FULL CONDITIONALS
% -----------------
% (1) beta_h | sigma2_h, tau, y   (exact, conjugate)
%       Pt  = Z_h'Z_h + W_h^{-1},
%       mp  = Pt^{-1} (Z_h'y_h + W_h^{-1} mu_h)      [sigma2-free]
%       beta_h ~ N(mp, sigma2_h Pt^{-1}).
%
% (2) sigma2_h | beta_h, tau, y   (exact)
%       sigma2_h ~ IG( a0_h + (N_h + m)/2,
%                      b0_h + [SSR_h + delta_h' W_h^{-1} delta_h]/2 ).
%
% (3) x_{g,h} | everything else   (single-site random-walk METROPOLIS;
%     the only non-exact step in the sampler)
%       log p(x) = -p_g x - S_{g,h} exp(-2x) / (2 lambda_h^2)
%                  + log p_prior(x | x_{g,h-1}, x_{g,h+1}, kappa_g),
%       S_{g,h}  = sum_{j in g} delta_{h,j}^2 / (sigma2_h d_{h,j}),
%     proposal x' = x + step * randn, accepted with the usual symmetric
%     Metropolis ratio.  `step` is adapted during BURN-IN ONLY (towards
%     a 0.30 acceptance rate), so the retained chain has a fixed,
%     valid transition kernel.
%     IMPLEMENTATION.  Sites at ODD horizons are conditionally
%     independent given the sites at EVEN horizons (the prior links only
%     neighbours, and blocks never interact), so the sweep is done as
%     two CHECKERBOARD half-sweeps, each a single vectorised
%     accept/reject over all blocks and all sites of one parity.  This
%     is exactly the same Markov kernel as updating the sites one at a
%     time in that order, at a fraction of the cost.
%
% (3b) LEVEL MOVE on the whole path (also Metropolis).  Under a
%     random-walk prior the LEVEL of x_{g,.} is the slow direction:
%     single-site moves have to drag H strongly dependent coordinates
%     one at a time, and the effective sample size of log tau collapses.
%     After each single-site sweep the sampler therefore proposes
%         x_{g,h} <- x_{g,h} + c   for ALL h simultaneously,  c ~ N(0, s_L^2),
%     accepted with the full path density ratio.  With phi = 1 the RW
%     increments are invariant to a common shift, so only the h = 1
%     anchor and the likelihood terms move; the code evaluates the
%     complete path density anyway, so it stays correct for phi < 1.
%     This is a symmetric proposal, so no Jacobian or proposal ratio is
%     needed, and s_L is likewise adapted during burn-in only.
%
% (4) kappa_g^2 | x   (exact, conjugate)
%       kappa_g^2 ~ IG( a_kappa + (H-1)/2,
%                       b_kappa + sum_{h=2}^{H} (x_{g,h} - phi x_{g,h-1})^2 / 2 ).
%
% INPUTS
% ------
% hdata : (1 x H) struct array, one entry per horizon, with fields
%   .y  (N_h x 1), .Z (N_h x m), .mu (m x 1), .d (m x 1),
%   .lambda scalar > 0, .a0, .b0 scalars (IG prior for sigma2_h).
% prior : struct
%   .block_id (m x 1) block labels, 0 = intercept (shared across h).
% opts  : struct
%   .n_burn, .n_keep      : burn-in and retained draws
%   .sample_tau           : true = pooled horseshoe/RW prior;
%                           false = tau fixed at opts.tau_fixed (used by
%                           the nesting tests)
%   .tau_fixed            : scalar, used when sample_tau = false
%   .phi                  : AR coefficient on log tau (default 1 = RW)
%   .a_kappa, .b_kappa    : IG prior for kappa_g^2 (default 2, 0.1 ->
%                           prior mean 0.1, i.e. a typical one-horizon
%                           move of about 0.3 in log tau)
%   .kappa_fixed          : OPTIONAL scalar; if present kappa is held
%                           there instead of sampled (sensitivity runs)
%   .mh_step              : initial single-site Metropolis step (default 0.5)
%   .mh_step_level        : initial LEVEL-move step (default 0.3)
%   .n_tau_sweeps         : Metropolis sweeps of the whole tau block per
%                           Gibbs iteration (default 5).  The beta step
%                           costs H Cholesky factorisations of an m x m
%                           matrix; a tau sweep costs G*(H+1) scalar
%                           density evaluations, i.e. almost nothing, so
%                           repeating it is much cheaper than lengthening
%                           the chain.  Raising it multiplies the
%                           effective sample size of log tau at roughly
%                           constant cost; the acceptance rates and ESS
%                           in out.diag say whether it is enough.
%   .x_min, .x_max        : clip bounds on log tau (default log of the
%                           tau2 bounds / 2, i.e. the same wide guard
%                           rails as the independent sampler)
%   .proj                 : OPTIONAL (m x 1) projection vector; if given,
%                           out.proj_draws(:, h) = beta_h' * proj per
%                           draw (used for IRF draws without storing the
%                           full beta path).
%   .fixed_blocks         : OPTIONAL block indices whose log tau path is
%                           HELD AT 0 (tau = 1, the global prior) at every
%                           horizon; the other blocks are sampled.  Same
%                           purpose as in gibbs_block_horseshoe.m (a block
%                           that violates the common-scale assumption,
%                           such as a white-noise instrument's lag block).
%   .seed                 : OPTIONAL; rng(opts.seed) is set here.
%
% OUTPUTS
% -------
% out : struct
%   .beta_mean   (m x H)             posterior mean coefficients
%   .beta_rb_mean (m x H)            RAO-BLACKWELLISED posterior mean: the
%                                    average over retained draws of
%                                    E[beta_h | tau, sigma2, y] instead of
%                                    of beta_h itself.  Same estimand as
%                                    .beta_mean with strictly less Monte
%                                    Carlo error; equals .beta_cond_mean
%                                    when tau is fixed.
%   .beta_cond_mean (m x H)          conditional mean at the LAST draw's
%                                    (tau, sigma2); with fixed tau this is
%                                    the exact closed-form posterior mean
%   .sig2_mean   (1 x H)
%   .tau2_draws  (n_keep x G x H)
%   .proj_draws  (n_keep x H)        empty unless opts.proj is given
%   .kappa_draws (n_keep x G)
%   .diag        struct: .acc_rate (G x H), .acc_rate_level (G x 1),
%                .lag1_acorr_logtau (G x H), .ess_logtau (G x H),
%                .lag1_acorr_beta (1 x H), .ess_proj (1 x H),
%                .n_tau_clip, .all_finite, .mh_step, .mh_step_level,
%                .n_tau_sweeps
%
% DIMENSIONS
% ----------
% H horizons, m regressors, G = max(block_id) blocks.
%
% NOTES
% -----
% Every random number routes through rand/randn (via utils/draw_ig.m),
% so the whole sampler is reproducible under a single rng(seed).

H = numel(hdata);
m = numel(prior.block_id);
G = max(prior.block_id);
assert(H >= 1 && m > 0 && G >= 1, 'gibbs_block_pooled_horizons: bad sizes.');

if isfield(opts, 'seed') && ~isempty(opts.seed)
    rng(opts.seed);
end

opts = set_default(opts, 'phi',        1.0);
opts = set_default(opts, 'a_kappa',    2.0);
opts = set_default(opts, 'b_kappa',    0.1);
opts = set_default(opts, 'mh_step',    0.5);
opts = set_default(opts, 'mh_step_level', 0.3);
opts = set_default(opts, 'n_tau_sweeps', 5);
opts = set_default(opts, 'x_min',     -11.5);   % tau  ~ 1e-5
opts = set_default(opts, 'x_max',       9.2);   % tau  ~ 1e4
opts = set_default(opts, 'sample_tau', true);
opts = set_default(opts, 'proj',       []);
opts = set_default(opts, 'fixed_blocks', []);
phi  = opts.phi;
free_block = true(G, 1);
if ~isempty(opts.fixed_blocks)
    fb = opts.fixed_blocks(:)';
    assert(all(fb >= 1 & fb <= G & fb == round(fb)), ...
        'gibbs_block_pooled_horizons: opts.fixed_blocks must index blocks 1..%d.', G);
    free_block(fb) = false;
end

is_int = (prior.block_id == 0);
blk    = max(prior.block_id, 1);        % index into tau for every column
% Bsel(g, j) = 1 if regressor j belongs to block g: turns the per-block
% sums of the tau conditional into one matrix product.
Bsel  = zeros(G, m);
p_g_vec = zeros(G, 1);
for g = 1:G
    Bsel(g, :) = (prior.block_id(:) == g)';
    p_g_vec(g) = sum(prior.block_id == g);
end

% --- precompute per-horizon cross-products ------------------------------
ZtZ = cell(1, H);  Zty = cell(1, H);  Nh = zeros(1, H);
lam2 = zeros(1, H);
for h = 1:H
    assert(numel(hdata(h).mu) == m && numel(hdata(h).d) == m, ...
        'gibbs_block_pooled_horizons: prior size mismatch at h = %d.', h);
    assert(hdata(h).lambda > 0 && all(hdata(h).d > 0), ...
        'gibbs_block_pooled_horizons: lambda/d must be positive at h = %d.', h);
    ZtZ{h} = hdata(h).Z' * hdata(h).Z;
    Zty{h} = hdata(h).Z' * hdata(h).y;
    Nh(h)  = numel(hdata(h).y);
    lam2(h) = hdata(h).lambda^2;
end

% --- initial values ------------------------------------------------------
beta = zeros(m, H);  sig2 = ones(1, H);
for h = 1:H
    r = ZtZ{h} + 1e-8 * trace(ZtZ{h}) / m * eye(m);
    beta(:, h) = r \ Zty{h};
    e = hdata(h).y - hdata(h).Z * beta(:, h);
    sig2(h) = max(e' * e / max(Nh(h) - m, 1), 1e-10);
end

if opts.sample_tau
    x = zeros(G, H);                       % log tau = 0  <=>  tau = 1
else
    assert(isfield(opts, 'tau_fixed'), ...
        'gibbs_block_pooled_horizons: tau_fixed required when sample_tau = false.');
    x = log(opts.tau_fixed) * ones(G, H);
end
if isfield(opts, 'kappa_fixed') && ~isempty(opts.kappa_fixed)
    kappa2 = opts.kappa_fixed^2 * ones(G, 1);
    sample_kappa = false;
else
    kappa2 = 0.1 * ones(G, 1);
    sample_kappa = true;
end

step   = opts.mh_step * ones(G, H);
stepL  = opts.mh_step_level * ones(G, 1);
n_acc  = zeros(G, H);   n_try = zeros(G, H);
n_acc_burn = zeros(G, H);  n_try_burn = zeros(G, H);
n_accL = zeros(G, 1);   n_tryL = zeros(G, 1);
n_accL_burn = zeros(G, 1);  n_tryL_burn = zeros(G, 1);
n_clip = 0;

n_total    = opts.n_burn + opts.n_keep;
tau2_draws = zeros(opts.n_keep, G, H);
kap_draws  = zeros(opts.n_keep, G);
beta_sum   = zeros(m, H);
if ~isempty(opts.proj)
    proj_draws = zeros(opts.n_keep, H);
else
    proj_draws = [];
end
x_draws = zeros(opts.n_keep, G, H);
mp_last = zeros(m, H);
mp_sum  = zeros(m, H);              % running sum for the Rao-Blackwell mean
if ~isempty(opts.proj)
    rb_draws = zeros(opts.n_keep, H);   % projected conditional means
else
    rb_draws = [];
end

for it = 1:n_total
    % ---- (1)-(2): per-horizon conjugate updates ------------------------
    for h = 1:H
        w = lam2(h) .* exp(2 * x(blk, h)) .* hdata(h).d;
        w(is_int) = hdata(h).d(is_int);
        Pt  = ZtZ{h} + diag(1 ./ w);
        rhs = Zty{h} + hdata(h).mu ./ w;
        L   = safe_chol_lower(Pt);
        mp  = L' \ (L \ rhs);
        mp_last(:, h) = mp;
        bh  = mp + sqrt(sig2(h)) * (L' \ randn(m, 1));
        beta(:, h) = bh;

        resid = hdata(h).y - hdata(h).Z * bh;
        delta = bh - hdata(h).mu;
        quad  = sum(delta.^2 ./ w);
        sig2(h) = draw_ig(hdata(h).a0 + (Nh(h) + m) / 2, ...
                          hdata(h).b0 + (resid' * resid + quad) / 2);
    end

    % ---- (3): Metropolis on the log-tau path ---------------------------
    % Repeated opts.n_tau_sweeps times per Gibbs iteration: the tau block
    % is the slow-mixing part of the chain and by far the cheapest to
    % update, so extra sweeps buy effective sample size almost free.
    if opts.sample_tau
        % Sufficient statistics S(g,h) = sum_{j in g} delta_{h,j}^2 /
        % (sigma2_h d_{h,j}); they depend on (beta, sigma2) only, so they
        % are computed ONCE per Gibbs iteration and reused by every sweep.
        S = zeros(G, H);
        for h = 1:H
            dl = beta(:, h) - hdata(h).mu;
            S(:, h) = Bsel * (dl.^2 ./ (sig2(h) * hdata(h).d));
        end
        lam2r = repmat(lam2, G, 1);             % (G x H)

        for sweep = 1:opts.n_tau_sweeps
            % --- two checkerboard half-sweeps over the horizons --------
            for parity = 1:2
                idx = parity:2:H;
                if isempty(idx), continue; end
                xc = x(:, idx);
                xp = xc + step(:, idx) .* randn(G, numel(idx));
                out_of_range = (xp < opts.x_min) | (xp > opts.x_max);
                n_clip = n_clip + sum(out_of_range(:));
                xp = min(max(xp, opts.x_min), opts.x_max);

                lp_c = site_logpost(xc, idx, x, p_g_vec, S(:, idx), ...
                                    lam2r(:, idx), phi, kappa2, H);
                lp_p = site_logpost(xp, idx, x, p_g_vec, S(:, idx), ...
                                    lam2r(:, idx), phi, kappa2, H);
                acc = log(rand(G, numel(idx))) < (lp_p - lp_c);
                acc(~free_block, :) = false;          % held blocks never move
                xc(acc) = xp(acc);
                x(:, idx) = xc;
                n_try(:, idx) = n_try(:, idx) + 1;
                n_acc(:, idx) = n_acc(:, idx) + acc;
                if it <= opts.n_burn
                    n_try_burn(:, idx) = n_try_burn(:, idx) + 1;
                    n_acc_burn(:, idx) = n_acc_burn(:, idx) + acc;
                end
            end

            % --- (3b) joint LEVEL move on the whole path, per block ----
            c = stepL .* randn(G, 1);
            xnew = x + c(:, ones(1, H));
            ok = all(xnew >= opts.x_min, 2) & all(xnew <= opts.x_max, 2);
            lp_c = path_logpost(x,    p_g_vec, S, lam2r, phi, kappa2);
            lp_p = path_logpost(xnew, p_g_vec, S, lam2r, phi, kappa2);
            accL = ok & (log(rand(G, 1)) < (lp_p - lp_c)) & free_block;
            x(accL, :) = xnew(accL, :);
            n_tryL = n_tryL + ok;
            n_accL = n_accL + accL;
            if it <= opts.n_burn
                n_tryL_burn = n_tryL_burn + ok;
                n_accL_burn = n_accL_burn + accL;
            end

            % --- (4) kappa_g^2 | x  (conjugate) ------------------------
            if sample_kappa && H > 1
                e = x(:, 2:H) - phi * x(:, 1:H-1);
                for g = 1:G
                    kappa2(g) = draw_ig(opts.a_kappa + (H - 1) / 2, ...
                                        opts.b_kappa + sum(e(g, :).^2) / 2);
                    kappa2(g) = min(max(kappa2(g), 1e-8), 1e4);
                end
            end
        end
        % Step-size adaptation during BURN-IN ONLY (kernel fixed after).
        % Held blocks are excluded: their acceptance is zero by
        % construction and must not shrink a step size.
        if it <= opts.n_burn && mod(it, 50) == 0
            acc = n_acc_burn ./ max(n_try_burn, 1);
            step(free_block, :) = step(free_block, :) .* ...
                exp(0.5 * (acc(free_block, :) - 0.30) / 0.30);
            step = min(max(step, 0.02), 5);
            n_acc_burn(:) = 0;  n_try_burn(:) = 0;
            accL = n_accL_burn ./ max(n_tryL_burn, 1);
            stepL(free_block) = stepL(free_block) .* ...
                exp(0.5 * (accL(free_block) - 0.30) / 0.30);
            stepL = min(max(stepL, 0.01), 5);
            n_accL_burn(:) = 0;  n_tryL_burn(:) = 0;
        end
    end

    % ---- store ----------------------------------------------------------
    if it > opts.n_burn
        keep = it - opts.n_burn;
        tau2_draws(keep, :, :) = reshape(exp(2 * x), [1, G, H]);
        x_draws(keep, :, :)    = reshape(x, [1, G, H]);
        kap_draws(keep, :)     = sqrt(kappa2)';
        beta_sum = beta_sum + beta;
        mp_sum   = mp_sum + mp_last;
        if ~isempty(opts.proj)
            rb_draws(keep, :) = opts.proj' * mp_last;
        end
        if ~isempty(opts.proj)
            proj_draws(keep, :) = opts.proj' * beta;
        end
    end
end

beta_mean = beta_sum / opts.n_keep;

% Lag-1 autocorrelation of beta, recomputed properly from the projection
% draws when available (cheap and the quantity that matters for the IRF).
lag1_beta = nan(1, H);  ess_proj = nan(1, H);  ess_rb = nan(1, H);
if ~isempty(proj_draws)
    for h = 1:H
        v = proj_draws(:, h) - mean(proj_draws(:, h));
        lag1_beta(h) = sum(v(1:end-1) .* v(2:end)) / (sum(v.^2) + eps);
        ess_proj(h)  = effective_sample_size(proj_draws(:, h));
        ess_rb(h)    = effective_sample_size(rb_draws(:, h));
    end
end

lag1_x = zeros(G, H);  ess_x = zeros(G, H);
for g = 1:G
    for h = 1:H
        v = x_draws(:, g, h);
        lag1_x(g, h) = lag1_autocorr(v);
        ess_x(g, h)  = effective_sample_size(v);
    end
end

out.beta_mean      = beta_mean;
out.beta_cond_mean = mp_last;
out.beta_rb_mean   = mp_sum / opts.n_keep;
out.rb_draws       = rb_draws;
out.sig2_mean      = sig2;
out.tau2_draws     = tau2_draws;
out.proj_draws     = proj_draws;
out.kappa_draws    = kap_draws;
out.diag.acc_rate          = n_acc ./ max(n_try, 1);
out.diag.acc_rate_level    = n_accL ./ max(n_tryL, 1);
out.diag.free_block        = free_block;    % false = tau held at 1
out.diag.lag1_acorr_logtau = lag1_x;
out.diag.ess_logtau        = ess_x;
out.diag.lag1_acorr_beta   = lag1_beta;
out.diag.ess_proj          = ess_proj;
out.diag.ess_rb            = ess_rb;
out.diag.n_tau_sweeps      = opts.n_tau_sweeps;
out.diag.n_tau_clip        = n_clip;
out.diag.mh_step           = step;
out.diag.mh_step_level     = stepL;
out.diag.all_finite = all(isfinite(beta_mean(:))) && ...
                      all(isfinite(tau2_draws(:))) && ...
                      all(isfinite(kap_draws(:)));

assert(out.diag.all_finite, 'gibbs_block_pooled_horizons: non-finite draws.');
end

% =====================================================================
function lp = site_logpost(xv, idx, xfull, p_g, S, lam2, phi, kap2, H)
% Log conditional density (up to a constant) of the sites at horizons
% idx for EVERY block, given the sites at the other parity.  xv is the
% (G x numel(idx)) matrix of candidate values; xfull supplies the fixed
% neighbours.  Vectorised over blocks and sites.
lp = -bsxfun(@times, p_g, xv) - S .* exp(-2 * xv) ./ (2 * lam2);

first = (idx == 1);
if any(first)
    % half-Cauchy(0,1) anchor on tau_{g,1}: p(x) ∝ e^x / (1 + e^{2x})
    z = xv(:, first);
    lp(:, first) = lp(:, first) + z - log1p_safe(exp(2 * z));
end
rest = ~first;
if any(rest)
    prev = xfull(:, idx(rest) - 1);
    e = xv(:, rest) - phi * prev;
    lp(:, rest) = lp(:, rest) - bsxfun(@rdivide, e.^2, 2 * kap2);
end

last = (idx == H);
notlast = ~last;
if any(notlast)
    nxt = xfull(:, idx(notlast) + 1);
    e2 = nxt - phi * xv(:, notlast);
    lp(:, notlast) = lp(:, notlast) - bsxfun(@rdivide, e2.^2, 2 * kap2);
end
end

function lp = path_logpost(xm, p_g, S, lam2, phi, kap2)
% Log density of the WHOLE log-tau path of every block (G x 1), up to an
% additive constant: every horizon's likelihood contribution, the
% half-Cauchy anchor at h = 1 and the AR(1)/random-walk increments.
% Used by the joint level move (3b).
H = size(xm, 2);
lp = sum(-bsxfun(@times, p_g, xm) - S .* exp(-2 * xm) ./ (2 * lam2), 2);
lp = lp + xm(:, 1) - log1p_safe(exp(2 * xm(:, 1)));
if H > 1
    e = xm(:, 2:H) - phi * xm(:, 1:H-1);
    lp = lp - sum(e.^2, 2) ./ (2 * kap2);
end
end

function v = log1p_safe(z)
% log(1 + z), elementwise, without overflow for very large z.
v = zeros(size(z));
big = z > 1e15;
v(big)  = log(z(big));
v(~big) = log1p(z(~big));
end

function r = lag1_autocorr(v)
v = v - mean(v);
r = sum(v(1:end-1) .* v(2:end)) / (sum(v.^2) + eps);
end

function ess = effective_sample_size(v)
% Initial-positive-sequence effective sample size (Geyer 1992), the
% standard MCMC diagnostic; returns n when the chain is uncorrelated.
n = numel(v);
v = v - mean(v);
den = sum(v.^2);
if den <= 0
    ess = n;  return;
end
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
ess = n / (1 + 2 * s);
ess = min(max(ess, 1), n);
end

function s = set_default(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
