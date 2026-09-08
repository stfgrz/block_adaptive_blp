function test_approximations()
% PURPOSE
% -------
% Pin the two structural claims that docs/APPROXIMATIONS.md rests on, so
% that they are verified rather than asserted, and check the machinery
% the sensitivity analysis uses to price them.
%
% CLAIM 1 -- the per-equation treatment of Sigma cannot move the tau = 1
% posterior mean.  With the conjugate (sigma2-scaled) prior, sigma2
% cancels from conditional (1'), so the conditional posterior mean of
% beta is sigma2-free.  Holding sigma2 at the SYSTEM NIW posterior value
% (cfg.blp.sigma2_mode = 'fixed_niw') instead of sampling it per
% equation must therefore leave the tau = 1 estimate EXACTLY unchanged.
% This is the precise sense in which "the sampler ignores cross-equation
% covariance" is harmless for the point estimate of the baseline -- and
% it stops being harmless once tau is sampled, which the test also
% shows, because the tau conditional divides by sigma2.
%
% CLAIM 2 -- the IRF is LINEAR in the impact vector.  theta_i(h) =
% b1n' * beta_block(i,:,h)'.  The sensitivity analysis relies on this to
% isolate impact-vector uncertainty EXACTLY, without re-estimating
% anything; if .beta_block ever stopped being the coefficient block that
% produces .theta_mean, that analysis would silently become wrong.
%
% Also checked: both band constructions are produced and finite, since
% the third approximation (quasi-Bayesian sandwich bands) is priced by
% comparing them.

fprintf('test_approximations:\n');

cfg = default_config();
cfg.mode = 'fmar';
cfg.fmar.h1_mode = 'lp';
cfg.H = 6;  cfg.T = 160;
cfg.gibbs.n_burn = 60;  cfg.gibbs.n_keep = 120;

rng(3101, 'twister');
dgp  = simulate_sparse_misspec_dgp(cfg);
bvar = estimate_bvar_niw(dgp.Y, cfg);
blpf = estimate_blp_fmar(dgp.Y, cfg, bvar);

% --- CLAIM 1, tau fixed: the two sigma2 treatments must agree exactly --
cfg_fix = cfg;  cfg_fix.blp.fix_tau = 1;
rng(3102, 'twister');
a = estimate_blp_blockadaptive(dgp.Y, cfg_fix, bvar, blpf.lambda);
cfg_fix2 = cfg_fix;  cfg_fix2.blp.sigma2_mode = 'fixed_niw';
rng(3102, 'twister');
b = estimate_blp_blockadaptive(dgp.Y, cfg_fix2, bvar, blpf.lambda);
d = max(max(abs(a.theta_cond - b.theta_cond)));
assert(d < 1e-12, ...
    ['with tau = 1 the posterior mean must be sigma2-free, so the two ' ...
     'sigma2 treatments must agree exactly (gap %.3g)'], d);
% and both must still equal the FMAR closed form
dF = max(max(abs(b.theta_cond(:, 2:end) - blpf.theta_mean(:, 2:end))));
assert(dF < 1e-10, 'fixed_niw broke the nesting (gap %.3g)', dF);
fprintf(['  tau = 1: sampling sigma2 per equation vs fixing it at the system\n' ...
         '           NIW value changes the estimate by %.1e (and nesting holds to %.1e)\n'], ...
        d, dF);

% --- CLAIM 1, tau sampled: the channel is real ------------------------
rng(3103, 'twister');
a2 = estimate_blp_blockadaptive(dgp.Y, cfg, bvar, blpf.lambda);
cfg2 = cfg;  cfg2.blp.sigma2_mode = 'fixed_niw';
rng(3103, 'twister');
b2 = estimate_blp_blockadaptive(dgp.Y, cfg2, bvar, blpf.lambda);
d2 = max(max(abs(a2.theta_mean(:, 3:end) - b2.theta_mean(:, 3:end))));
rel = d2 / max(max(a2.diag.post_sd(:, 2:end)));
assert(d2 > 0, ...
    'with tau sampled the sigma2 treatment should matter at least a little');
fprintf(['  tau sampled: it does matter, but by %.4f in IRF units\n' ...
         '           (%.2f of the largest posterior sd) -- the number ' ...
         'run_sensitivity_approximations reports\n'], d2, rel);

% --- CLAIM 2: theta is linear in b1n ----------------------------------
K = cfg.K;  H = cfg.H;
recon = zeros(K, H);
for h = 1:H
    recon(:, h) = a2.beta_block(:, :, h) * bvar.b1n;
end
dl = max(max(abs(recon - a2.theta_mean(:, 2:end))));
assert(dl < 1e-12, ...
    ['.beta_block does not reproduce .theta_mean (gap %.3g); the ' ...
     'impact-vector sensitivity analysis depends on this identity'], dl);
% a different impact vector must give a different, still finite, IRF
b1alt = bvar.b1n_draws(:, 1);
alt = zeros(K, H);
for h = 1:H, alt(:, h) = a2.beta_block(:, :, h) * b1alt; end
assert(all(isfinite(alt(:))) && max(max(abs(alt - recon))) > 0, ...
    'the IRF did not respond to a different impact vector');
fprintf('  theta = b1n'' * beta_block reproduces the reported IRF to %.1e: OK\n', dl);

% --- both band constructions exist and are finite ----------------------
for f = {'lo', 'hi', 'lo_post', 'hi_post'}
    v = a2.(f{1})(:, 2:end);
    assert(all(isfinite(v(:))), 'band field .%s has non-finite entries', f{1});
end
w_sand = mean(mean(a2.hi(:, 3:end) - a2.lo(:, 3:end)));
w_post = mean(mean(a2.hi_post(:, 3:end) - a2.lo_post(:, 3:end)));
assert(w_sand > 0 && w_post > 0, 'a band construction produced zero width');
fprintf(['  both band constructions present and finite (mean width: ' ...
         'sandwich %.3f, posterior %.3f)\n'], w_sand, w_post);

fprintf('PASS: test_approximations\n\n');
end
