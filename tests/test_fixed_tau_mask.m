function test_fixed_tau_mask()
% PURPOSE
% -------
% Pin the behaviour of the per-cell hold cfg.blocks.fixed_tau_mask and of
% the equation subset cfg.blp.equations, both added for the out-of-sample
% block ablation (empirical/montecarlo/run_block_ablation.m, CH7_REDESIGN
% Section 8), and smoke-test that script on the synthetic fixture.
%
% WHAT IS CHECKED, for BOTH adaptive estimators
% ---------------------------------------------
% 1. DEFAULT UNCHANGED.  mask = [] and a configuration lacking the new
%    fields give identical output.
% 2. ALL-TRUE MASK == fixed_tau = 1:G: theta_mean and tau_mean isequal
%    under the same seed, every tau summary exactly 1.
% 3. ONE CELL HELD.  That cell reports tau_mean = tau_med = tau_q = 1 and
%    P(tau > 1) = 0 at every h; the same block in the OTHER equations and
%    the other blocks of the SAME equation still move; the mask is
%    recorded in blp.fixed_tau_mask.
% 4. UNION with fixed_tau.  mask(2,1) plus fixed_tau = 3 holds column 3
%    and cell (2,1), and is isequal to passing the union mask directly.
% 5. EQUATION SUBSET.  cfg.blp.equations = 2 leaves equations 1 and 3 NaN
%    (tau included, i.e. NaN not 1), equation 2 finite, h = 0 and
%    prior_theta filled for all, tau summaries with the full shape.
%    (Skipped equations consume no draws, so equation 2 is NOT draw-
%    identical to the full run; that is documented, not tested.)
% 6. VALIDATION.  A wrong-sized mask, a non-logical mask and an out-of-
%    range equation index are refused.
% 7. SMOKE.  run_block_ablation on the synthetic fixture in quick mode:
%    shapes and finiteness only (the numbers mean nothing by design).

fprintf('test_fixed_tau_mask:\n');

cfg = default_config();
cfg.mode = 'fmar';
cfg.fmar.h1_mode = 'lp';
cfg.H = 8;  cfg.T = 180;
cfg.gibbs.n_burn = 80;  cfg.gibbs.n_keep = 160;
cfg.fmar.n_niw_draws = 40;

rng(6201, 'twister');
dgp  = simulate_sparse_misspec_dgp(cfg);
bvar = estimate_bvar_niw(dgp.Y, cfg);
blpf = estimate_blp_fmar(dgp.Y, cfg, bvar);
K = cfg.K;  G = K;  H = cfg.H;
nQ = numel(cfg.blp.tau_probs);

ests = {'independent', @(c, s) run_ind(dgp.Y, c, bvar, blpf.lambda, s); ...
        'pooled',      @(c, s) run_pool(dgp.Y, c, bvar, blpf.lambda, s)};
for e = 1:size(ests, 1)
    nm = ests{e, 1};  f = ests{e, 2};

    % --- 1. default unchanged -------------------------------------------
    c_nofield = cfg;
    c_nofield.blocks = rmfield(c_nofield.blocks, 'fixed_tau_mask');
    c_nofield.blp    = rmfield(c_nofield.blp, 'equations');
    c_empty = cfg;  c_empty.blocks.fixed_tau_mask = [];  c_empty.blp.equations = [];
    a = f(c_nofield, 6202);  b = f(c_empty, 6202);
    assert(isequal(a.theta_mean, b.theta_mean) && isequal(a.tau_mean, b.tau_mean) ...
        && isequaln(a.lo, b.lo) && isequal(a.beta_mean, b.beta_mean), ...   % lo(:,1) is NaN by design
        '%s: an empty fixed_tau_mask / equations changed the output', nm);
    assert(isequal(a.fixed_tau_mask, false(K, G)) && isequal(a.equations, 1:K), ...
        '%s: default fixed_tau_mask / equations not reported as all-false / 1:K', nm);

    % --- 2. all-true mask == fixed_tau = 1:G ----------------------------
    c_all  = cfg;  c_all.blocks.fixed_tau = 1:G;
    c_mask = cfg;  c_mask.blocks.fixed_tau_mask = true(K, G);
    r_all = f(c_all, 6203);  r_mask = f(c_mask, 6203);
    assert(isequal(r_all.theta_mean, r_mask.theta_mean) && isequal(r_all.tau_mean, r_mask.tau_mean), ...
        '%s: all-true mask differs from fixed_tau = 1:G', nm);
    assert(all(r_mask.tau_mean(:) == 1) && all(r_mask.tau_med(:) == 1) && ...
           all(r_mask.tau_q(:) == 1) && all(r_mask.p_tau_gt1(:) == 0), ...
        '%s: all-true mask does not report tau = 1 everywhere', nm);
    d = max(max(abs(r_mask.theta_cond(:, 2:end) - blpf.theta_mean(:, 2:end))));
    assert(d < 1e-10, '%s: all-true mask is not the FMAR closed form (%.2e)', nm, d);

    % --- 3. one cell held --------------------------------------------------
    M1 = false(K, G);  M1(2, 1) = true;
    c_one = cfg;  c_one.blocks.fixed_tau_mask = M1;
    r_one = f(c_one, 6204);
    assert(all(r_one.tau_mean(2, 1, :) == 1) && all(r_one.tau_med(2, 1, :) == 1) && ...
           all(all(r_one.tau_q(2, 1, :, :) == 1)) && all(r_one.p_tau_gt1(2, 1, :) == 0), ...
        '%s: the held cell (2,1) moved', nm);
    same_block_other_eq = r_one.tau_mean([1 3], 1, :);
    assert(any(same_block_other_eq(:) ~= 1), '%s: block 1 stopped moving in the other equations', nm);
    same_eq_other_blocks = r_one.tau_mean(2, 2:end, :);
    assert(any(same_eq_other_blocks(:) ~= 1), '%s: the other blocks of equation 2 stopped moving', nm);
    assert(isequal(r_one.fixed_tau_mask, M1) && isempty(r_one.fixed_tau_blocks), ...
        '%s: fixed_tau_mask not recorded', nm);
    d1 = max(abs(r_one.theta_mean(2, 2:end) - a.theta_mean(2, 2:end)));
    assert(d1 > 1e-8, '%s: holding cell (2,1) left equation 2 unchanged (%.1e)', nm, d1);

    % --- 4. union with fixed_tau ---------------------------------------------
    c_un = cfg;  c_un.blocks.fixed_tau_mask = M1;  c_un.blocks.fixed_tau = 3;
    Mu = M1;  Mu(:, 3) = true;
    c_mu = cfg;  c_mu.blocks.fixed_tau_mask = Mu;
    r_un = f(c_un, 6205);  r_mu = f(c_mu, 6205);
    assert(isequal(r_un.fixed_tau_mask, Mu) && isequal(r_un.fixed_tau_blocks, 3), ...
        '%s: union mask not recorded', nm);
    col3 = r_un.tau_mean(:, 3, :);  c21 = r_un.tau_mean(2, 1, :);
    assert(all(col3(:) == 1) && all(c21(:) == 1), ...
        '%s: union did not hold column 3 and cell (2,1)', nm);
    free_un = r_un.tau_mean(:, 2, :);
    assert(any(free_un(:) ~= 1), '%s: union held a block it should not have', nm);
    assert(isequal(r_un.theta_mean, r_mu.theta_mean) && isequal(r_un.tau_mean, r_mu.tau_mean), ...
        '%s: union of mask and fixed_tau differs from the explicit union mask', nm);

    % --- 5. equation subset ----------------------------------------------------
    c_eq = cfg;  c_eq.blp.equations = 2;
    r_eq = f(c_eq, 6206);
    assert(isequal(r_eq.equations, 2), '%s: equations not recorded', nm);
    for fld = {'theta_mean', 'theta_med', 'theta_cond', 'lo', 'hi', 'lo_post', 'hi_post'}
        v = r_eq.(fld{1});
        assert(isequal(size(v), [K, H + 1]), '%s: %s has the wrong shape', nm, fld{1});
        assert(all(isnan(reshape(v([1 3], 2:end), [], 1))), '%s: %s not NaN for skipped equations', nm, fld{1});
        assert(all(isfinite(v(2, 2:end))), '%s: %s not finite where estimated', nm, fld{1});
        if strncmp(fld{1}, 'theta', 5)      % bands carry no h = 0 entry (NaN by design)
            assert(all(isfinite(v(:, 1))), '%s: %s h = 0 column not filled for all equations', nm, fld{1});
        end
    end
    assert(all(isfinite(r_eq.prior_theta(:))), '%s: prior_theta not filled for all equations', nm);
    assert(isequal(size(r_eq.tau_mean), [K, G, H]) && isequal(size(r_eq.tau_q), [K, G, H, nQ]) && ...
           isequal(size(r_eq.p_tau_gt1), [K, G, H]), '%s: tau summaries have the wrong shape', nm);
    tsk = r_eq.tau_mean([1 3], :, :);  psk = r_eq.p_tau_gt1([1 3], :, :);  qsk = r_eq.tau_q([1 3], :, :, :);
    assert(all(isnan(tsk(:))) && all(isnan(psk(:))) && all(isnan(qsk(:))), ...
        '%s: tau summaries of skipped equations are not NaN', nm);
    t2 = r_eq.tau_mean(2, :, :);  q2 = r_eq.tau_q(2, :, :, :);
    assert(all(isfinite(t2(:))) && all(isfinite(q2(:))) && all(t2(:) > 0), ...
        '%s: tau summaries of equation 2 not finite/positive', nm);
    bsk = r_eq.beta_mean(:, [1 3], :);  b2 = r_eq.beta_mean(:, 2, :);
    assert(all(isnan(bsk(:))) && all(isfinite(b2(:))), '%s: beta_mean subset pattern wrong', nm);
    esk = r_eq.diag.ess_rb([1 3], :);
    assert(all(isnan(esk(:))) && all(isfinite(r_eq.diag.ess_rb(2, :))), ...
        '%s: diag entries subset pattern wrong', nm);
    if strcmp(nm, 'independent')
        bbs = r_eq.beta_block([1 3], :, :);  lsk = r_eq.lambda([1 3], :);
        assert(all(isnan(bbs(:))) && all(isnan(lsk(:))), ...
            '%s: beta_block / lambda subset pattern wrong', nm);
    else
        ksk = r_eq.kappa_mean([1 3], :);
        assert(all(isnan(ksk(:))) && all(isfinite(r_eq.kappa_mean(2, :))), ...
            '%s: kappa_mean subset pattern wrong', nm);
    end
    % subset combined with a held cell in that equation
    c_eqm = c_eq;  c_eqm.blocks.fixed_tau_mask = M1;
    r_eqm = f(c_eqm, 6206);
    h21 = r_eqm.tau_mean(2, 1, :);  s11 = r_eqm.tau_mean(1, 1, :);
    assert(all(h21(:) == 1) && all(isnan(s11(:))), ...
        '%s: subset + mask: held cell not 1 or skipped cell not NaN', nm);
    fprintf('  %-11s: default unchanged; all-true = fixed_tau (FMAR to %.1e); one cell; union; subset: OK\n', nm, d);
end

% --- 6. validation ------------------------------------------------------------
bad = {@() setmask(cfg, true(K, G + 1)), 'wrong-sized mask'; ...
       @() setmask(cfg, 2 * ones(K, G)),  'non-logical mask'; ...
       @() seteq(cfg, K + 1),             'out-of-range equation'; ...
       @() seteq(cfg, 0),                 'zero equation index'};
for r = 1:size(bad, 1)
    for e = 1:size(ests, 1)
        threw = false;
        try
            ests{e, 2}(bad{r, 1}(), 1);
        catch
            threw = true;
        end
        assert(threw, '%s: %s was accepted', ests{e, 1}, bad{r, 2});
    end
end
fprintf('  wrong-sized / non-logical mask and bad equation index refused: OK\n');

% --- 7. smoke run of the block ablation on the synthetic fixture -------------
ds = make_synthetic_fixture(struct('T', 120, 'seed', 6207, 'out_mat', ''));
out_dir = fullfile(tempname, 'ablation_smoke');
o = struct('dataset', ds, 'cells', {{'hicp', 'hicp'}}, 'control_cells', {{'ip', 'stoxx'}}, ...
           'quick', true, 'allow_synthetic', true, 'seed', 6208, ...
           'out_dir', out_dir, 'out_stem', 'smoke', 'verbose', false);
res = run_block_ablation(o);
nO = numel(res.origins);  Hq = res.cfg.H;  Kq = size(ds.Y, 2);
assert(nO == 3 && isequal(res.h_eval, [2 4]) && Hq == 4 && res.cfg.p == 2, ...
    'ablation smoke: quick settings not applied');
assert(res.synthetic, 'ablation smoke: synthetic flag not carried');
assert(isequal(size(res.err.free), [nO, Hq, 2]) && isequal(size(res.err.ablated), [nO, Hq, 2]) && ...
       isequal(size(res.err.fmar), [nO, Hq, Kq]) && isequal(size(res.err.bvar), [nO, Hq, Kq]), ...
    'ablation smoke: error arrays have the wrong shape');
assert(all(isfinite(res.err.free(:))) && all(isfinite(res.err.ablated(:))) && ...
       all(isfinite(res.err.fmar(:))) && all(isfinite(res.err.bvar(:))), ...
    'ablation smoke: non-finite forecast errors');
assert(all(isfinite(res.tau_free(:))) && all(res.p_gt1_free(:) >= 0 & res.p_gt1_free(:) <= 1), ...
    'ablation smoke: tau summaries of the free run not finite');
assert(any(res.err.free(:) ~= res.err.ablated(:)), 'ablation smoke: ablation changed nothing');
assert(numel(res.table) == 2 * (numel(res.h_eval) + 1), 'ablation smoke: wrong number of table rows');
assert(isequal({res.table.role}, {'escape', 'escape', 'escape', 'control', 'control', 'control'}), ...
    'ablation smoke: roles wrong');
for r = 1:numel(res.table)
    tr = res.table(r);
    assert(tr.n == nO && isfinite(tr.msfe_free) && isfinite(tr.msfe_ablated) && ...
           isfinite(tr.ratio) && isfinite(tr.dm_stat) && isfinite(tr.cw_stat), ...
        'ablation smoke: non-finite statistic in row %d (%s, h = %d)', r, tr.cell, tr.h);
end
assert(exist(res.files.mat, 'file') == 2 && exist(res.files.csv, 'file') == 2, ...
    'ablation smoke: outputs not written');
% a refused synthetic dataset
threw = false;
try
    o2 = o;  o2.allow_synthetic = false;
    run_block_ablation(o2);
catch
    threw = true;
end
assert(threw, 'ablation smoke: synthetic dataset accepted without allow_synthetic');
fprintf('  run_block_ablation quick smoke (3 origins, 2 cells): shapes, finiteness, synthetic guard: OK\n');

fprintf('PASS: test_fixed_tau_mask\n\n');
end

% =====================================================================
function r = run_ind(Y, cfg, bvar, lam, seed)
rng(seed, 'twister');
r = estimate_blp_blockadaptive(Y, cfg, bvar, lam);
end

function r = run_pool(Y, cfg, bvar, lam, seed)
rng(seed, 'twister');
r = estimate_blp_blockpooled(Y, cfg, bvar, lam);
end

function c = setmask(c, M)
c.blocks.fixed_tau_mask = M;
end

function c = seteq(c, v)
c.blp.equations = v;
end
