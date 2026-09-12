% RUN_EMPIRICAL  Chapter 7 driver: block-adaptive BLP on euro-area data.
%
% Design document: docs/CH7_DESIGN.md.  Integration notes, data status and
% verified repo signatures: README_EMPIRICAL.md.
%
% PREREQUISITES (run once, in this order; each is checked below and the
% script stops with an actionable message if one is missing)
%   1. build_shock_series();        (needs data/raw/ea_empd_events.csv,
%                                    shipped with the package)
%   2. fetch_outcome_data();        (must report "all four outcome files
%                                    ready and validated"; if the downloads
%                                    fail it prints the manual route)
%   3. assemble_dataset();          (baseline variant mps_gc_1y)
%   4. SMOKE_TEST_EMPIRICAL         (seconds; verifies the estimator
%                                    signatures before you commit hours)
% Step D builds its own JK-variant dataset if it is not already there, so
% assemble_dataset(struct('shock_variant','mps_gc_jk')) is optional.
%
% STEPS (toggles below; each can be preset in the workspace before running,
% e.g.  PSI_FLOOR = false; DO_C = false; RUN_EMPIRICAL )
%   A  baseline estimation at p = 12: BVAR, BLP-FMAR, BLP-block,
%      BLP-pooled, LP-LA, plus the instrument-relevance diagnostics
%   B  tau heatmap export for both adaptive estimators  [needs A's output]
%   C  dose-response in p in {2, 6, 12} (Leg 3; see CH7_DESIGN Sec. 5.2)
%   D  JK information-effect comparison                 [needs A's output]
%   E  coherence table (Sec. 5.3) and, when the matching null calibration
%      exists, the pre-registered READING PROTOCOL of Sec. 5.1 applied
%      mechanically to the step-A tau maps            [needs A's output]
%   F  the four-variable LEVEL system without the surprise (i1y, ip, hicp,
%      stoxx; recursive innovation to i1y): the same diagnostic on a
%      conventional monetary VAR, as a companion exhibit and as the check
%      that the surprise block is what drives the tightness path (below)
%   The null calibration itself (empirical/montecarlo/run_null_calibration.m)
%   is a SEPARATE, longer job; run it after A-D, then rerun with only
%   DO_E = true to produce the protocol table.
%
% ESTIMATORS.  Both block-adaptive estimators of the simulation study are
% run: the INDEPENDENT one (tau_{i,g,h} free at every horizon,
% estimate_blp_blockadaptive) and the HORIZON-POOLED one (log tau smoothed
% across h, estimate_blp_blockpooled), which the Monte Carlo found to give
% a markedly sharper tau map (docs/RESULTS.md).  DO_POOLED switches the
% second off.
%
% h = 1 CONVENTION.  cfg.fmar.h1_mode = 'lp', as in every reported Monte
% Carlo (README Sec. 3a): the global baseline estimates h = 1 as a local
% projection with its own marginal-likelihood tightness, so the adaptive
% estimators inherit an LP tightness at h = 1 rather than the BVAR's
% Minnesota lambda.  The BVAR is still reported separately at every h.
%
% PRIOR-SCALE FLOOR (PSI_FLOOR, default true).  The published FMAR prior
% scale psi_v(h) is a Newey-West long-run variance.  For the policy
% surprise -- a mean-reverting, near-white-noise series -- that estimate
% falls to a fifth of the residual variance at long horizons, which makes
% the Minnesota-scaled prior on every coefficient attached to lagged
% surprises far too loose; the marginal likelihood then compensates by
% driving the GLOBAL lambda_h to ~0.01 for h <= 28, and its objective
% becomes bimodal (modes near 0.01 and 0.4 of almost equal height) from
% h ~ 27, so the selected lambda_h flips between horizons and the tau
% surface shows a vertical stripe.  With the floor (psi never below the
% residual variance; priors/fmar_prior_scale.m) the objective is unimodal
% at every horizon and lambda_h rises smoothly from 0.01 to ~0.4.  The
% floor never binds for the four persistent series, and on the simulation
% designs it binds by a few percent at h = 1 only (max IRF change 1e-3,
% tests/test_psi_floor.m).  Set PSI_FLOOR = false
% to reproduce the pure published scale; every output then carries the
% suffix _nofloor so the two runs coexist.
%
% THE INSTRUMENT BLOCK IS HELD AT tau = 1 (INSTR_FIXED, default true).
% The lag block of the surprise holds one large coefficient -- the
% contemporaneous one, which carries the identification of every IRF --
% and eleven near-zero ones, so it violates the group prior's common-scale
% assumption: with every block free, the surprise block's tau collapsed to
% 0.03-0.4 in every equation and dragged the identifying coefficient to
% the VAR centre, which put the adaptive IRFs on top of the BVAR while the
% global BLP followed the LP.  With cfg.blocks.fixed_tau = 1 the surprise
% block stays under the global FMAR prior, exactly like the intercept, and
% the diagnostic is read on the four LEVEL blocks of each equation.  Set
% INSTR_FIXED = false to let every block adapt; outputs then carry the
% suffix _allblocks.
%
% RESTARTING: every step saves its own .mat into results/, so if the run
% dies part-way you can set the completed steps' DO_ toggles to false and
% rerun.  B, D and E read results/empirical_<tag>.mat, which step A writes;
% if you switch A off before it has ever run, they stop with a clear
% message instead of a cryptic load error.
%
% RUNTIME: about five minutes in MATLAB for A-F at p = 12, H = 48 (a few
% hours in Octave).  Plotting is wrapped in try/catch, so a headless
% session loses the pngs but keeps every .mat and .csv -- the real outputs.

if ~exist('DO_A', 'var'), DO_A = true; end
if ~exist('DO_B', 'var'), DO_B = true; end
if ~exist('DO_C', 'var'), DO_C = true; end
if ~exist('DO_D', 'var'), DO_D = true; end
if ~exist('DO_E', 'var'), DO_E = true; end
if ~exist('DO_F', 'var'), DO_F = true; end
if ~exist('DO_POOLED', 'var'), DO_POOLED = true; end
if ~exist('PSI_FLOOR', 'var'), PSI_FLOOR = true; end
if ~exist('INSTR_FIXED', 'var'), INSTR_FIXED = true; end

% --- path bootstrap (works from any directory, and when pasted) ----------
ea_this = mfilename('fullpath');
if isempty(ea_this), ea_this = fullfile(pwd, 'RUN_EMPIRICAL'); end
addpath(genpath(fileparts(fileparts(ea_this))));
P = ea_paths();
if exist(P.results, 'dir') ~= 7, mkdir(P.results); end

% --- configuration -----------------------------------------------------------
assert(exist(P.dataset, 'file') == 2, ...
       ['RUN_EMPIRICAL: %s not found.\n' ...
        'Run the prerequisites first:\n' ...
        '    build_shock_series(); fetch_outcome_data(); assemble_dataset();\n' ...
        'then SMOKE_TEST_EMPIRICAL before launching this.'], P.dataset);
ds  = load(P.dataset);

% --- refuse a SYNTHETIC dataset -----------------------------------------
% empirical/tests/make_synthetic_fixture.m can build a simulated dataset
% with the right shape so the interfaces can be exercised on a machine
% with no euro-area data.  Nothing computed from it is an empirical
% result, so it must never reach this driver: the fixture stamps
% ds.synthetic = true and this is the check that stops it.
if isfield(ds, 'synthetic') && ds.synthetic
    error(['%s: this dataset is a SYNTHETIC FIXTURE (ds.synthetic = true), ' ...
           'built for interface testing only. It is simulated data, not ' ...
           'euro-area data, and no result computed from it is an empirical ' ...
           'finding. Build the real dataset first:\n' ...
           '    build_shock_series(); fetch_outcome_data(); assemble_dataset();'], ...
          mfilename);
end

cfg = default_config();
cfg.mode      = 'fmar';
cfg.fmar.h1_mode = 'lp';           % see the header: h = 1 is an LP, as in the MC
cfg.fmar.psi_floor = logical(PSI_FLOOR);   % see the header
if INSTR_FIXED
    cfg.blocks.fixed_tau = 1;      % the surprise block stays under the global prior
else
    cfg.blocks.fixed_tau = [];
end
cfg.shock_var = 1;                 % surprise ordered first (internal instrument)
cfg.p         = 12;                % monthly baseline lag order
cfg.H         = 48;                % 4-year horizon
cfg.ci_level  = 0.90;
cfg.fmar.isrw = ds.isrw;           % [0 1 1 1 1]; needs the vectorised isrw
                                   % prior mean in estimate_bvar_niw
K = size(ds.Y, 2);
cfg.K = K;                         % keep cfg self-consistent in saved output
                                   % (the estimators read K from size(Y,2);
                                   %  default_config ships the K = 3 of the
                                   %  simulation DGPs)

% Early-horizon window of the tau_bar statistic (CH7_DESIGN Sec. 5.1).
% CLAMPED to cfg.H exactly as run_null_calibration clamps it, so the two
% always compute the same statistic and a reduced-H run cannot die
% part-way through with an out-of-bound index.
h_early = [2, min(12, cfg.H)];
assert(h_early(2) >= h_early(1), ...
       'RUN_EMPIRICAL: cfg.H = %d is too small for the early-h window.', cfg.H);

% Output filenames carry the baseline lag order and the prior-scale
% convention, so runs at another p, or with the published scale, cannot
% overwrite the baseline the write-up cites.
tag = sprintf('p%d', cfg.p);
if ~PSI_FLOOR,   tag = [tag '_nofloor']; end
if ~INSTR_FIXED, tag = [tag '_allblocks']; end
step_a_file = fullfile(P.results, sprintf('empirical_%s.mat', tag));
% the level system (step F) has no instrument: its files carry only the
% lag order and the prior-scale convention
tag_lev = sprintf('p%d', cfg.p);
if ~PSI_FLOOR, tag_lev = [tag_lev '_nofloor']; end

% What the dataset was built from travels with every result file.
ds_info = struct('varnames', {ds.varnames}, 'ym', ds.ym, ...
                 'shock_variant', ds.shock_variant, 'meta', ds.meta);

fprintf('RUN_EMPIRICAL: %s\n  T = %d, K = %d, p = %d, H = %d, early-h window [%d %d], variant = %s\n  h1_mode = %s, psi_floor = %d, instrument block held at tau = 1: %d, pooled = %d, tag = %s\n', ...
        P.dataset, size(ds.Y, 1), K, cfg.p, cfg.H, h_early(1), h_early(2), ...
        ds.shock_variant, cfg.fmar.h1_mode, cfg.fmar.psi_floor, INSTR_FIXED, DO_POOLED, tag);

% ==============================================================================
if DO_A
fprintf('\n[A] Baseline estimation, p = %d, H = %d ...\n', cfg.p, cfg.H);
t_a = tic;
rng(20260107, 'twister');
bvar = estimate_bvar_niw(ds.Y, cfg);
fprintf('  BVAR done (%.1f s): lambda = %.3f, max |eig| = %.4f\n', toc(t_a), bvar.lambda, bvar.max_eig);
blpf = estimate_blp_fmar(ds.Y, cfg, bvar);
fprintf('  BLP-FMAR done (%.1f s)\n', toc(t_a));
fprintf('  lambda_h, h = 1..%d:\n   ', cfg.H);
fprintf(' %.3f', blpf.lambda(1, :));  fprintf('\n');
fprintf('  lambda objective multimodal at %d of %d horizons; psi floor bound at %d (variable, horizon) cells\n', ...
        blpf.diag.lambda_multimodal, cfg.H, sum(blpf.diag.psi_floor_bound));
if blpf.diag.lambda_multimodal > 0
    fprintf('  horizons with a multimodal lambda objective:');
    fprintf(' %d', find(blpf.diag.n_local_max > 1));  fprintf('\n');
end
blpb = estimate_blp_blockadaptive(ds.Y, cfg, bvar, blpf.lambda);
fprintf('  BLP-block done (%.1f s)\n', toc(t_a));
if DO_POOLED
    blpp = estimate_blp_blockpooled(ds.Y, cfg, bvar, blpf.lambda);
    fprintf('  BLP-pooled done (%.1f s)\n', toc(t_a));
else
    blpp = [];
end
lpla = estimate_lp_lagaug(ds.Y, cfg, bvar);

% Instrument relevance and the 25 bp normalisation: common factor from the
% BVAR impact response of the 1Y rate (variable 2) to the surprise
% innovation.  Applied at REPORTING time only; all saved IRFs are in
% unit-surprise-innovation scale.
rel = ea_relevance(ds.Y, bvar, cfg);
imp_i1y = bvar.theta(2, 1);
k25 = rel.k25;
fprintf('  relevance: impact of a 1pp surprise innovation on i1y = %.4f pp (robust t = %.2f, F = %.2f);\n', ...
        rel.impact_b, rel.impact_t, rel.impact_F);
fprintf('             naive same-month d i1y on mps: b = %.3f (t = %.2f); 25bp scale k = %.3f\n', ...
        rel.naive_b, rel.naive_t, k25);
% The normalisation is only meaningful if a positive surprise raises the
% 1Y rate on impact by a non-trivial amount.  A near-zero or negative
% impact response silently turns every reported IRF into nonsense (huge
% or sign-flipped), so say so loudly rather than plotting it.
if ~(imp_i1y > 0)
    warning('RUN_EMPIRICAL:normalisation', ...
            ['impact response of i1y to the surprise is %.4g (<= 0). The 25bp ' ...
             'normalisation FLIPS THE SIGN of every reported IRF. Check the ' ...
             'shock series and the variable ordering before reading step A.'], imp_i1y);
elseif imp_i1y < 0.01
    warning('RUN_EMPIRICAL:normalisation', ...
            ['impact response of i1y to the surprise is only %.4g pp, so the ' ...
             '25bp scale factor is %.1f. Reported IRFs are amplified by that ' ...
             'factor; treat them with care.'], imp_i1y, k25);
end
if rel.impact_F < 10
    fprintf(['  NOTE: first-stage F = %.2f < 10 -- the surprise is a WEAK instrument for the\n' ...
             '        monthly-average 1Y rate.  The tau diagnostic is a statement about the\n' ...
             '        dynamic specification of the system and does not depend on this; the\n' ...
             '        structural reading of the IRFs and the 25bp scale do.  Report it.\n'], rel.impact_F);
end
if bvar.max_eig >= 1
    fprintf('  note: fitted BVAR has max |eig| = %.4f >= 1 (expected with RW\n', bvar.max_eig);
    fprintf('        prior means on levels; the null bootstrap conditions on\n');
    fprintf('        the actual first p observations, see simulate_fitted_bvar_dgp).\n');
end

% The NIW companion draws (Kp x Kp x n_draws, ~14 MB at p = 12) are only
% used by the simulation-side sensitivity analysis; drop them from the
% saved copy so the result file stays small enough to keep in the repo.
bvar = slim_bvar(bvar);
save(step_a_file, 'bvar', 'blpf', 'blpb', 'blpp', 'lpla', 'cfg', 'k25', 'rel', ...
     'h_early', 'ds_info', 'tag');
fprintf('  step A done in %.1f min -> %s\n', toc(t_a) / 60, step_a_file);

% quick-look IRF figure (one panel per variable, 25bp-normalised)
try
    figure('visible', 'off', 'Position', [50 50 1100 800]);
    hs = 0:cfg.H;
    for i = 1:K
        subplot(ceil(K / 2), 2, i);  hold on
        plot(hs, k25 * bvar.theta(i, :), 'r-');
        plot(hs, k25 * blpf.theta_mean(i, :), 'b-');
        plot(hs, k25 * blpb.theta_mean(i, :), 'k-');
        if DO_POOLED, plot(hs, k25 * blpp.theta_mean(i, :), 'm-'); end
        plot(hs, k25 * lpla.theta(i, :), 'g--');
        plot(hs, k25 * blpb.lo(i, :), 'k:');  plot(hs, k25 * blpb.hi(i, :), 'k:');
        plot(hs, zeros(size(hs)), '-', 'Color', [0.6 0.6 0.6]);
        title(sprintf('%s (25bp surprise, k = %.2f)', ds.varnames{i}, k25));
        if i == 1
            if DO_POOLED
                legend('BVAR', 'BLP-FMAR', 'BLP-block', 'BLP-pooled', 'LP-LA', 'Location', 'best');
            else
                legend('BVAR', 'BLP-FMAR', 'BLP-block', 'LP-LA', 'Location', 'best');
            end
        end
    end
    print(fullfile(P.results, sprintf('fig_irf_%s.png', tag)), '-dpng', '-r110');
    % the tightness path, the object behind every tau in the heatmaps
    figure('visible', 'off', 'Position', [50 50 700 380]);
    semilogy(1:cfg.H, blpf.lambda(1, :), 'k.-');  grid on
    xlabel('horizon h');  ylabel('\lambda_h (log scale)');
    title(sprintf('selected tightness path (%s); multimodal at %d horizons', ...
                  tag, blpf.diag.lambda_multimodal), 'Interpreter', 'none');
    print(fullfile(P.results, sprintf('fig_lambda_%s.png', tag)), '-dpng', '-r110');
catch err
    fprintf('  (plotting skipped: %s)\n', err.message);
end
end

% ==============================================================================
if DO_B
fprintf('\n[B] Tau heatmap export ...\n');
assert(exist(step_a_file, 'file') == 2, ...
       ['RUN_EMPIRICAL step B needs %s, which step A writes.\n' ...
        'Set DO_A = true (once), or rerun with DO_B = false.'], step_a_file);
L = load(step_a_file);
ests = {'block', L.blpb};
if isfield(L, 'blpp') && ~isempty(L.blpp), ests(end + 1, :) = {'pooled', L.blpp}; end
out_csv = fullfile(P.results, sprintf('tau_heatmap_%s.csv', tag));
fid = fopen(out_csv, 'w');
assert(fid > 0, 'RUN_EMPIRICAL: cannot write %s', out_csv);
fprintf(fid, 'scale_estimator,equation,block,block_held_at_1,h,lambda_h,tau_mean,tau_median,p_tau_gt1\n');
for e = 1:size(ests, 1)
    est = ests{e, 2};
    held = false(1, K);
    if isfield(est, 'fixed_tau_blocks'), held(est.fixed_tau_blocks) = true; end
    for i = 1:K
        for g = 1:K
            for h = 1:size(est.tau_mean, 3)
                fprintf(fid, '%s,%s,%s,%d,%d,%.6g,%.5f,%.5f,%.4f\n', ests{e, 1}, ...
                        ds.varnames{i}, ds.varnames{g}, held(g), h, L.blpf.lambda(1, h), ...
                        est.tau_mean(i, g, h), est.tau_med(i, g, h), est.p_tau_gt1(i, g, h));
            end
        end
    end
end
fclose(fid);
try
    for e = 1:size(ests, 1)
        est = ests{e, 2};
        figure('visible', 'off', 'Position', [50 50 1100 800]);
        for i = 1:K
            subplot(ceil(K / 2), 2, i);
            imagesc(log(squeeze(est.tau_mean(i, :, :))));  colorbar
            title(sprintf('eq %s: log tau (%s) by block x h', ds.varnames{i}, ests{e, 1}));
            set(gca, 'YTick', 1:K, 'YTickLabel', ds.varnames);
            xlabel('horizon h');
        end
        print(fullfile(P.results, sprintf('fig_tau_heatmap_%s_%s.png', ests{e, 1}, tag)), '-dpng', '-r110');
    end
catch err
    fprintf('  (plotting skipped: %s)\n', err.message);
end
fprintf('  wrote %s\n', out_csv);
end

% ==============================================================================
if DO_C
fprintf('\n[C] Dose-response in p (Leg 3) ...\n');
plist = [2 6 12];
h1 = h_early(1);  h2 = h_early(2);                 % protocol window, clamped
dose = struct('p', plist, 'h_early', h_early, 'H', cfg.H, 'varnames', {ds.varnames}, ...
              'tag', tag, 'psi_floor', cfg.fmar.psi_floor, ...
              'tau_bar_block', nan(numel(plist), K, K), ...
              'tau_bar_pooled', nan(numel(plist), K, K), ...
              'pgt1_bar_block', nan(numel(plist), K, K), ...
              'pgt1_bar_pooled', nan(numel(plist), K, K), ...
              'tau_block', {cell(1, numel(plist))}, 'tau_pooled', {cell(1, numel(plist))}, ...
              'lambda', {cell(1, numel(plist))}, ...
              'theta_bvar', {cell(1, numel(plist))}, 'theta_fmar', {cell(1, numel(plist))}, ...
              'theta_block', {cell(1, numel(plist))}, 'theta_pooled', {cell(1, numel(plist))}, ...
              'bvar_lambda', nan(1, numel(plist)), 'bvar_max_eig', nan(1, numel(plist)), ...
              'lambda_multimodal', nan(1, numel(plist)), 'rel', {cell(1, numel(plist))});
for ip = 1:numel(plist)
    cfp = cfg;  cfp.p = plist(ip);
    rng(20260200 + plist(ip), 'twister');
    t_c = tic;
    bv = estimate_bvar_niw(ds.Y, cfp);
    bf = estimate_blp_fmar(ds.Y, cfp, bv);
    bb = estimate_blp_blockadaptive(ds.Y, cfp, bv, bf.lambda);
    dose.tau_bar_block(ip, :, :)  = mean(bb.tau_mean(:, :, h1:h2), 3);
    dose.pgt1_bar_block(ip, :, :) = mean(bb.p_tau_gt1(:, :, h1:h2), 3);
    dose.tau_block{ip} = bb.tau_mean;
    dose.lambda{ip} = bf.lambda(1, :);
    dose.lambda_multimodal(ip) = bf.diag.lambda_multimodal;
    dose.theta_bvar{ip} = bv.theta;  dose.theta_fmar{ip} = bf.theta_mean;
    dose.theta_block{ip} = bb.theta_mean;
    if DO_POOLED
        bp = estimate_blp_blockpooled(ds.Y, cfp, bv, bf.lambda);
        dose.tau_bar_pooled(ip, :, :)  = mean(bp.tau_mean(:, :, h1:h2), 3);
        dose.pgt1_bar_pooled(ip, :, :) = mean(bp.p_tau_gt1(:, :, h1:h2), 3);
        dose.tau_pooled{ip} = bp.tau_mean;
        dose.theta_pooled{ip} = bp.theta_mean;
    end
    dose.bvar_lambda(ip) = bv.lambda;  dose.bvar_max_eig(ip) = bv.max_eig;
    dose.rel{ip} = ea_relevance(ds.Y, bv, cfp);
    fprintf('  p = %2d done in %.1f min (BVAR lambda %.3f, max|eig| %.4f, impact t %.2f, lambda multimodal at %d h)\n', ...
            plist(ip), toc(t_c) / 60, bv.lambda, bv.max_eig, dose.rel{ip}.impact_t, ...
            bf.diag.lambda_multimodal);
end
save(fullfile(P.results, sprintf('dose_response_%s.mat', tag)), '-struct', 'dose');
out_csv = fullfile(P.results, sprintf('dose_response_%s.csv', tag));
fid = fopen(out_csv, 'w');
assert(fid > 0, 'RUN_EMPIRICAL: cannot write %s', out_csv);
fprintf(fid, 'p,equation,block,tau_bar_early_block,pgt1_early_block,tau_bar_early_pooled,pgt1_early_pooled\n');
for ip = 1:numel(plist)
    for i = 1:K
        for g = 1:K
            fprintf(fid, '%d,%s,%s,%.5f,%.4f,%.5f,%.4f\n', plist(ip), ds.varnames{i}, ...
                    ds.varnames{g}, dose.tau_bar_block(ip, i, g), dose.pgt1_bar_block(ip, i, g), ...
                    dose.tau_bar_pooled(ip, i, g), dose.pgt1_bar_pooled(ip, i, g));
        end
    end
end
fclose(fid);
fprintf('  wrote %s\n', out_csv);

% The dose-response READ THROUGH THE PROTOCOL: tau is measured relative to
% lambda_h, and lambda_h is much tighter at p = 12 than at p = 2 (a
% richer fitted VAR), so the LEVEL of tau_bar is not comparable across p.
% What is comparable is how many cells exceed their own null q95 at each
% p.  Uses the null calibrations null_calibration_p<p>[...].mat if present.
n_esc = nan(numel(plist), 2);  n_cells = nan(numel(plist), 1);
for ip = 1:numel(plist)
    ptag = strrep(tag, sprintf('p%d', cfg.p), sprintf('p%d', plist(ip)));
    nf = fullfile(P.results, sprintf('null_calibration_%s.mat', ptag));
    if exist(nf, 'file') ~= 2, continue; end
    nl = load(nf);
    if ~isequal(nl.h_early(:)', h_early(:)'), continue; end
    tb = squeeze(dose.tau_bar_block(ip, :, :));
    free = true(1, K);  free(cfg.blocks.fixed_tau) = false;
    n_cells(ip) = K * sum(free);
    n_esc(ip, 1) = sum(sum(tb(:, free) > nl.block.q95(:, free)));
    if DO_POOLED && isfield(nl, 'pooled')
        tp = squeeze(dose.tau_bar_pooled(ip, :, :));
        n_esc(ip, 2) = sum(sum(tp(:, free) > nl.pooled.q95(:, free)));
    end
end
if any(~isnan(n_esc(:)))
    out_csv = fullfile(P.results, sprintf('dose_response_protocol_%s.csv', tag));
    fid = fopen(out_csv, 'w');
    fprintf(fid, 'p,n_free_cells,n_escapes_q95_block,n_escapes_q95_pooled\n');
    for ip = 1:numel(plist)
        fprintf(fid, '%d,%g,%g,%g\n', plist(ip), n_cells(ip), n_esc(ip, 1), n_esc(ip, 2));
        fprintf('  p = %2d: %g of %g free cells above their null q95 (block), %g (pooled)\n', ...
                plist(ip), n_esc(ip, 1), n_cells(ip), n_esc(ip, 2));
    end
    fclose(fid);
    fprintf('  wrote %s\n', out_csv);
else
    fprintf('  (no null calibration found for p in {%s}: the protocol count is skipped)\n', num2str(plist));
end
end

% ==============================================================================
if DO_D
fprintf('\n[D] Information-effect comparison (JK policy-only surprises) ...\n');
assert(exist(step_a_file, 'file') == 2, ...
       ['RUN_EMPIRICAL step D compares against %s, which step A writes.\n' ...
        'Set DO_A = true (once), or rerun with DO_D = false.'], step_a_file);
[dpth, dnm, dext] = fileparts(P.dataset);
jk_file = fullfile(dpth, sprintf('%s_mps_gc_jk%s', dnm, dext));
if exist(jk_file, 'file') ~= 2
    fprintf('  building the JK-variant dataset ...\n');
    assemble_dataset(struct('shock_variant', 'mps_gc_jk'));
end
dj = load(jk_file);
cfj = cfg;  cfj.fmar.isrw = dj.isrw;
rng(20260300, 'twister');
bvj = estimate_bvar_niw(dj.Y, cfj);
bfj = estimate_blp_fmar(dj.Y, cfj, bvj);
bbj = estimate_blp_blockadaptive(dj.Y, cfj, bvj, bfj.lambda);
if DO_POOLED, bpj = estimate_blp_blockpooled(dj.Y, cfj, bvj, bfj.lambda); else, bpj = []; end
relj = ea_relevance(dj.Y, bvj, cfj);
fprintf('  JK variant: BVAR lambda %.3f; impact of the surprise on i1y %.4f (t %.2f); lambda multimodal at %d h\n', ...
        bvj.lambda, relj.impact_b, relj.impact_t, bfj.diag.lambda_multimodal);
bvj = slim_bvar(bvj);
save(fullfile(P.results, sprintf('empirical_jk_%s.mat', tag)), 'bvj', 'bfj', 'bbj', 'bpj', 'cfj', 'relj', 'tag');
% side-by-side early-h tau maps, baseline vs JK
L = load(step_a_file);
assert(L.cfg.p == cfg.p && L.cfg.H == cfg.H, ...
       ['RUN_EMPIRICAL step D: %s was written at p = %d, H = %d but this run ' ...
        'uses p = %d, H = %d. Rerun step A before comparing.'], ...
       step_a_file, L.cfg.p, L.cfg.H, cfg.p, cfg.H);
tb_base = mean(L.blpb.tau_mean(:, :, h_early(1):h_early(2)), 3);
tb_jk   = mean(bbj.tau_mean(:, :, h_early(1):h_early(2)), 3);
tp_base = nan(K);  tp_jk = nan(K);
if DO_POOLED && isfield(L, 'blpp') && ~isempty(L.blpp)
    tp_base = mean(L.blpp.tau_mean(:, :, h_early(1):h_early(2)), 3);
    tp_jk   = mean(bpj.tau_mean(:, :, h_early(1):h_early(2)), 3);
end
out_csv = fullfile(P.results, sprintf('tau_jk_comparison_%s.csv', tag));
fid = fopen(out_csv, 'w');
assert(fid > 0, 'RUN_EMPIRICAL: cannot write %s', out_csv);
fprintf(fid, 'equation,block,tau_bar_baseline_block,tau_bar_jk_block,tau_bar_baseline_pooled,tau_bar_jk_pooled\n');
for i = 1:K
    for g = 1:K
        fprintf(fid, '%s,%s,%.5f,%.5f,%.5f,%.5f\n', ds.varnames{i}, ds.varnames{g}, ...
                tb_base(i, g), tb_jk(i, g), tp_base(i, g), tp_jk(i, g));
    end
end
fclose(fid);
fprintf('  wrote %s\n', out_csv);
end

% ==============================================================================
if DO_E
fprintf('\n[E] Coherence table and reading protocol ...\n');
assert(exist(step_a_file, 'file') == 2, ...
       ['RUN_EMPIRICAL step E needs %s, which step A writes.\n' ...
        'Set DO_A = true (once), or rerun with DO_E = false.'], step_a_file);
L = load(step_a_file);
have_pooled = isfield(L, 'blpp') && ~isempty(L.blpp);

% --- E1. mechanical coherence (CH7_DESIGN Sec. 5.3) ----------------------
% The block prior acts on COEFFICIENT BLOCKS, so that is where coherence
% is checked: for every (equation i, block g, horizon h) the distance of
% the block's posterior-mean coefficients from their UNSHRUNK value (OLS
% on the same detrended regressors, the tau -> infinity limit) under the
% adaptive estimator relative to the same distance under the global FMAR
% estimator,  ratio = ||b_adaptive,g - b_ols,g|| / ||b_fmar,g - b_ols,g||.
% A block that escapes (tau > 1) must have ratio < 1 (released toward the
% data); a block that is tightened (tau < 1) must have ratio > 1 (pulled
% further into the VAR centre); a held block must have ratio = 1.  An
% IRF-level version of this check is not informative: the IRF mixes all
% blocks' contemporaneous coefficients through the impact vector.
dtE = var_deterministic_trend(ds.Y, L.bvar.B, cfg.p);
xE  = dtE.x;  TxE = size(xE, 1);
ZallE = build_lp_regressors(xE, cfg.p);
bpE = build_block_prior(ds.Y, cfg);
out_csv = fullfile(P.results, sprintf('coherence_%s.csv', tag));
fid = fopen(out_csv, 'w');
assert(fid > 0, 'RUN_EMPIRICAL: cannot write %s', out_csv);
fprintf(fid, ['equation,block,block_held_at_1,h,lambda_h,tau_block,tau_pooled,' ...
              'd_fmar_ols,d_block_ols,d_pooled_ols,ratio_block,ratio_pooled\n']);
held = false(1, K);  held(cfg.blocks.fixed_tau) = true;
R_b = [];  T_b = [];  R_p = [];  T_p = [];
for h = 1:cfg.H
    Zh = ZallE(1:end - h, :);  Yh = xE(cfg.p + h:TxE, :);
    B_ols = Zh \ Yh;                                    % (m x K)
    for i = 1:K
        for g = 1:K
            idx = bpE.block_id == g;
            d_f = norm(L.blpf.beta_mean(idx, i, h) - B_ols(idx, i));
            d_b = norm(L.blpb.beta_mean(idx, i, h) - B_ols(idx, i));
            tb  = L.blpb.tau_mean(i, g, h);
            if have_pooled
                d_p = norm(L.blpp.beta_mean(idx, i, h) - B_ols(idx, i));
                tp  = L.blpp.tau_mean(i, g, h);
            else
                d_p = NaN;  tp = NaN;
            end
            fprintf(fid, '%s,%s,%d,%d,%.6g,%.5f,%.5f,%.6g,%.6g,%.6g,%.5f,%.5f\n', ...
                    ds.varnames{i}, ds.varnames{g}, held(g), h, L.blpf.lambda(1, h), ...
                    tb, tp, d_f, d_b, d_p, d_b / d_f, d_p / d_f);
            if ~held(g) && h >= 2 && isfinite(d_f) && d_f > 0
                R_b(end + 1) = d_b / d_f;  T_b(end + 1) = tb;         %#ok<SAGROW>
                if have_pooled, R_p(end + 1) = d_p / d_f;  T_p(end + 1) = tp; end  %#ok<SAGROW>
            end
        end
    end
end
fclose(fid);
fprintf('  wrote %s\n', out_csv);
coh_summary(R_b, T_b, 'BLP-block');
if have_pooled, coh_summary(R_p, T_p, 'BLP-pooled'); end

% --- E2. the reading protocol (CH7_DESIGN Sec. 5.1) ----------------------
null_file = fullfile(P.results, sprintf('null_calibration_%s.mat', tag));
if exist(null_file, 'file') ~= 2
    fprintf(['  no %s yet: run the null calibration\n' ...
             '      null = run_null_calibration(%s);\n' ...
             '  and rerun RUN_EMPIRICAL with only DO_E = true for the protocol table.\n'], ...
            null_file, tern(PSI_FLOOR, '', 'struct(''null'', struct(''psi_floor'', false))'));
else
    nl = load(null_file);
    if ~isfield(nl, 'block'), nl.block = nl; end          % pre-revision layout
    if nl.p ~= cfg.p
        fprintf(2, ['  %s was produced at p = %d but step A ran at p = %d;\n' ...
                    '  the thresholds do not apply.\n'], null_file, nl.p, cfg.p);
    else
        assert(isequal(nl.h_early(:)', h_early(:)'), ...
               'RUN_EMPIRICAL step E: null window [%d %d] differs from [%d %d].', ...
               nl.h_early(1), nl.h_early(2), h_early(1), h_early(2));
        sets = {'block', L.blpb};
        if have_pooled && isfield(nl, 'pooled'), sets(end + 1, :) = {'pooled', L.blpp}; end
        apply_protocol(sets, nl, ds.varnames, h_early, cfg.blocks.fixed_tau, ...
                       fullfile(P.results, sprintf('tau_protocol_%s.csv', tag)));
    end
end
end

% ==============================================================================
if DO_F
fprintf('\n[F] Four-variable level system without the surprise ...\n');
t_f = tic;
Y4 = ds.Y(:, 2:5);
vn4 = ds.varnames(2:5);
cf4 = cfg;  cf4.K = 4;  cf4.fmar.isrw = [1 1 1 1];  cf4.shock_var = 1;   % i1y innovation
cf4.blocks.fixed_tau = [];         % no instrument here: every block adapts
tag = tag_lev;                     % level-system files: lag order + scale convention only
rng(20260400, 'twister');
bv4 = estimate_bvar_niw(Y4, cf4);
bf4 = estimate_blp_fmar(Y4, cf4, bv4);
bb4 = estimate_blp_blockadaptive(Y4, cf4, bv4, bf4.lambda);
if DO_POOLED, bp4 = estimate_blp_blockpooled(Y4, cf4, bv4, bf4.lambda); else, bp4 = []; end
lp4 = estimate_lp_lagaug(Y4, cf4, bv4);
fprintf('  BVAR lambda %.3f, max|eig| %.4f; lambda_h multimodal at %d of %d horizons; psi floor bound at %d cells\n', ...
        bv4.lambda, bv4.max_eig, bf4.diag.lambda_multimodal, cfg.H, sum(bf4.diag.psi_floor_bound));
fprintf('  lambda_h, h = 1..%d:\n   ', cfg.H);
fprintf(' %.3f', bf4.lambda(1, :));  fprintf('\n');
bv4 = slim_bvar(bv4);
ds4_info = struct('varnames', {vn4}, 'ym', ds.ym, 'meta', ds.meta);
save(fullfile(P.results, sprintf('empirical_levels_%s.mat', tag)), ...
     'bv4', 'bf4', 'bb4', 'bp4', 'lp4', 'cf4', 'h_early', 'ds4_info', 'tag');
% a dataset file for a level-system null calibration
ds4 = struct('Y', Y4, 'varnames', {vn4}, 'ym', ds.ym, 'isrw', [1 1 1 1], ...
             'shock_variant', 'none (level system, i1y innovation)', ...
             'meta', ds.meta, 'synthetic', false);
save(fullfile(P.data, 'ea_dataset_levels.mat'), '-struct', 'ds4');
out_csv = fullfile(P.results, sprintf('tau_heatmap_levels_%s.csv', tag));
fid = fopen(out_csv, 'w');
assert(fid > 0, 'RUN_EMPIRICAL: cannot write %s', out_csv);
fprintf(fid, 'scale_estimator,equation,block,h,lambda_h,tau_mean,tau_median,p_tau_gt1\n');
ests4 = {'block', bb4};
if DO_POOLED, ests4(end + 1, :) = {'pooled', bp4}; end
for e = 1:size(ests4, 1)
    est = ests4{e, 2};
    for i = 1:4
        for g = 1:4
            for h = 1:cfg.H
                fprintf(fid, '%s,%s,%s,%d,%.6g,%.5f,%.5f,%.4f\n', ests4{e, 1}, vn4{i}, vn4{g}, h, ...
                        bf4.lambda(1, h), est.tau_mean(i, g, h), est.tau_med(i, g, h), est.p_tau_gt1(i, g, h));
            end
        end
    end
end
fclose(fid);
try
    for e = 1:size(ests4, 1)
        est = ests4{e, 2};
        figure('visible', 'off', 'Position', [50 50 1100 700]);
        for i = 1:4
            subplot(2, 2, i);
            imagesc(log(squeeze(est.tau_mean(i, :, :))));  colorbar
            title(sprintf('LEVELS, eq %s: log tau (%s) by block x h', vn4{i}, ests4{e, 1}));
            set(gca, 'YTick', 1:4, 'YTickLabel', vn4);
            xlabel('horizon h');
        end
        print(fullfile(P.results, sprintf('fig_tau_heatmap_levels_%s_%s.png', ests4{e, 1}, tag)), '-dpng', '-r110');
    end
catch err
    fprintf('  (plotting skipped: %s)\n', err.message);
end
% the reading protocol for the level system, if its null has been run
nf4 = fullfile(P.results, sprintf('null_calibration_%s_levels.mat', tag_lev));
if exist(nf4, 'file') == 2
    nl4 = load(nf4);
    sets4 = {'block', bb4};
    if DO_POOLED && isfield(nl4, 'pooled'), sets4(end + 1, :) = {'pooled', bp4}; end
    apply_protocol(sets4, nl4, vn4, h_early, [], ...
                   fullfile(P.results, sprintf('tau_protocol_levels_%s.csv', tag_lev)));
else
    fprintf(['  (no %s yet: run\n' ...
             '     run_null_calibration(struct(''null'', struct(''dataset'', ''%s'', ' ...
             '''fixed_tau'', [], ''stem'', ''%s_levels'')))\n   for the level-system protocol)\n'], ...
            nf4, fullfile(P.data, 'ea_dataset_levels.mat'), tag_lev);
end
fprintf('  step F done in %.1f min; wrote %s\n', toc(t_f) / 60, out_csv);
end

fprintf('\nRUN_EMPIRICAL finished.  Outputs are in %s (tag %s)\n', P.results, tag);
fprintf('Reminder: the null calibration (empirical/montecarlo/run_null_calibration.m)\n');
fprintf('is a SEPARATE, longer job.  Do a timing run first:\n');
fprintf('  null = run_null_calibration(struct(''null'', struct(''n_rep'', 10)));\n');
fprintf('then scale n_rep from the printed seconds-per-rep, and rerun this script\n');
fprintf('with only DO_E = true to apply the reading protocol.\n');

% -------------------------------------------------------------------------
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function apply_protocol(sets, nl, vn, h_early, held_blocks, out_csv)
% The pre-registered reading protocol (CH7_DESIGN Sec. 5.1) applied to one
% or two scale estimators against a stored null calibration.  Blocks held
% at tau = 1 are excluded from every comparison (they carry no signal).
K = numel(vn);
free = true(1, K);  free(held_blocks) = false;
fid = fopen(out_csv, 'w');
assert(fid > 0, 'RUN_EMPIRICAL: cannot write %s', out_csv);
fprintf(fid, ['scale_estimator,equation,block,block_held_at_1,tau_bar,null_mean,null_q90,null_q95,' ...
              'null_q99,escape_q95,pgt1_bar,null_pgt1_q95,is_argmax,null_argmax_freq,' ...
              'eq_maxstat,eq_maxstat_q95,eq_flag\n']);
for e = 1:size(sets, 1)
    nm = sets{e, 1};  est = sets{e, 2};  S = nl.(nm);
    tb = mean(est.tau_mean(:, :, h_early(1):h_early(2)), 3);        % K x G
    pb = mean(est.p_tau_gt1(:, :, h_early(1):h_early(2)), 3);
    fprintf('  [%s] reading protocol (null R = %d, p = %d, window h = %d..%d):\n', ...
            nm, nl.n_rep, nl.p, h_early(1), h_early(2));
    for i = 1:K
        tbi = tb(i, :);  tbi(~free) = -inf;
        [mx, am] = max(tbi);
        flag_eq = mx > S.maxstat_q95(i);
        fprintf('    eq %-6s max tau_bar %.3f (null q95 %.3f) %s | argmax block %s (null freq %.2f)', ...
                vn{i}, mx, S.maxstat_q95(i), tern(flag_eq, 'ESCAPE', 'quiet '), ...
                vn{am}, S.argmax_freq(i, am));
        esc = find(free & tb(i, :) > S.q95(i, :));
        if isempty(esc), fprintf(' | no cell above q95\n');
        else, fprintf(' | cells above q95: %s\n', strjoin(vn(esc), ', ')); end
        for g = 1:K
            pq = NaN;
            if isfield(S, 'pgt1_q95'), pq = S.pgt1_q95(i, g); end
            fprintf(fid, '%s,%s,%s,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%d,%.4f,%.4f,%d,%.4f,%.5f,%.5f,%d\n', ...
                    nm, vn{i}, vn{g}, ~free(g), tb(i, g), S.tau_bar_mean(i, g), ...
                    S.q90(i, g), S.q95(i, g), S.q99(i, g), free(g) && tb(i, g) > S.q95(i, g), ...
                    pb(i, g), pq, g == am, S.argmax_freq(i, g), mx, S.maxstat_q95(i), flag_eq);
        end
    end
end
fclose(fid);
fprintf('  wrote %s\n', out_csv);
end

function b = slim_bvar(b)
% Drop the per-draw companion matrices before saving (see step A).
if isfield(b, 'F_draws'), b = rmfield(b, 'F_draws'); end
end

function coh_summary(R, T, nm)
% Print the coherence check: the distance ratio by tau class, and the rank
% correlation between log tau and log ratio over free blocks, h >= 2.
if isempty(R), return; end
cls = {T > 2, T > 1 & T <= 2, T <= 1};
lbl = {'tau > 2', '1 < tau <= 2', 'tau <= 1'};
fprintf('  coherence, %s (free blocks, h >= 2, %d cells): ratio = |b_adapt - b_ols| / |b_fmar - b_ols|\n', ...
        nm, numel(R));
for c = 1:3
    if any(cls{c})
        fprintf('    %-13s n = %4d  median ratio %.3f  share < 1: %.2f\n', lbl{c}, ...
                sum(cls{c}), median(R(cls{c})), mean(R(cls{c}) < 1));
    end
end
% Spearman correlation via ranks (no toolbox)
[~, o1] = sort(log(T));  r1 = zeros(size(T));  r1(o1) = 1:numel(T);
[~, o2] = sort(log(R));  r2 = zeros(size(R));  r2(o2) = 1:numel(R);
rs = corrcoef(r1, r2);
fprintf('    Spearman corr(log tau, log ratio) = %.3f  (coherent = strongly negative)\n', rs(1, 2));
end
