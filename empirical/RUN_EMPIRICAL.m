% RUN_EMPIRICAL  Chapter 7 driver: block-adaptive BLP on euro-area data.
%
% Design document: docs/CH7_DESIGN.md.  Integration notes and assumed repo
% signatures: README_EMPIRICAL.md.
%
% PREREQUISITES (run once, in this order; each is checked below and the
% script stops with an actionable message if one is missing)
%   1. build_shock_series();        (needs data/raw/ea_empd_events.csv,
%                                    shipped with the package)
%   2. fetch_outcome_data();        (must report "all four outcome files
%                                    ready and validated"; if the downloads
%                                    fail it prints the manual route)
%   3. assemble_dataset();          (baseline variant mps_gc_1y)
%   4. SMOKE_TEST_EMPIRICAL         (5 minutes; verifies the estimator
%                                    signatures before you commit a night
%                                    to them)
% Step D builds its own JK-variant dataset if it is not already there, so
% assemble_dataset(struct('shock_variant','mps_gc_jk')) is optional.
%
% STEPS (toggle below)
%   A  baseline estimation at p = 12: BVAR, BLP-FMAR, BLP-block, LP-LA
%   B  tau heatmap export (the Leg-1 finding)          [needs A's output]
%   C  dose-response in p in {2, 6, 12} (Leg 3; see CH7_DESIGN Sec. 5.2)
%   D  JK information-effect comparison                [needs A's output]
%   E  null calibration pointer (empirical/montecarlo/run_null_calibration.m
%      is a SEPARATE, longer job -- do not queue it in the same session)
%
% RESTARTING: every step saves its own .mat into results/, so if the run
% dies part-way you can set the completed steps' DO_ toggles to false and
% rerun.  B and D read results/empirical_p12.mat, which step A writes; if
% you switch A off before it has ever run, they stop with a clear message
% instead of a cryptic load error.
%
% RUNTIME: steps A+C+D are five full estimations (one BVAR + K*H lambda
% searches + K*H Gibbs chains each).  Expect a few hours at p = 12, H = 48.
% Plotting is wrapped in try/catch, so a headless session loses the pngs
% but keeps every .mat and .csv -- those are the real outputs.

DO_A = true;  DO_B = true;  DO_C = true;  DO_D = true;

% --- path bootstrap (works from any directory, and when pasted) ----------
ea_this = mfilename('fullpath');
if isempty(ea_this), ea_this = fullfile(pwd, 'RUN_EMPIRICAL'); end
addpath(genpath(fileparts(fileparts(ea_this))));
P = ea_paths();
if exist(P.results, 'dir') ~= 7, mkdir(P.results); end

% (defined after cfg below, since the filename carries cfg.p)

% --- configuration -----------------------------------------------------------
assert(exist(P.dataset, 'file') == 2, ...
       ['RUN_EMPIRICAL: %s not found.\n' ...
        'Run the prerequisites first:\n' ...
        '    build_shock_series(); fetch_outcome_data(); assemble_dataset();\n' ...
        'then SMOKE_TEST_EMPIRICAL before launching this.'], P.dataset);
ds  = load(P.dataset);
cfg = default_config();
cfg.mode      = 'fmar';
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

% Output filenames carry the baseline lag order, so a run at a different p
% cannot silently overwrite the p = 12 results the write-up cites.
tag = sprintf('p%d', cfg.p);
step_a_file = fullfile(P.results, sprintf('empirical_%s.mat', tag));

fprintf('RUN_EMPIRICAL: %s\n  T = %d, K = %d, p = %d, H = %d, early-h window [%d %d], variant = %s\n', ...
        P.dataset, size(ds.Y, 1), K, cfg.p, cfg.H, h_early(1), h_early(2), ds.shock_variant);

% ==============================================================================
if DO_A
fprintf('\n[A] Baseline estimation, p = %d, H = %d ...\n', cfg.p, cfg.H);
t_a = tic;
rng(20260107, 'twister');
bvar = estimate_bvar_niw(ds.Y, cfg);
blpf = estimate_blp_fmar(ds.Y, cfg, bvar);
blpb = estimate_blp_blockadaptive(ds.Y, cfg, bvar, blpf.lambda);
lpla = estimate_lp_lagaug(ds.Y, cfg, bvar);

% 25 bp normalisation: common factor from the BVAR impact response of the
% 1Y rate (variable 2) to the surprise innovation.  Applied at REPORTING
% time only; all saved IRFs are in unit-surprise-innovation scale.
imp_i1y = bvar.theta(2, 1);
k25 = 0.25 / imp_i1y;
fprintf('  BVAR lambda = %.3f, max |eig| = %.4f; 25bp scale k = %.3f\n', ...
        bvar.lambda, bvar.max_eig, k25);
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
if bvar.max_eig >= 1
    fprintf('  note: fitted BVAR has max |eig| = %.4f >= 1 (expected with RW\n', bvar.max_eig);
    fprintf('        prior means on levels; the null bootstrap conditions on\n');
    fprintf('        the actual first p observations, see simulate_fitted_bvar_dgp).\n');
end

save(step_a_file, 'bvar', 'blpf', 'blpb', 'lpla', 'cfg', 'k25');
fprintf('  step A done in %.1f min -> %s\n', toc(t_a) / 60, step_a_file);

% quick-look IRF figure (one panel per variable, 25bp-normalised)
try
    figure('visible', 'off');
    hs = 0:cfg.H;
    for i = 1:K
        subplot(ceil(K / 2), 2, i);  hold on
        plot(hs, k25 * bvar.theta(i, :), 'r-');
        plot(hs, k25 * blpf.theta_mean(i, :), 'b-');
        plot(hs, k25 * blpb.theta_mean(i, :), 'k-');
        plot(hs, k25 * lpla.theta(i, :), 'g--');
        plot(hs, k25 * blpb.lo(i, :), 'k:');  plot(hs, k25 * blpb.hi(i, :), 'k:');
        title(ds.varnames{i});  if i == 1, legend('BVAR', 'BLP-FMAR', 'BLP-block', 'LP-LA'); end
    end
    print(fullfile(P.results, sprintf('fig_irf_%s.png', tag)), '-dpng');
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
tau = L.blpb.tau_mean;                                   % K x G x H
out_csv = fullfile(P.results, sprintf('tau_heatmap_%s.csv', tag));
fid = fopen(out_csv, 'w');
assert(fid > 0, 'RUN_EMPIRICAL: cannot write %s', out_csv);
fprintf(fid, 'equation,block,h,tau_mean\n');
for i = 1:K
    for g = 1:K
        for h = 1:size(tau, 3)
            fprintf(fid, '%s,%s,%d,%.5f\n', ds.varnames{i}, ds.varnames{g}, ...
                    h, tau(i, g, h));
        end
    end
end
fclose(fid);
try
    figure('visible', 'off');
    for i = 1:K
        subplot(ceil(K / 2), 2, i);
        imagesc(squeeze(tau(i, :, :)));  colorbar
        title(sprintf('eq %s: tau by block x h', ds.varnames{i}));
        set(gca, 'YTick', 1:K, 'YTickLabel', ds.varnames);
    end
    print(fullfile(P.results, sprintf('fig_tau_heatmap_%s.png', tag)), '-dpng');
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
dose = struct('p', plist, 'tau_bar', nan(numel(plist), K, K), ...
              'h_early', h_early, 'H', cfg.H);
for ip = 1:numel(plist)
    cfp = cfg;  cfp.p = plist(ip);
    rng(20260200 + plist(ip), 'twister');
    t_c = tic;
    bv = estimate_bvar_niw(ds.Y, cfp);
    bf = estimate_blp_fmar(ds.Y, cfp, bv);
    bb = estimate_blp_blockadaptive(ds.Y, cfp, bv, bf.lambda);
    dose.tau_bar(ip, :, :) = mean(bb.tau_mean(:, :, h1:h2), 3);
    fprintf('  p = %2d done in %.1f min (BVAR lambda %.3f)\n', ...
            plist(ip), toc(t_c) / 60, bv.lambda);
end
save(fullfile(P.results, 'dose_response.mat'), '-struct', 'dose');
out_csv = fullfile(P.results, 'dose_response.csv');
fid = fopen(out_csv, 'w');
assert(fid > 0, 'RUN_EMPIRICAL: cannot write %s', out_csv);
fprintf(fid, 'p,equation,block,tau_bar_early\n');
for ip = 1:numel(plist)
    for i = 1:K
        for g = 1:K
            fprintf(fid, '%d,%s,%s,%.5f\n', plist(ip), ds.varnames{i}, ...
                    ds.varnames{g}, dose.tau_bar(ip, i, g));
        end
    end
end
fclose(fid);
fprintf('  wrote %s\n', out_csv);
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
save(fullfile(P.results, 'empirical_jk.mat'), 'bvj', 'bfj', 'bbj', 'cfj');
% side-by-side early-h tau maps, baseline vs JK
L = load(step_a_file);
assert(L.cfg.p == cfg.p && L.cfg.H == cfg.H, ...
       ['RUN_EMPIRICAL step D: %s was written at p = %d, H = %d but this run ' ...
        'uses p = %d, H = %d. Rerun step A before comparing.'], ...
       step_a_file, L.cfg.p, L.cfg.H, cfg.p, cfg.H);
tb_base = mean(L.blpb.tau_mean(:, :, h_early(1):h_early(2)), 3);
tb_jk   = mean(bbj.tau_mean(:, :, h_early(1):h_early(2)), 3);
out_csv = fullfile(P.results, 'tau_jk_comparison.csv');
fid = fopen(out_csv, 'w');
assert(fid > 0, 'RUN_EMPIRICAL: cannot write %s', out_csv);
fprintf(fid, 'equation,block,tau_bar_baseline,tau_bar_jk\n');
for i = 1:K
    for g = 1:K
        fprintf(fid, '%s,%s,%.5f,%.5f\n', ds.varnames{i}, ds.varnames{g}, ...
                tb_base(i, g), tb_jk(i, g));
    end
end
fclose(fid);
fprintf('  wrote %s\n', out_csv);
end

fprintf('\nRUN_EMPIRICAL finished.  Outputs are in %s\n', P.results);
fprintf('Reminder: the null calibration (empirical/montecarlo/run_null_calibration.m)\n');
fprintf('is a SEPARATE, longer job.  Do a timing run first:\n');
fprintf('  null = run_null_calibration(struct(''null'', struct(''n_rep'', 10)));\n');
fprintf('then scale n_rep from the printed seconds-per-rep.\n');
