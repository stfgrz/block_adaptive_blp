% RUN_EMPIRICAL  Chapter 7 driver: block-adaptive BLP on euro-area data.
%
% Design document: docs/CH7_DESIGN.md.  Integration notes and assumed repo
% signatures: README_EMPIRICAL.md.  Run AFTER:
%   1. build_shock_series();        (needs data/raw/ea_empd_events.csv)
%   2. fetch_outcome_data();        (or manual downloads, see that file)
%   3. assemble_dataset();          (baseline variant mps_gc_1y)
%   4. assemble_dataset(struct('shock_variant', 'mps_gc_jk'));   % step D
%
% STEPS (toggle below)
%   A  baseline estimation at p = 12: BVAR, BLP-FMAR, BLP-block, LP-LA
%   B  tau heatmap export (the Leg-1 finding)
%   C  dose-response in p in {2, 6, 12} (Leg 3; see CH7_DESIGN Sec. 5.2)
%   D  JK information-effect comparison (policy-only vs baseline surprises)
%   E  null calibration pointer (run montecarlo/run_null_calibration.m
%      separately -- it is the long job)

DO_A = true;  DO_B = true;  DO_C = true;  DO_D = true;

addpath(genpath(pwd));
if ~exist('results', 'dir'), mkdir('results'); end

% --- configuration -----------------------------------------------------------
ds  = load(fullfile('data', 'ea_dataset.mat'));
cfg = default_config();
cfg.mode      = 'fmar';
cfg.shock_var = 1;                 % surprise ordered first (internal instrument)
cfg.p         = 12;                % monthly baseline lag order
cfg.H         = 48;                % 4-year horizon
cfg.ci_level  = 0.90;
cfg.fmar.isrw = ds.isrw;           % [0 1 1 1 1] -- REQUIRES the vectorised
                                   % isrw patch, see README_EMPIRICAL.md
K = size(ds.Y, 2);

% ==============================================================================
if DO_A
fprintf('\n[A] Baseline estimation, p = %d, H = %d ...\n', cfg.p, cfg.H);
rng(20260107, 'twister');
bvar = estimate_bvar_niw(ds.Y, cfg);
blpf = estimate_blp_fmar(ds.Y, cfg, bvar);
blpb = estimate_blp_blockadaptive(ds.Y, cfg, bvar, blpf.lambda);
lpla = estimate_lp_lagaug(ds.Y, cfg, bvar);

% 25 bp normalisation: common factor from the BVAR impact response of the
% 1Y rate (variable 2) to the surprise innovation.  Applied at REPORTING
% time only; all saved IRFs are in unit-surprise-innovation scale.
k25 = 0.25 / bvar.theta(2, 1);
fprintf('  BVAR lambda = %.3f, max |eig| = %.4f; 25bp scale k = %.3f\n', ...
        bvar.lambda, bvar.max_eig, k25);

save(fullfile('results', 'empirical_p12.mat'), ...
     'bvar', 'blpf', 'blpb', 'lpla', 'cfg', 'k25');

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
    print(fullfile('results', 'fig_irf_p12.png'), '-dpng');
catch err
    fprintf('  (plotting skipped: %s)\n', err.message);
end
end

% ==============================================================================
if DO_B
fprintf('\n[B] Tau heatmap export ...\n');
L = load(fullfile('results', 'empirical_p12.mat'));
tau = L.blpb.tau_mean;                                   % K x G x H
fid = fopen(fullfile('results', 'tau_heatmap_p12.csv'), 'w');
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
    print(fullfile('results', 'fig_tau_heatmap_p12.png'), '-dpng');
catch err
    fprintf('  (plotting skipped: %s)\n', err.message);
end
fprintf('  wrote results/tau_heatmap_p12.csv\n');
end

% ==============================================================================
if DO_C
fprintf('\n[C] Dose-response in p (Leg 3) ...\n');
plist = [2 6 12];
h1 = 2;  h2 = 12;                                  % early-h window (protocol)
dose = struct('p', plist, 'tau_bar', nan(numel(plist), K, K));
for ip = 1:numel(plist)
    cfp = cfg;  cfp.p = plist(ip);
    rng(20260200 + plist(ip), 'twister');
    bv = estimate_bvar_niw(ds.Y, cfp);
    bf = estimate_blp_fmar(ds.Y, cfp, bv);
    bb = estimate_blp_blockadaptive(ds.Y, cfp, bv, bf.lambda);
    dose.tau_bar(ip, :, :) = mean(bb.tau_mean(:, :, h1:h2), 3);
    fprintf('  p = %2d done (BVAR lambda %.3f)\n', plist(ip), bv.lambda);
end
save(fullfile('results', 'dose_response.mat'), '-struct', 'dose');
fid = fopen(fullfile('results', 'dose_response.csv'), 'w');
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
fprintf('  wrote results/dose_response.csv\n');
end

% ==============================================================================
if DO_D
fprintf('\n[D] Information-effect comparison (JK policy-only surprises) ...\n');
jk_file = fullfile('data', 'ea_dataset_mps_gc_jk.mat');
if ~exist(jk_file, 'file')
    fprintf('  building the JK-variant dataset ...\n');
    assemble_dataset(struct('shock_variant', 'mps_gc_jk'));
end
dj = load(jk_file);
cfj = cfg;  cfj.fmar.isrw = dj.isrw;
rng(20260300, 'twister');
bvj = estimate_bvar_niw(dj.Y, cfj);
bfj = estimate_blp_fmar(dj.Y, cfj, bvj);
bbj = estimate_blp_blockadaptive(dj.Y, cfj, bvj, bfj.lambda);
save(fullfile('results', 'empirical_jk.mat'), 'bvj', 'bfj', 'bbj', 'cfj');
% side-by-side early-h tau maps, baseline vs JK
L = load(fullfile('results', 'empirical_p12.mat'));
tb_base = mean(L.blpb.tau_mean(:, :, 2:12), 3);
tb_jk   = mean(bbj.tau_mean(:, :, 2:12), 3);
fid = fopen(fullfile('results', 'tau_jk_comparison.csv'), 'w');
fprintf(fid, 'equation,block,tau_bar_baseline,tau_bar_jk\n');
for i = 1:K
    for g = 1:K
        fprintf(fid, '%s,%s,%.5f,%.5f\n', ds.varnames{i}, ds.varnames{g}, ...
                tb_base(i, g), tb_jk(i, g));
    end
end
fclose(fid);
fprintf('  wrote results/tau_jk_comparison.csv\n');
end

fprintf('\nRUN_EMPIRICAL finished.  Reminder: the null-calibration job\n');
fprintf('(montecarlo/run_null_calibration.m) is separate; do a quick run first:\n');
fprintf('  null = run_null_calibration(struct(''null'', struct(''n_rep'', 10)));\n');
