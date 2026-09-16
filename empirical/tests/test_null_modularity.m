function test_null_modularity()
% PURPOSE
% -------
% Checks the null-modularity layer of the Ch. 7 protocol on fakes built
% in memory (no data, no sampler, well under 20 s):
%
%   1. ea_design_key is deterministic, changes with every design-relevant
%      choice (p, T, a variable name, h_early, psi_floor, fixed_tau) and
%      does NOT change with the identification label;
%   2. ea_apply_protocol on a fake null with the exact field layout of
%      run_null_calibration's null_stats: planted escapes are counted,
%      percentiles lie in [0,1], Holm rejects no more than q95 does, the
%      held block is excluded, a mismatched design key is refused (and
%      only warned about under opts.force), the flat legacy layout works;
%   3. ea_cross_p_table on two temporary result/null pairs, one with a
%      matching key and one with a mismatched key.
%
% Every number here comes from a lognormal fake and is discarded.

fprintf('test_null_modularity:\n');
here = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

% --- 1. the design key ----------------------------------------------------
vn = {'a', 'b', 'c'};
cfg = default_config();
cfg.mode = 'fmar';  cfg.p = 2;  cfg.H = 4;
cfg.fmar.h1_mode = 'lp';  cfg.fmar.psi_floor = true;  cfg.fmar.isrw = [0 1 1];
cfg.blocks.fixed_tau = 1;
ex = struct('T', 200, 'h_early', [2 4]);
[k0, info0] = ea_design_key(vn, cfg, ex);
k0b = ea_design_key(vn, cfg, ex);
assert(strcmp(k0, k0b), 'design key is not deterministic');
assert(numel(info0.hash) == 8 && all(ismember(info0.hash, '0123456789abcdef')), ...
       'hash is not 8 hex characters: %s', info0.hash);
assert(strcmp(k0(end - 7:end), info0.hash), 'key does not end with its hash');
% a dataset struct gives the same key as the cellstr + T + isrw route
ds = struct('varnames', {vn}, 'Y', zeros(200, 3), 'isrw', [0 1 1]);
assert(strcmp(ea_design_key(ds, cfg, struct('h_early', [2 4])), k0), ...
       'dataset-struct and cellstr routes disagree');

variants = cell(0, 2);
c = cfg;  c.p = 3;                        variants(end + 1, :) = {'p', ea_design_key(vn, c, ex)};
variants(end + 1, :) = {'T', ea_design_key(vn, cfg, struct('T', 201, 'h_early', [2 4]))};
variants(end + 1, :) = {'varname', ea_design_key({'a', 'b', 'd'}, cfg, ex)};
variants(end + 1, :) = {'h_early', ea_design_key(vn, cfg, struct('T', 200, 'h_early', [2 3]))};
c = cfg;  c.fmar.psi_floor = false;       variants(end + 1, :) = {'psi_floor', ea_design_key(vn, c, ex)};
c = cfg;  c.blocks.fixed_tau = [];        variants(end + 1, :) = {'fixed_tau', ea_design_key(vn, c, ex)};
% H is deliberately NOT part of the key (the statistic lives on h <= h_early(2);
% the null is run at a shorter H than the IRF run): the key must be blind to it
c = cfg;  c.H = 6;
assert(strcmp(ea_design_key(vn, c, ex), k0), 'design key must not change with H');
for k = 1:size(variants, 1)
    assert(~strcmp(variants{k, 2}, k0), 'design key did not change with %s', variants{k, 1});
end
for k = 1:size(variants, 1)
    for j = k + 1:size(variants, 1)
        assert(~strcmp(variants{k, 2}, variants{j, 2}), ...
               'variants %s and %s collide', variants{k, 1}, variants{j, 1});
    end
end
% identification is recorded but not keyed; chains only on request
exi = ex;  exi.ident = 'proxy:mps_gc_1y';
[ki, infoi] = ea_design_key(vn, cfg, exi);
assert(strcmp(ki, k0) && strcmp(infoi.ident, 'proxy:mps_gc_1y'), ...
       'identification label must be recorded but must not enter the key');
exc = ex;  exc.include_chains = true;
assert(~strcmp(ea_design_key(vn, cfg, exc), k0), 'include_chains had no effect');
% fixed_tau order does not matter; a live mask does
c = cfg;  c.blocks.fixed_tau = [3 1];  ka = ea_design_key(vn, c, ex);
c.blocks.fixed_tau = [1 3];            assert(strcmp(ka, ea_design_key(vn, c, ex)), 'fixed_tau order leaked into the key');
c = cfg;  c.blocks.fixed_tau_mask = false(3);  assert(strcmp(ea_design_key(vn, c, ex), k0), 'an all-false mask changed the key');
c.blocks.fixed_tau_mask(2, 3) = true;         assert(~strcmp(ea_design_key(vn, c, ex), k0), 'a live mask did not change the key');
% the hash is real FNV-1a: the local copy matches the published vectors,
% and ea_design_key's hash is that of its own canonical string
assert(strcmp(fnv_of(''), '811c9dc5') && strcmp(fnv_of('a'), 'e40c292c'), 'FNV-1a reference vectors fail');
assert(strcmp(fnv_of(info0.canonical), info0.hash), 'ea_design_key hash is not FNV-1a of the canonical string');
fprintf('  ea_design_key: deterministic, sensitive to %d design fields, blind to ident: OK\n', ...
        size(variants, 1));

% --- 2. the protocol on a fake null -----------------------------------------
K = 3;  G = 3;  H = 4;  R = 200;  h_early = [2 4];
rng(4242, 'twister');
tau_all = exp(0.3 * randn(R, K, G, H));
pg_all  = rand(R, K, G, H);
cover_all = rand(R, K, H - 1) < 0.9;
nl = struct();
nl.block = fake_null_stats(tau_all, pg_all, cover_all, h_early(1), h_early(2));
nl.pooled = fake_null_stats(exp(0.2 * randn(R, K, G, H)), rand(R, K, G, H), cover_all, h_early(1), h_early(2));
nl.h_early = h_early;  nl.n_rep = R;  nl.p = cfg.p;  nl.H = H;  nl.varnames = vn;
nl.h1_mode = 'lp';  nl.psi_floor = true;  nl.fixed_tau = 1;  nl.T = 200;
nl.design_key = k0;
for f = {'q90', 'q95', 'q99', 'tau_bar_mean', 'argmax_freq', 'maxstat_q95', ...
         'coverage_under_null', 'tau_bar_draws'}
    nl.(f{1}) = nl.block.(f{1});
end

% fake estimator output: at the null median everywhere, two planted escapes
% in free cells, one huge value in the HELD block (must be ignored)
med = exp(0);                                  % lognormal(0, .) median
S = nl.block;
big = 5 * max(S.q99(:));
est = struct('tau_mean', med * ones(K, G, H), 'p_tau_gt1', 0.5 * ones(K, G, H));
est.tau_mean(1, 2, :) = big;                   % eq a <- block b
est.tau_mean(2, 3, :) = big;                   % eq b <- block c
est.tau_mean(3, 1, :) = 10 * big;              % eq c <- HELD block a
held = 1;
estp = struct('tau_mean', med * ones(K, G, H), 'p_tau_gt1', 0.5 * ones(K, G, H));

csv = fullfile(tempdir, sprintf('tau_protocol_test_%d.csv', round(rand * 1e6)));
prot = ea_apply_protocol({'block', est; 'pooled', estp}, nl, vn, h_early, held, csv, ...
                         struct('design_key', k0, 'quiet', true));
assert(strcmp(prot.key_status, 'match'), 'key_status = %s, expected match', prot.key_status);
assert(isequal(prot.names, {'block', 'pooled'}), 'both estimators should be processed');
B = prot.block;
assert(B.n_free_cells == K * (K - 1), 'n_free_cells = %d', B.n_free_cells);
assert(B.n_escape_q95 == 2, 'n_escape_q95 = %d, expected 2', B.n_escape_q95);
assert(B.escape_q95(1, 2) && B.escape_q95(2, 3) && ~B.escape_q95(3, 1), 'escape pattern wrong');
assert(~any(B.holm_reject(:, held)), 'held block was Holm-rejected');
assert(B.n_escape_holm <= B.n_escape_q95, 'Holm (%d) rejected more than q95 (%d)', ...
       B.n_escape_holm, B.n_escape_q95);
assert(B.n_escape_holm == 2, 'planted p = 0 cells should survive Holm, got %d', B.n_escape_holm);
assert(B.n_eq_flag == 2 && B.eq_flag(1) && B.eq_flag(2) && ~B.eq_flag(3), ...
       'equation flags wrong (held block must not flag eq c)');
pf = B.null_pct(:, prot.free);
assert(all(pf(:) >= 0 & pf(:) <= 1) && all(B.p_cell(:) >= 0 & B.p_cell(:) <= 1), ...
       'percentiles outside [0,1]');
assert(B.null_pct(1, 2) == 1 && B.p_cell(1, 2) == 0, 'planted escape should sit at percentile 1');
assert(abs(B.null_pct(2, 2) - 0.5) < 0.12, 'null median cell has percentile %.2f', B.null_pct(2, 2));
assert(all(abs(B.ratio_q95(:) - B.tau_bar(:) ./ S.q95(:)) < 1e-12), 'ratio_q95 definition');
% Monte Carlo p-values: never exactly 0, planted cells at the floor 1/(R+1);
% the family-wise max-statistic test rejects exactly the two planted cells
assert(all(B.p_mc(:) > 0) && abs(B.p_mc(1, 2) - 1 / (R + 1)) < 1e-12, 'p_mc convention');
assert(B.n_escape_fwer == 2 && B.fwer_reject(1, 2) && B.fwer_reject(2, 3) && ~B.fwer_reject(3, 1), 'FWER max-stat test');
assert(abs(B.p_fwer_min - 1 / (R + 1)) < 1e-12 && B.null_R == R && strcmp(B.null_method, 'resample'), 'p_fwer_min / null_R / null_method');
assert(all(B.p_fwer(:) >= B.p_mc(:) - 1e-12), 'the family-wise p-value must not be below the per-cell one');
Pp = prot.pooled;
assert(Pp.n_escape_q95 == 0 && Pp.n_escape_holm == 0 && Pp.n_eq_flag == 0, ...
       'the null-median pooled fake should escape nowhere');
% csv: header carries the new columns, the held block is marked
txt = fileread(csv);
lines = strsplit(strtrim(txt), sprintf('\n'));
assert(~isempty(strfind(lines{1}, ',ratio_q95,null_pct,p_cell,holm_reject,p_mc,p_fwer,fwer_reject')), 'csv header lacks the new columns');
assert(numel(lines) == 1 + 2 * K * G, 'csv has %d lines, expected %d', numel(lines), 1 + 2 * K * G);
cols = strsplit(lines{2}, ',');
assert(numel(cols) == 24, 'csv row has %d columns, expected 24', numel(cols));
held_rows = 0;
for k = 2:numel(lines)
    cc = strsplit(lines{k}, ',');
    if strcmp(cc{3}, vn{held}), assert(strcmp(cc{4}, '1')); held_rows = held_rows + 1; end
end
assert(held_rows == 2 * K, 'held-block rows = %d', held_rows);
delete(csv);
fprintf('  ea_apply_protocol: counts, percentiles, Holm, held block, csv: OK\n');

% design-key refusal, and the forced path
threw = false;
try
    ea_apply_protocol({'block', est}, nl, vn, h_early, held, '', ...
                      struct('design_key', variants{1, 2}, 'quiet', true));
catch err
    threw = true;
    assert(~isempty(strfind(err.message, k0)) && ~isempty(strfind(err.message, variants{1, 2})), ...
           'the refusal must print both keys');
end
assert(threw, 'a mismatched design key was NOT refused');
lastwarn('');
ws = warning('off', 'ea_apply_protocol:design_key');
pf2 = ea_apply_protocol({'block', est}, nl, vn, h_early, held, '', ...
                        struct('design_key', variants{1, 2}, 'force', true, 'quiet', true));
warning(ws);
[~, wid] = lastwarn();
assert(strcmp(wid, 'ea_apply_protocol:design_key'), 'force = true must warn');
assert(strcmp(pf2.key_status, 'mismatch-forced'));
% legacy null (no key) and flat layout (no substructs)
nl_flat = rmfield(nl, {'block', 'pooled', 'design_key'});
pf3 = ea_apply_protocol({'block', est; 'pooled', estp}, nl_flat, vn, h_early, held, '', ...
                        struct('design_key', k0, 'quiet', true));
assert(strcmp(pf3.key_status, 'legacy') && isequal(pf3.names, {'block'}), ...
       'flat legacy null: status %s, names %s', pf3.key_status, strjoin(pf3.names, ','));
assert(pf3.block.n_escape_q95 == 2, 'flat layout gives different counts');
% a null without draws: percentiles NaN, counts still there
nl_nod = nl;  nl_nod.block = rmfield(nl_nod.block, 'tau_bar_draws');
pf4 = ea_apply_protocol({'block', est}, nl_nod, vn, h_early, held, '', struct('quiet', true));
assert(all(isnan(pf4.block.null_pct(:))) && isnan(pf4.block.n_escape_holm) && isnan(pf4.block.n_escape_fwer) && ...
       pf4.block.n_escape_q95 == 2 && strcmp(pf4.key_status, 'unchecked'), 'null without draws');
fprintf('  ea_apply_protocol: key refusal, forced warning, legacy/flat null, no-draw null: OK\n');

% --- 3. the cross-design table ------------------------------------------------
tdir = fullfile(tempdir, sprintf('ea_crossp_%d', round(rand * 1e6)));
mkdir(tdir);
blpb = est;  blpp = estp;                                          %#ok<NASGU>
varnames = vn;                                                     %#ok<NASGU>
f_res1 = fullfile(tdir, 'res1.mat');  f_nul1 = fullfile(tdir, 'nul1.mat');
save(f_res1, 'blpb', 'blpp', 'h_early', 'varnames', 'cfg');
save(f_nul1, '-struct', 'nl');
% second design: p = 3 in the run, but the null was keyed at p = 2
cfg2 = cfg;  cfg2.p = 3;
f_res2 = fullfile(tdir, 'res2.mat');  f_nul2 = fullfile(tdir, 'nul2.mat');
blpb = estp;  cfg_keep = cfg;  cfg = cfg2;                          %#ok<NASGU>
ds_info = struct('varnames', {vn}, 'T', 200);                      %#ok<NASGU>
save(f_res2, 'blpb', 'h_early', 'ds_info', 'cfg');
cfg = cfg_keep;
save(f_nul2, '-struct', 'nl');
spec = struct('label', {'fake p=2', 'fake p=3'}, ...
              'result_file', {f_res1, f_res2}, 'null_file', {f_nul1, f_nul2}, ...
              'held_blocks', {1, []});
f_csv = fullfile(tdir, 'cross.csv');
tab = ea_cross_p_table(spec, f_csv);
assert(numel(tab) == 3, 'expected 3 rows (2 estimators + 1), got %d', numel(tab));
assert(strcmp(tab(1).key_status, 'match') && strcmp(tab(3).key_status, 'mismatch-forced'), ...
       'key status per design: %s / %s', tab(1).key_status, tab(3).key_status);
assert(tab(1).n_escape_q95 == 2 && tab(1).n_free_cells == 6 && tab(3).n_free_cells == 9, 'counts in the table');
assert(~isempty(strfind(tab(1).top1, '(')) && ~isempty(strfind(tab(1).top1, '<-')), ...
       'top-cell string: %s', tab(1).top1);
% the top cells are the two planted escapes, never the held block, never NaN
tops = {tab(1).top1, tab(1).top2, tab(1).top3};
for t = 1:3
    assert(isempty(strfind(tops{t}, 'NaN')), 'top%d carries a NaN ratio: %s', t, tops{t});
    assert(isempty(strfind(tops{t}, ['<-' vn{held} '('])), 'top%d names the held block: %s', t, tops{t});
end
assert(~isempty(strfind(tops{1}, sprintf('(%.2f)', tab(1).max_ratio_q95))), ...
       'top1 (%s) should carry max_ratio_q95 = %.2f', tops{1}, tab(1).max_ratio_q95);
assert(any(strncmp(tops(1:2), 'a<-b(', 5)) && any(strncmp(tops(1:2), 'b<-c(', 5)), ...
       'top-2 should be the planted escapes, got %s | %s', tops{1}, tops{2});
assert(tab(1).max_ratio_q95 >= tab(1).median_ratio_q95, 'max < median');
assert(tab(3).p == 3 && tab(1).p == 2, 'p column');
ct = fileread(f_csv);
assert(strncmp(ct, '#', 1) && ~isempty(strfind(ct, 'not comparable across designs')), ...
       'csv header must state that raw tau is not compared');
rmdir(tdir, 's');
fprintf('  ea_cross_p_table: 2 designs, key status per row, csv: OK\n');

fprintf('PASS: test_null_modularity\n\n');
end

% -------------------------------------------------------------------------
function S = fake_null_stats(tau_all, pg_all, cover_all, h1, h2)
% Replicates null_stats of run_null_calibration (same field names).
[R, K, G, ~] = size(tau_all);
tb = mean(tau_all(:, :, :, h1:h2), 4);
pb = mean(pg_all(:, :, :, h1:h2), 4);
S.q90 = zeros(K, G);  S.q95 = zeros(K, G);  S.q99 = zeros(K, G);
S.pgt1_q95 = zeros(K, G);
S.tau_bar_mean = reshape(mean(tb, 1), [K, G]);
S.pgt1_mean    = reshape(mean(pb, 1), [K, G]);
S.argmax_freq  = zeros(K, G);
S.maxstat_q95  = zeros(K, 1);
S.pflag_q95    = zeros(K, 1);
for i = 1:K
    mx = reshape(max(tb(:, i, :), [], 3), [R, 1]);
    S.maxstat_q95(i) = empirical_quantile(mx, 0.95);
    px = reshape(max(pb(:, i, :), [], 3), [R, 1]);
    S.pflag_q95(i) = empirical_quantile(px, 0.95);
    for g = 1:G
        q3 = empirical_quantile(reshape(tb(:, i, g), [R, 1]), [0.90, 0.95, 0.99]);
        S.q90(i, g) = q3(1);  S.q95(i, g) = q3(2);  S.q99(i, g) = q3(3);
        S.pgt1_q95(i, g) = empirical_quantile(reshape(pb(:, i, g), [R, 1]), 0.95);
    end
    [~, am] = max(reshape(tb(:, i, :), [R, G]), [], 2);
    for g = 1:G
        S.argmax_freq(i, g) = mean(am == g);
    end
end
S.coverage_under_null = reshape(mean(mean(cover_all, 1), 3), [K, 1]);
S.tau_bar_draws = tb;
S.pgt1_draws = pb;
end

function hx = fnv_of(s)
% Independent 32-bit FNV-1a, checked against the published test vectors
% and then against ea_design_key's own hash of its canonical string.
b = double(uint8(s));
h = 2166136261;  prime = 16777619;
for k = 1:numel(b)
    h = bitxor(h, b(k));
    lo = mod(h, 65536);  hi = (h - lo) / 65536;
    h = mod(mod(hi * prime, 65536) * 65536 + lo * prime, 4294967296);
end
hx = lower(dec2hex(h, 8));
end
