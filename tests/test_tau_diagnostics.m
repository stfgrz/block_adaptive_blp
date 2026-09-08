function test_tau_diagnostics()
% PURPOSE
% -------
% Check the tau diagnostic accounting in montecarlo/tau_diagnostics.m.
% These numbers are the evidence for the "specification diagnostic"
% reading of the thesis, so the arithmetic behind them is pinned rather
% than eyeballed.
%
% WHAT IS CHECKED
% ---------------
% 1. Every documented field is present, with the documented shapes, for
%    BOTH scale estimators ('block' and 'pooled').
% 2. LOCALISATION arithmetic recomputed by hand on the stored tau array:
%    the aggregate and per-equation argmax frequencies, and the detection
%    probability, must equal what the module reports.
% 3. FREQUENCY IDENTITIES: each argmax-frequency row sums to 1, and the
%    detection probability equals the entry of the frequency vector at
%    the truly misspecified block.
% 4. RANK SANITY: the average rank of the true block lies in [1, G], and
%    .rank_chance is (G+1)/2 -- the value pure noise would deliver, which
%    is what a reported rank must be compared against.
% 5. FALSE-POSITIVE SETTING: on the CORRECT DGP the struct is stamped
%    .false_positive = true, detection probabilities are NaN (there is no
%    true block to detect) and the flag rates are, by construction, false
%    positive rates.  This is the guard against quoting a flag rate from
%    a misspecified DGP as if it were a false-positive rate.
% 6. P(tau > 1) is a probability in [0, 1] and is populated by both
%    estimators.
% 7. LEGACY RESULTS: a struct carrying only the old top-level
%    mc.tau_mean (results saved before mc.tau existed) is still handled,
%    with the probability-based fields returned as NaN rather than
%    silently zero.

fprintf('test_tau_diagnostics:\n');

cfg = mc_preset('smoke');
cfg.mc.n_rep = 6;  cfg.H = 5;

mc_s = run_montecarlo(cfg, 'sparse');
mc_c = run_montecarlo(cfg, 'correct');

for which = {'block', 'pooled'}
    d = tau_diagnostics(mc_s, which{1});
    need = {'tau_bar', 'tau_med_bar', 'p_gt1_bar', 'tau_rep_q', ...
            'argmax_freq', 'argmax_freq_by_eq', 'argmax_freq_by_eq_early', ...
            'detect_prob', 'detect_prob_by_eq', 'detect_prob_by_eq_early', ...
            'mean_rank_true', 'mean_rank_true_by_eq', 'rank_chance', ...
            'flag_ratio', 'flag_prob', 'joint_detect_ratio', ...
            'joint_detect_prob', 'flag_rate', 'false_positive', ...
            'max_argmax_freq', 'max_argmax_freq_early', 'argmax_cell', ...
            'argmax_cell_early', 'argmax_chance'};
    for k = 1:numel(need)
        assert(isfield(d, need{k}), '%s diagnostics lack .%s', which{1}, need{k});
    end
    assert(isequal(size(d.tau_bar), [d.K, d.G, d.H]), ...
        '%s tau_bar has size %s', which{1}, mat2str(size(d.tau_bar)));
    assert(all(d.p_gt1_bar(:) >= 0 & d.p_gt1_bar(:) <= 1), ...
        '%s P(tau>1) outside [0,1]', which{1});
    assert(any(d.p_gt1_bar(:) > 0), '%s P(tau>1) is identically zero', which{1});
    assert(abs(sum(d.argmax_freq) - 1) < 1e-12, ...
        '%s argmax frequencies do not sum to 1', which{1});
    for i = 1:d.K
        assert(abs(sum(d.argmax_freq_by_eq(i, :)) - 1) < 1e-12, ...
            '%s per-equation frequencies do not sum to 1 in eq %d', which{1}, i);
    end
    assert(abs(d.detect_prob - d.argmax_freq(d.misspec_block)) < 1e-12, ...
        '%s detect_prob does not match the argmax frequency', which{1});
    assert(d.mean_rank_true >= 1 && d.mean_rank_true <= d.G, ...
        '%s mean rank %g outside [1, %d]', which{1}, d.mean_rank_true, d.G);
    assert(abs(d.rank_chance - (d.G + 1) / 2) < 1e-12, ...
        '%s rank_chance is not (G+1)/2', which{1});
    assert(~d.false_positive, ...
        '%s: the sparse DGP must not be flagged as a false-positive setting', which{1});
end
fprintf('  fields, shapes and frequency identities (block and pooled): OK\n');

% --- 2. recompute the localisation by hand -----------------------------
d = tau_diagnostics(mc_s, 'block');
tm = mc_s.tau.block.mean;                 % (K x G x H x R)
[K, G, H, R] = size(tm);
tau_rep = zeros(G, R);
for r = 1:R
    tau_rep(:, r) = squeeze(mean(mean(tm(:, :, :, r), 3), 1));
end
[~, win] = max(tau_rep, [], 1);
for g = 1:G
    assert(abs(d.argmax_freq(g) - mean(win == g)) < 1e-12, ...
        'aggregate argmax frequency for block %d recomputed differently', g);
end
for i = 1:K
    w = zeros(1, R);
    for r = 1:R
        [~, w(r)] = max(squeeze(mean(tm(i, :, :, r), 3)));
    end
    for g = 1:G
        assert(abs(d.argmax_freq_by_eq(i, g) - mean(w == g)) < 1e-12, ...
            'per-equation frequency (eq %d, block %d) recomputed differently', i, g);
    end
end
fprintf('  localisation frequencies recomputed by hand from mc.tau: OK\n');

% --- concentration statistic: definition and calibrated null -----------
% This is the statistic that discriminates (the aggregate flag rules do
% not, because they average tau over equations before comparing blocks).
% Its null is calibrated: under the correct DGP every cell of the
% equation x block winner map has expectation 1/G.
assert(abs(d.argmax_chance - 1 / d.G) < 1e-12, 'argmax_chance is not 1/G');
assert(abs(d.max_argmax_freq - max(d.argmax_freq_by_eq(:))) < 1e-12, ...
    'max_argmax_freq is not the max of the winner map');
assert(abs(d.max_argmax_freq_early - max(d.argmax_freq_by_eq_early(:))) < 1e-12, ...
    'max_argmax_freq_early is not the max of the early winner map');
assert(d.argmax_freq_by_eq(d.argmax_cell(1), d.argmax_cell(2)) == d.max_argmax_freq, ...
    'argmax_cell does not point at the maximum');
assert(d.argmax_freq_by_eq_early(d.argmax_cell_early(1), d.argmax_cell_early(2)) ...
       == d.max_argmax_freq_early, 'argmax_cell_early does not point at the maximum');
assert(d.max_argmax_freq >= d.argmax_chance - 1e-12, ...
    'the maximum of the winner map is below chance, which is impossible');
fprintf('  concentration statistic and its 1/G null: OK (max %.2f at eq %d, block %d)\n', ...
        d.max_argmax_freq, d.argmax_cell(1), d.argmax_cell(2));

% --- 5. the correct DGP is the false-positive setting ------------------
dc = tau_diagnostics(mc_c, 'block');
assert(dc.false_positive, 'the correct DGP must be a false-positive setting');
assert(isnan(dc.detect_prob), ...
    'there is no block to detect under the correct DGP: detect_prob must be NaN');
assert(all(isnan(dc.detect_prob_by_eq)), ...
    'per-equation detection must be NaN under the correct DGP');
assert(all(dc.flag_ratio >= 0 & dc.flag_ratio <= 1), 'flag rates outside [0,1]');
fprintf(['  correct DGP: flagged as a false-positive setting, detection NaN, ' ...
         'flag rates %.2f/%.2f/%.2f at thresholds %.2f/%.2f/%.2f: OK\n'], ...
        dc.flag_ratio(1), dc.flag_ratio(2), dc.flag_ratio(3), ...
        dc.flag_thresholds(1), dc.flag_thresholds(2), dc.flag_thresholds(3));

% --- 6. "nothing wrong" vs "nothing localisable" -------------------------
% The dense DGP IS misspecified but has no single block to find.  Marking
% it a false-positive setting would report its flag rate as a
% false-positive rate, which it is not.
cfg_d = cfg;  cfg_d.mc.n_rep = 4;
mc_d = run_montecarlo(cfg_d, 'dense');
dd = tau_diagnostics(mc_d, 'block');
assert(~dd.false_positive, ...
    'the dense DGP is misspecified and must not be a false-positive setting');
assert(dd.no_unique_block, ...
    'the dense DGP has no unique block and must be marked as such');
assert(~dc.no_unique_block, ...
    'the correct DGP must not be marked "no unique block" -- nothing is wrong');
assert(~d.no_unique_block, ...
    'the sparse DGP has a unique block and must not be marked otherwise');
% ...and a fitted lag order that nests the truth IS a false-positive setting
cfg_p = cfg;  cfg_p.p = 3;  cfg_p.mc.n_rep = 4;
mc_p = run_montecarlo(cfg_p, 'sparse');
dp = tau_diagnostics(mc_p, 'block');
assert(dp.false_positive && dp.fitted_p_nests_truth, ...
    ['a fitted VAR(3) nests the sparse truth, so that run must be a ' ...
     'false-positive setting']);
assert(isnan(dp.detect_prob), ...
    'there is no block to detect when the fitted model nests the truth');
fprintf(['  dense = misspecified but not localisable; fitted p >= 3 = nothing\n' ...
         '          wrong at all; the two are kept apart: OK\n']);

% --- 7. legacy result structs --------------------------------------------
legacy = rmfield(mc_s, 'tau');
dl = tau_diagnostics(legacy, 'block');
assert(all(isnan(dl.p_gt1_bar(:))), ...
    'a legacy result has no stored P(tau>1); it must come back NaN, not 0');
assert(all(isnan(dl.flag_prob)), ...
    'probability-based flags must be NaN when the input has no P(tau>1)');
assert(max(abs(dl.argmax_freq - d.argmax_freq)) < 1e-12, ...
    'the legacy path changed the localisation numbers');
fprintf('  legacy results (mc.tau_mean only) handled, probabilities NaN: OK\n');

fprintf('PASS: test_tau_diagnostics\n\n');
end
